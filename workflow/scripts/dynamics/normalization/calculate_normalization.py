import re
import collections
import pysam
import logging
import sys
from snakemake.script import snakemake
import numpy as np

################################################################################
# Constants
################################################################################

LOG_FORMAT = '[%(asctime)s] [%(levelname)s] %(message)s'
LOG_DATE_FORMAT = '%Y-%m-%d %H:%M:%S (%Z)'

################################################################################
# Global Variables
################################################################################

fle_and_scg_length_dictionary = {}
fle_and_scg_mapped_reads_dictionary = {}
scg_names_list = []
fle_names_list = []
unknown_names_list = []

################################################################################
# Functions
################################################################################

def classify_reference(name: str):
    if name.endswith("_fle"):
        return "fle"
    elif name.endswith("_scg"):
        return "scg"
    else:
        return "?"

def get_weight_from_read(read):
    cigar = read.cigarstring
    readlen = read.query_length
    mapq = read.mapping_quality
    nm = read.get_tag("NM") if read.has_tag("NM") else None
    nh = read.get_tag("NH") if read.has_tag("NH") else 1
    reference_category = classify_reference(read.reference_name)

    aligned = 0.0
    softclip = 0

    for length, op in re.findall(r"(\d+)([A-Z])", cigar):
        length = int(length)
        if op == "=":
            aligned += length
        elif op == "M":
            aligned += 0.5 * length
        elif op == "X":
            aligned += 0.1 * length
        elif op == "S":
            softclip += length

    if readlen == 0:
        return 0.0

    aligned_frac = aligned / readlen
    softclip_frac = softclip / readlen

    if softclip_frac > 0.30:
        return 0.0

    nm_penalty = 1.0 / (1.0 + nm) if nm is not None else 1.0

    if reference_category == "scg":
        if mapq is None:
            mapq_penalty = 1
        elif mapq >= 30:
            mapq_penalty = 1.0
        elif mapq >= 20:
            mapq_penalty = 0.5
        else:
            return 0.0

    if reference_category == "fle":
        if mapq is None:
            mapq_penalty = 0.8
        elif mapq >= 30:
            mapq_penalty = 1.0
        elif mapq >= 20:
            mapq_penalty = 0.5
        else:
            return 0.2

    nh_penalty = 1.0 / nh if nh > 0 else 1.0

    weight = aligned_frac * nm_penalty * mapq_penalty * nh_penalty

    return weight


################################################################################
# Main
################################################################################

bam_file_path = snakemake.input.bam
ignore_duplicates = snakemake.params.get("ignore_duplicates", True)

summary_filename = snakemake.output.summary
coverage_filename = snakemake.output.coverage
log_filename = snakemake.log[0]

logging.basicConfig(
    level=logging.INFO,
    format=LOG_FORMAT,
    datefmt=LOG_DATE_FORMAT,
    handlers=[
        logging.StreamHandler(sys.stderr),
        logging.FileHandler(log_filename)
    ],
    force=True
)

logger = logging.getLogger(__name__)

logger.info(f"Starting normalization")
logger.info(f"BAM file: {bam_file_path}")
logger.info(f"Summary output file: {summary_filename}")
logger.info(f"Coverage output file: {coverage_filename}")
logger.info(f"Ignore duplicates: {ignore_duplicates}")

fle_weights_col = collections.defaultdict(float)
scg_weights_col = collections.defaultdict(float)

count_reads = 0
count_mapped_reads = 0
count_reads_passed_filter = 0
count_reads_passed_fle = 0
count_reads_passed_scg = 0
weighted_sum_passed_reads = 0.0
sum_fle_weights = 0.0
sum_scg_weights = 0.0

logger.info("Opening BAM file and reading header...")

with pysam.AlignmentFile(bam_file_path, "rb") as bam_file_object:

    logger.info("Processing references in BAM header (lengths, names)...")
    for name, length in zip(bam_file_object.references, bam_file_object.lengths):
        fle_and_scg_length_dictionary[name] = length
        fle_and_scg_mapped_reads_dictionary[name] = 0

        category = classify_reference(name)
        if category == "fle":
            fle_names_list.append(name)
        elif category == "scg":
            scg_names_list.append(name)
        else:
            unknown_names_list.append(name)
            logger.warning(f"Unknown reference category for {name}; storing under '?'")

    logger.info("Processing read alignments...")
    for read in bam_file_object:
        count_reads += 1

        if read.is_unmapped:
            continue

        if read.is_duplicate and ignore_duplicates:
            continue

        count_mapped_reads += 1

        ref_sequence_name = bam_file_object.get_reference_name(read.reference_id)
        category = classify_reference(ref_sequence_name)

        if read.query_sequence is None:
            continue

        if read.cigarstring is None:
            continue

        count_reads_passed_filter += 1

        read_weight = get_weight_from_read(read)
        weighted_sum_passed_reads += read_weight
        fle_and_scg_mapped_reads_dictionary[ref_sequence_name] += 1

        if category == "fle":
            fle_weights_col[ref_sequence_name] += read_weight
            sum_fle_weights += read_weight
            count_reads_passed_fle += 1
        elif category == "scg":
            scg_weights_col[ref_sequence_name] += read_weight
            sum_scg_weights += read_weight
            count_reads_passed_scg += 1
        else:
            logger.warning(f"Unknown reference '{ref_sequence_name}' encountered; categorizing as '?'")

if not fle_weights_col:
    logger.error("No FLE reads found in BAM file")
    raise Exception("No FLE reads found in BAM file")

if not scg_weights_col:
    logger.error("No SCG reads found in BAM file")
    raise Exception("No SCG reads found in BAM file")

logger.info("Calculating median SCG coverage for normalization using all SCGs...")

scg_coverages = []
for scg in scg_names_list:
    length = fle_and_scg_length_dictionary[scg]
    weight = scg_weights_col[scg]
    cov = weight / float(length)
    scg_coverages.append(cov)

if len(scg_coverages) == 0:
    logger.error("No SCG coverages found. Cannot normalize.")
    raise Exception("No SCG coverages found. Cannot normalize.")

median_scg_coverage = float(np.median(scg_coverages))
mean_scg_coverage = float(np.mean(scg_coverages))

if median_scg_coverage == 0.0:
    logger.error("Median SCG coverage is zero, cannot normalize.")
    raise Exception("Median SCG coverage is zero, cannot normalize.")

logger.info("Writing summary statistics to file...")
with open(summary_filename, 'w') as summary_file:
    summary_file.write("{0}\t{1}\t{2}\n".format("parameter", "ignore_duplicates", ignore_duplicates))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "number_of_scgs_provided", len(scg_names_list)))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "number_of_scgs_used_for_normalization", len(scg_names_list)))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "all_reads", count_reads))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "mapped_reads", count_mapped_reads))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "reads_with_mapq", count_reads_passed_filter))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "weighted_reads_with_mapq", round(weighted_sum_passed_reads, 2)))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "mapping_to_fles", count_reads_passed_fle))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "mapping_to_scgs", count_reads_passed_scg))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "mapping_to_fle_weighted", round(sum_fle_weights, 2)))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "mapping_to_scg_weighted", round(sum_scg_weights, 2)))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "mean_scg_coverage", round(mean_scg_coverage, 4)))
    summary_file.write("{0}\t{1}\t{2}\n".format("summary", "median_scg_coverage", round(median_scg_coverage, 4)))

logger.info(f"Number of SCGs used for normalization: {len(scg_names_list)}")
logger.info(f"Number of reads: {count_reads}")
logger.info(f"Number of mapped reads: {count_mapped_reads}")
logger.info(f"Number of mapped reads passing filters: {count_reads_passed_filter}")
logger.info(f"Number of mapped reads to FLEs: {count_reads_passed_fle}")
logger.info(f"Number of mapped reads to SCGs: {count_reads_passed_scg}")
logger.info(f"Mean SCG coverage: {mean_scg_coverage}")
logger.info(f"Median SCG coverage: {median_scg_coverage}")
logger.info("Summary statistics written.")

logger.info("Writing normalized coverage results to file...")
with open(coverage_filename, 'w') as coverage_file:
    coverage_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
        "type", "name", "length", "weighted_mapped_reads", "coverage", "normalized_coverage", "mapped_reads"))

    for scg in scg_names_list:
        weight = scg_weights_col[scg]
        length = fle_and_scg_length_dictionary[scg]
        coverage = float(weight) / float(length)
        normalized_coverage = coverage / median_scg_coverage
        mapped_reads = fle_and_scg_mapped_reads_dictionary[scg]
        coverage_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
            "scg", scg, length, round(weight, 2), round(coverage, 2), round(normalized_coverage, 2), mapped_reads))

    for te in fle_names_list:
        length = fle_and_scg_length_dictionary[te]
        weight = fle_weights_col[te]
        coverage = float(weight) / float(length)
        normalized_coverage = coverage / median_scg_coverage
        mapped_reads = fle_and_scg_mapped_reads_dictionary[te]
        coverage_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
            "te", te, length, round(weight, 2), round(coverage, 2), round(normalized_coverage, 2), mapped_reads))

logger.info("Normalized coverage results written.")
logger.info("Script completed successfully.")