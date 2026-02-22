#!/usr/bin/env python
import argparse
import logging
import pysam
from modules import SequenceEntryBuilder, FileWriter, load_fasta


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = argparse.ArgumentParser(description="""           
summarize coverage for diverse features
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
    Sarah Saadain
""")
parser.add_argument('--infile', type=str, dest="infile", required=True, help="Input BAM or SAM file path")
parser.add_argument("--fasta", type=str, required=True, dest="fasta", default=None, help="the fasta file to which reads were mapped")
parser.add_argument("--mapqth", type=int, required=False, dest="min_mapq", default=5, help="mapping quality threshold; below ambiguous")
parser.add_argument("--mc-snp", type=int, required=False, dest="mcsnp", default=5, help="minimum count of SNPs")
parser.add_argument("--mf-snp", type=float, required=False, dest="mfsnp", default=0.1, help="minimum frequency of SNPs")
parser.add_argument("--mc-indel", type=int, required=False, dest="mcindel", default=3, help="minimum count of indels")
parser.add_argument("--mf-indel", type=float, required=False, dest="mfindel", default=0.01, help="minimum frequency of indels")
parser.add_argument("--outfile", type=str, required=False, dest="outfile", default=None, help="output file in so format; if none is provided output will be screen")
parser.add_argument("--log-level", type=str, required=False, dest="loglevel", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)")


args = parser.parse_args()
logging.getLogger().setLevel(args.loglevel)

#if no output file is provided, don't write log to screen, otherwise it will mess up the output
if args.outfile is None:
    logging.getLogger().setLevel("ERROR")

writer = FileWriter(args.outfile)

# load fasta from file into dict
reference_dict = load_fasta(args.fasta)

sequence_entry_builder=None

infile_path = args.infile

mode = 'rb' if infile_path.lower().endswith('.bam') else 'r'
logging.info(f"Processing file: {infile_path} with mode: {mode}")

samfile = pysam.AlignmentFile(infile_path, mode)

for read in samfile:

    if read.is_unmapped:
        continue
    if read.is_secondary:
        continue
    if read.is_supplementary:
        continue

    ref_name = read.reference_name
    alignment_start_pos = read.reference_start + 1
    mapping_quality = read.mapping_quality if read.mapping_quality is not None else 0
    cigar_string = read.cigarstring
    read_sequence = read.query_sequence.upper() if read.query_sequence is not None else ''
    
    if cigar_string is None or cigar_string == '*':
        continue
    
    if sequence_entry_builder is None:
        ref_sequence = reference_dict[ref_name]
        sequence_entry_builder = SequenceEntryBuilder(ref_sequence, ref_name, args.min_mapq)

    if ref_name != sequence_entry_builder.ref_sequence_name:
        seq_entry = sequence_entry_builder.to_SequenceEntry(
            snp_min_count=args.mcsnp, 
            snp_min_frequency=args.mfsnp, 
            indel_min_count=args.mcindel, 
            indel_min_frequency=args.mfindel)
        
        writer.write(str(seq_entry))
        
        ref_sequence = reference_dict[ref_name]
        sequence_entry_builder = SequenceEntryBuilder(ref_sequence, ref_name, args.min_mapq)

    sequence_entry_builder.add_read(alignment_start_pos, cigar_string, mapping_quality, read_sequence)

samfile.close()

# process the last one as well
seq_entry = None

if sequence_entry_builder is not None:
    seq_entry = sequence_entry_builder.to_SequenceEntry(args.mcsnp, args.mfsnp, args.mcindel, args.mfindel)
    writer.write(str(seq_entry))



 