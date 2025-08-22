# pipeline_steps_enum.py

from enum import Enum

class ConfigSettings(Enum):
    PIPELINE_CONFIG = 'pipeline'
    SPECIES_CONFIG = 'species'
    TOOLS_CONFIG = 'tools'

# Define Enums for pipeline stages
class PipelineStages(Enum):
    RAW_READS_PROCESSING = 'raw_reads_processing'
    REFERENCE_GENOME_PROCESSING = 'reference_genome_processing' # Corrected name
    POST_PROCESSING = 'post_processing'

# Define Enums for pipeline steps within each stage
class RawReadsProcessingSteps(Enum):
    QUALITY_CHECK_RAW = 'quality_checking_raw'
    QUALITY_CHECK_ADAPTER_REMOVED = 'quality_checking_adapter_removed'
    QUALITY_CHECK_QUALITY_FILTERED = 'quality_checking_quality_filtered'
    QUALITY_CHECK_DUPLICATES_REMOVED = 'quality_checking_duplicates_removed'
    ADAPTER_REMOVAL = 'adapter_removal' # Corrected name
    QUALITY_FILTERING = 'quality_filtering' # Corrected name
    DEDUPLICATION = 'deduplication' # Corrected name
    MERGE_READS_BY_INDIVIDUAL = 'merge_reads_by_individual'
    GENERATE_QUALITY_CHECK_REPORT = 'generate_quality_check_report'
    READS_PROCESSING_ANALYSIS = 'reads_processing_analysis'
    READ_LENGTH_DISTRIBUTION_ANALYSIS = 'read_length_distribution_analysis'
    CONTAMINATION_ANALYSIS = 'contamination_analysis'
    GENERATE_RAW_READS_PLOTS = 'generate_raw_reads_plots'


class ReferenceGenomeProcessingSteps(Enum):
    PREPARE_REFERENCE_GENOME = 'prepare_reference_genome' # Corrected name
    MAP_READS_TO_REFERENCE_GENOME = 'map_reads_to_reference_genome' # Corrected name
    ENDOGENOUS_READS_ANALYSIS = 'endogenous_reads_analysis'
    CREATE_CONSENSUS_SEQUENCE = 'create_consensus_sequence'
    COVERAGE_ANALYSIS = 'coverage_analysis'
    DAMAGE_ANALYSIS = 'damage_analysis'
    EXTRACT_SPECIAL_SEQUENCES = 'extract_special_sequences'
    GENERATE_REF_GENOME_PLOTS = 'generate_reference_genome_plots'


class PostProcessingSteps(Enum):
    # Added 'mtdna_analysis' key
    MTDNA_ANALYSIS = 'mtdna_analysis'
    GENERATE_SPECIES_COMPARISON_PLOTS = 'generate_species_comparison_plots'

# Define Enums for mtDNA analysis steps specifically under post_processing.mtdna_analysis
class MtdnaAnalysisSteps(Enum):
    MTDNA_MAP_TO_REF_GENOME = 'mtdna_map_to_ref_genome'
    MTDNA_DETERMINE_REGIONS = 'mtdna_determine_regions'
    MTDNA_CREATE_AND_MAP_CONSENSUS = 'mtdna_create_and_map_consensus'
    MTDNA_EXTRACT_COI_REGIONS = 'mtdna_extract_coi_regions'
    MTDNA_CHECK_EXTRACTED_REGIONS = 'mtdna_check_extracted_regions'

class AdapterRemovalSettings(Enum):
    ADAPTERS = 'adapters'
    ADAPTERS_R1 = 'r1'
    ADAPTERS_R2 = 'r2'

class ContaminationCheckSettings(Enum):
    CENTRIFUGE_DB = 'centrifuge_db'    
    KRAKEN_DB = 'kraken_db'

class QualityControlSettings(Enum):
    THREADS = 'threads'

# You might also want a mapping from stage key strings to their step Enums
STAGE_STEP_ENUM_MAP = {
    PipelineStages.RAW_READS_PROCESSING.value: RawReadsProcessingSteps,
    PipelineStages.REFERENCE_GENOME_PROCESSING.value: ReferenceGenomeProcessingSteps,
    PipelineStages.POST_PROCESSING.value: PostProcessingSteps,
}