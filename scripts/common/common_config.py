import os
from datetime import datetime
import logging
import common.common_constants as common_constants
import common.common_logging as common_logging
import common.config_manager as config_manager 
from common.common_config_enumerations import ConfigSettings

from common.config_manager import ConfigManager

from enum import Enum



# Load the configuration file (only once)
try:
    common_logging.print_step_execution("Loading configuration file ...")
    # Path to your config file. we assume it is in the same directory where this script is located
    config_file_path = 'config.yaml'  

    common_logging.print_info(f"Loading config file: {config_file_path}")

    config_dict = config_manager.load_config(config_file_path)  # Or however you specify the path

    config = ConfigManager(config_dict)
    
    common_logging.print_info("Config loaded successfully.")
    
    common_logging.print_info(f"Project path: {config.get_project_path()}")
    common_logging.print_info(f"Species: {config.get_species_list()}")
    common_logging.print_info(f"Threads default: {config.get_threads()}")

    # Set up logging based on the config
    # Assuming the config has a 'log_level' key
    # and you want to set the logging level accordingly#
    # --- Update Log Level (if in config) ---
    log_level_config = config.get_log_level()
    log_level = getattr(logging, log_level_config, common_logging.print_info)

    # update the logging level
    logging.getLogger().setLevel(log_level)
    common_logging.print_info(f"Log level set to {log_level_config}")

    # --- Add File Handler ---
    log_dir = os.path.join(config.get_project_path(), common_constants.FOLDER_LOGS)
    os.makedirs(log_dir, exist_ok=True)

    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    log_filename = os.path.join(log_dir, f'{timestamp}_pipeline.log')

    common_logging.print_info(f"Log file: {log_filename}")

    file_handler = logging.FileHandler(log_filename)
    file_handler.setFormatter(logging.Formatter(common_logging.LOG_FORMAT, datefmt=common_logging.LOG_DATE_FORMAT))
    
    logging.getLogger().addHandler(file_handler)

    common_logging.print_debug(f"Loaded Config: {config_dict}") 

except FileNotFoundError:
    common_logging.print_error("Config file not found.  Exiting.")
    exit(1) #Or handle more gracefully


PATH_ADNA_PROJECT = config.get_project_path()

# config
THREADS_DEFAULT = config.get_threads()

# species folders
FOLDER_SPECIES = config.get_species_list()

# program related constants
PROGRAM_PATH_FASTP = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'fastp', default='fastp')
PROGRAM_PATH_MULTIQC = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'multiqc', default='multiqc')
PROGRAM_PATH_FASTQC = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'fastqc', default='fastqc')
PROGRAM_PATH_BEDTOOLS = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'bedtools', default='bedtools')
PROGRAM_PATH_BAMTOBED = "bamtobed"
PROGRAM_PATH_BWA = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'bwa', default='bwa')
PROGRAM_PATH_BWA_MEM = "mem"
PROGRAM_PATH_BWA_INDEX = "index"
PROGRAM_PATH_SAMTOOLS = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'samtools', default='samtools')
PROGRAM_PATH_SAMTOOLS_FAIDX =  "faidx"
PROGRAM_PATH_SAMTOOLS_VIEW =  "view"
PROGRAM_PATH_SAMTOOLS_SORT = "sort"
PROGRAM_PATH_SAMTOOLS_INDEX = "index"
PROGRAM_PATH_SAMTOOLS_DEPTH = "depth"
PROGRAM_PATH_ANGSD = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'angsd', default='angsd')
PROGRAM_PATH_SEQKIT = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'seqkit', default='seqkit')
PROGRAM_PATH_SEQKIT_STATS = "stats"
PROGRAM_PATH_CENTRIFUGE = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'centrifuge', default='centrifuge')
PROGRAM_PATH_KRAKEN = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'kraken', default='kraken2')
PROGRAM_PATH_MAPDAMAGE = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'mapdamage', default='mapDamage')
PROGRAM_PATH_ECMSD = config.get_value(ConfigSettings.TOOLS_CONFIG.value, 'ecmsd', default='ECMSD')
