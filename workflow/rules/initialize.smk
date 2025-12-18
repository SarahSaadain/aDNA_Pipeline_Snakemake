# =================================================================================================
#     Dependencies and Environment Setup
# =================================================================================================
# Import required Python modules for system, platform, logging, and workflow management
import os
import sys  
import pwd
import re
import socket, platform
import subprocess
from datetime import datetime
import logging
from pathlib import Path
import yaml
import json

# Import Snakemake plugin settings for executor modes
from snakemake_interface_executor_plugins.settings import ExecMode

# --- Logging Setup (EARLY) ---
# Configure logging format and output for workflow debugging and status reporting
LOG_FORMAT = '[%(asctime)s] [%(levelname)s] %(message)s'
LOG_DATE_FORMAT = '%Y-%m-%d %H:%M:%S (%Z)'

logging.basicConfig(  # Basic config ASAP (for fallback)
    level=logging.INFO,
    format=LOG_FORMAT,
    datefmt=LOG_DATE_FORMAT,
    handlers=[logging.StreamHandler()]  # Only console for now
)

# Add a description of the workflow to the final report
report: os.path.join(workflow.basedir, "reports/workflow.rst")

# =================================================================================================
#     Snakemake Version Check
# =================================================================================================
# Ensure the minimum required Snakemake version is available for compatibility
snakemake.utils.min_version("9.9.0")
basedir = workflow.basedir

aDNA_Pipeline_version = "0.0.1" 

# =================================================================================================
#     Configuration Files and Reporting
# =================================================================================================
# Specify the main configuration file for the workflow
configfile: "config/config.yaml"

# Get abs paths of all config files
cfgfiles = []
for cfg in workflow.configfiles:
    cfgfiles.append(os.path.abspath(cfg))
cfgfiles = "\n                        ".join(cfgfiles)

# =================================================================================================
#     Platform and OS Information
# =================================================================================================
# Gather platform and OS version information for reproducibility and debugging
pltfrm = platform.platform() + "; " + platform.version()
try:
    # Not available in all versions, so we need to catch this
    ld = platform.linux_distribution()
    if len(ld):
        pltfrm += "; " + ld
    del ld
except:
    pass

try:
    # Mac OS version comes back as a nested tuple?!
    # Need to merge the tuples...
    def merge_osx_tuple(x, bases=(tuple, list)):
        for e in x:
            if type(e) in bases:
                for e in merge_osx_tuple(e, bases):
                    yield e
            else:
                yield e

    mv = " ".join(merge_osx_tuple(platform.mac_ver()))
    if not mv.isspace():
        pltfrm += "; " + mv
    del mv, merge_osx_tuple
except:
    pass

# =================================================================================================
#     User and Host Information
# =================================================================================================
# Get the current user and hostname for provenance tracking
username = pwd.getpwuid(os.getuid())[0]
hostname = socket.gethostname()
hostname = hostname + ("; " + platform.node() if platform.node() != socket.gethostname() else "")

# =================================================================================================
#     Git Version Information
# =================================================================================================
# Get the git hash, if available, to track the exact code version used
try:
    process = subprocess.Popen(
        ["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = process.communicate()
    out = out.decode("ascii")
    aDNA_Pipeline_git_hash = out.strip()
    if aDNA_Pipeline_git_hash:
        aDNA_Pipeline_version += "-" + aDNA_Pipeline_git_hash
    del process, out, err, aDNA_Pipeline_git_hash
except:
    pass

# =================================================================================================
#     Conda Environment Information
# =================================================================================================
# Get the conda version, if available, for environment reproducibility
try:
    process = subprocess.Popen(["conda", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = out.decode("ascii")
    conda_ver = out[out.startswith("conda") and len("conda") :].strip()
    del process, out, err
    if not conda_ver:
        conda_ver = "n/a"
except:
    conda_ver = "n/a"

# Get the conda env name, if available.
# See https://stackoverflow.com/a/42660674/4184258
conda_env = os.environ["CONDA_DEFAULT_ENV"] + " (" + os.environ["CONDA_PREFIX"] + ")"
if conda_env == " ()":
    conda_env = "n/a"

# Get nicely wrapped command line for reproducibility
cmdline = sys.argv[0]
for i in range(1, len(sys.argv)):
    cmdline += " " + sys.argv[i]

# =================================================================================================
#     Workflow Header Logging
# =================================================================================================
# Main aDNA Pipeline header, helping with debugging etc for user issues
logger.info("aDNA Pipeline " + aDNA_Pipeline_version + " run:")

logger.info("\tDate:               " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
logger.info("\tPlatform:           " + pltfrm)
logger.info("\tHost:               " + hostname)
logger.info("\tUser:               " + username)
logger.info("\tConda:              " + str(conda_ver))
logger.info("\tPython:             " + str(sys.version.split(" ")[0]))
logger.info("\tSnakemake:          " + str(snakemake.__version__))
logger.info("\tConda env:          " + str(conda_env))
logger.info("\tCommand:            " + cmdline)
logger.info("\tBase directory:     " + workflow.basedir)
logger.info("\tWorking directory:  " + os.getcwd())
logger.info("\tConfig file(s):     " + cfgfiles)

# Print config values of pipline 
# Convert config dict to a YAML-style string for clean formatting
config_str = yaml.dump(config["pipeline"], sort_keys=False, default_flow_style=False)

# Log it
logging.info("Loaded configuration:\n%s", config_str)

# Extract and log species with names and keys
species_section = config.get("species", {})
species_list = [
    f"{sdata.get('name', sname)} [{sname}]" for sname, sdata in species_section.items()
]

logging.info("Detected species:\n- %s", "\n- ".join(species_list))


# Clean up variables that are no longer needed to avoid polluting the namespace
# No need to have these output vars available in the rest of the snakefiles
del pltfrm, hostname, username
del conda_ver, conda_env
del cmdline, cfgfiles
# =================================================================================================
# End of initialize.smk
# =================================================================================================