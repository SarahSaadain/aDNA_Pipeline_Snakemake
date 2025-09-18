# =================================================================================================
#     Dependencies
# =================================================================================================

import os, sys, pwd, re
import socket, platform
import subprocess
from datetime import datetime
import logging
from pathlib import Path

from snakemake_interface_executor_plugins.settings import ExecMode

# --- Logging Setup (EARLY) ---
LOG_FORMAT = '[%(asctime)s] [%(levelname)s] %(message)s'
LOG_DATE_FORMAT = '%Y-%m-%d %H:%M:%S (%Z)'

logging.basicConfig(  # Basic config ASAP (for fallback)
    level=logging.INFO,
    format=LOG_FORMAT,
    datefmt=LOG_DATE_FORMAT,
    handlers=[logging.StreamHandler()]  # Only console for now
)

# Ensure min Snakemake version
snakemake.utils.min_version("9.9.0")
basedir = workflow.basedir

aDNA_Pipeline_version = "0.0.1" 

configfile: "config.yaml"

# Add a description of the workflow to the final report
report: os.path.join(workflow.basedir, "reports/workflow.rst")

# Get the OS and version
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


# Get the user and hostname
username = pwd.getpwuid(os.getuid())[0]
hostname = socket.gethostname()
hostname = hostname + ("; " + platform.node() if platform.node() != socket.gethostname() else "")


# Get the git hash, if available. 
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

# Get the conda version, if available. 
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

# Get nicely wrapped command line
cmdline = sys.argv[0]
for i in range(1, len(sys.argv)):
    cmdline += " " + sys.argv[i]

# Get abs paths of all config files
cfgfiles = []
for cfg in workflow.configfiles:
    cfgfiles.append(os.path.abspath(cfg))
cfgfiles = "\n                        ".join(cfgfiles)

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

# No need to have these output vars available in the rest of the snakefiles
del pltfrm, hostname, username
del conda_ver, conda_env
del cmdline, cfgfiles