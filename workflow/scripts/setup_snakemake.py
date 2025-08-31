# scripts/setup_snakemake.py
import os
import sys

if sys.platform == "darwin":  # Only true on macOS
    # On macOS, conda will prefer osx-64 packages (even on Apple Silicon, 
    # it will pull x86_64 binaries, which often avoids compilation headaches).
    # On Linux or Windows, the variable is never set, so nothing breaks.
    # Note: forcing osx-64 on Apple Silicon means you’ll run Intel binaries under 
    # Rosetta 2 emulation. That usually works fine, but can be slower than arm64-native 
    # builds. If speed matters, you might instead want to try fixing the compiler setup 
    # (so arm64 builds succeed).
    os.environ["CONDA_SUBDIR"] = "osx-64"

# Absolute path to project root
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
scripts_dir = os.path.join(project_root, "scripts")
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)