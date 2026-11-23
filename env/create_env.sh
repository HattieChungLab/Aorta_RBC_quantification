#!/bin/bash
# Build the RBC pipeline environment (Python from conda, EVERYTHING else from pip).
# Pip-only native libraries avoid the conda/pip OpenMP clash that crashed HistoPLUS (handover §6).
#
#   bash "/nfs/roberts/project/pi_hc859/swf23/TAA_histology/Mina_RBC_segmentation/RBC_segmentation/env/create_env.sh"
#
# Needs internet (pip, GitHub). Exports the resolved env to rbc_env.yml next to this script.
set -euo pipefail

ENV=/nfs/roberts/scratch/pi_hc859/swf23/rbc_env
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# cu128 / torch 2.7: needed for Blackwell GPUs (sm_120); also runs on older GPUs (V100..H100)
TORCH_INDEX=https://download.pytorch.org/whl/cu128

module load miniconda
conda create -y -p "$ENV" python=3.11 pip
PIP="$ENV/bin/python -m pip"

$PIP install --upgrade pip
# torch + a prebuilt xformers wheel built against exactly this torch (no source build, no shim)
$PIP install torch==2.7.0 torchvision==0.22.0 xformers==0.0.30 --index-url "$TORCH_INDEX"
# HistoPLUS (not on PyPI) - pulls timm==1.0.8, openslide-python/-bin, shapely, scikit-image, opencv-headless
$PIP install "git+https://github.com/owkin/histoplus.git"
# Cellpose-SAM + blob-detector pipeline deps
$PIP install cellpose scipy pillow pyyaml pandas matplotlib
# extra segmenters for analysis/benchmark_segmenters.py; -c keeps torch/numpy at the versions above
$PIP freeze | grep -iE "^(torch|torchvision|numpy|xformers|timm|scikit-image)==" > /tmp/rbc_pins.txt
$PIP install instanseg-torch micro-sam -c /tmp/rbc_pins.txt
$PIP install "git+https://github.com/facebookresearch/dinov3" -c /tmp/rbc_pins.txt   # Cellpose cpdino models

# opencv: cellpose may pull the GUI build; keep only headless (one cv2 in the env)
if $PIP show opencv-python >/dev/null 2>&1; then
  $PIP uninstall -y opencv-python
  $PIP install --force-reinstall --no-deps opencv-python-headless
fi

conda env export -p "$ENV" > "$HERE/rbc_env.yml"
$PIP freeze > "$HERE/rbc_env_pip_freeze.txt"
echo "ENV READY: $ENV"
