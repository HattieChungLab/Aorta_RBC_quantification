#!/bin/bash
# ONE job: step 2 for all slides (parallel on CPUs), then all 26 method variants in every lesion (GPU).
# Resumable: resubmit the same script if it stops; finished lesions are skipped and step 2 is not redone.
#
#   sbatch "/nfs/roberts/project/pi_hc859/swf23/TAA_histology/Mina_RBC_segmentation/RBC_segmentation/scripts/submit_all_methods.sh"
#
#SBATCH --job-name=rbc_all_methods_rerun
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=2:00:00
#SBATCH --output=/nfs/roberts/scratch/pi_hc859/swf23/aorta_histology/data/Mina_RBC_segmentation/logs/all_methods_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sam.fenske@yale.edu

# Slides to run: leave empty for all slides, or list slide IDs (codes, e.g. "S03 S07") separated by
# spaces to re-run just those. Their step 2 is always recomputed. Move their old outputs/all_methods/<slide>/
# folder away first, or the lesions already in it are skipped.
SLIDES=""

CONFIG="/nfs/roberts/project/pi_hc859/swf23/TAA_histology/Mina_RBC_segmentation/RBC_segmentation/config/base.yaml"
PY=/nfs/roberts/scratch/pi_hc859/swf23/rbc_env/bin/python
HERE=/nfs/roberts/project/pi_hc859/swf23/TAA_histology/Mina_RBC_segmentation/RBC_segmentation

cd /tmp
nvidia-smi -L
"$PY" "$HERE/scripts/run_all_methods.py" --config "$CONFIG" ${SLIDES:+--slides $SLIDES} || exit 1
"$PY" "$HERE/analysis/report_all_methods.py" --config "$CONFIG"
