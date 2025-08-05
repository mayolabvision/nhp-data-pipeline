#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=0-08:00:00
#SBATCH --gres=gpu:2
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --constraint=40g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name=si
#SBATCH --error=/ix1/pmayo/outfiles/spikesort/out_%A.out
#SBATCH --output=/ix1/pmayo/outfiles/spikesort/out_%A.out
#SBATCH --mail-type=done,fail
#SBATCH --mail-user=knoneman@pitt.edu

# INPUTS = SESSION_NAME, PROBE_ID
module purge
module load python/ondemand-jupyter-python3.10
source activate /ihome/pmayo/knoneman/.conda/envs/si_env

#export OPENBLAS_NUM_THREADS=1

SESSION_NAME="${1}"
PROBE_ID=${2}

SORTER_TYPE="kilosort4"
DRIFT_CORRECT="dredge"
RAW_DATA_PATH="/ix1/pmayo/lab_NHPdata/"

echo "SESSION_NAME = $SESSION_NAME"
echo "PROBE_ID = $PROBE_ID"
echo "SORTER_TYPE = $SORTER_TYPE"
echo "DRIFT_CORRECT = $DRIFT_CORRECT"
echo "RAW_DATA_PATH = $RAW_DATA_PATH"

# Run the Python function with environment variable conversion
$CONDA_PREFIX/bin/python -c "
from run_si_pipeline import run_si_spikesort

run_si_spikesort(
    '$SESSION_NAME',
    probe_id=int('$PROBE_ID'),
    sorter_type='$SORTER_TYPE',
    drift_correct='$DRIFT_CORRECT',
    raw_data_path='$RAW_DATA_PATH'
)
"

echo "DONE"

crc-job-stats
