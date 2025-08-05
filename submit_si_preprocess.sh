#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name=si
#SBATCH --error=/ix1/pmayo/outfiles/spikesort/out_%A.out
#SBATCH --output=/ix1/pmayo/outfiles/spikesort/out_%A.out
#SBATCH --mail-type=done,fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-08:00:00

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

# INPUTS = SESSION_NAME, PROBE_ID
module purge
module load python/ondemand-jupyter-python3.10
source activate /ihome/pmayo/knoneman/.conda/envs/si_env

SESSION_NAME="${1}"
PROBE_ID=${2}

DRIFT_CORRECT="dredge"
REMOVE_BAD_CHANNELS=true
RAW_DATA_PATH="/ix1/pmayo/lab_NHPdata/"

echo "SESSION_NAME = $SESSION_NAME"
echo "PROBE_ID = $PROBE_ID"
echo "DRIFT_CORRECT = $DRIFT_CORRECT"
echo "REMOVE_BAD_CHANNELS = $REMOVE_BAD_CHANNELS"
echo "RAW_DATA_PATH = $RAW_DATA_PATH"

# Run the Python function with environment variable conversion
$CONDA_PREFIX/bin/python -c "
from run_si_pipeline import run_si_preprocess

def str2bool(val):
    return str(val).lower() in ('true', '1', 'yes')

run_si_preprocess(
    '$SESSION_NAME',
    probe_id=int('$PROBE_ID'),
    drift_correct='$DRIFT_CORRECT',
    remove_bad_channels=str2bool('$REMOVE_BAD_CHANNELS'),
    raw_data_path='$RAW_DATA_PATH'
)
"

echo "DONE"

crc-job-stats
