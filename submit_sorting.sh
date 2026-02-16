#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=0-09:59:59
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=a100_nvlink
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=sorting
#SBATCH --error=/ix1/pmayo/outfiles/out_%A_%a.out
#SBATCH --output=/ix1/pmayo/outfiles/out_%A_%a.out
#SBATCH --mail-type=done,fail
#SBATCH --mail-user=mayolab@pitt.edu
#SBATCH --array=0-1

# ----- Load environment -----
module purge
module load python/ondemand-jupyter-python3.11

ENV_PATH=$(python -c "import config; print(config.ENV_PATH)")
source activate "$ENV_PATH"

# ----- Specify inputs -----
echo "======================================================"

SESSION="${1}"
PROTOCOL="${2:-np-ks4}"
PROBE_ID=$SLURM_ARRAY_TASK_ID

echo "SESSION    =  $SESSION"
echo "PROBE_ID   =  $PROBE_ID"
echo "PROTOCOL   =  $PROTOCOL"

echo "======================================================"

#################################################################
####################### RUN SI Sorting ##########################

echo "Running spike sorting pipeline........................"
$CONDA_PREFIX/bin/python -c "
from main_pipeline import run_sorting

run_sorting(
    '${SESSION}', 
    probe_id=int('$PROBE_ID'),
    protocol='${PROTOCOL}.json' 
)"

echo "======================================================"

#################################################################
echo "DONE"

crc-job-stats
