#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=widgets
#SBATCH --error=/ix1/pmayo/outfiles/out_%A_%a.out 
#SBATCH --output=/ix1/pmayo/outfiles/out_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-02:59:59
#SBATCH --array=0-1

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "My SLURM_ARRAY_TASK_COUNT is $SLURM_ARRAY_TASK_COUNT"
echo "Job started at $(date)"

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

echo "SESSION    =  '$SESSION'"
echo "PROBE_ID   =  $PROBE_ID"
echo "PROTOCOL   =  $PROTOCOL"
echo "ENV_PATH   =  '$ENV_PATH'"

echo "======================================================"

#################################################################
####################### RUN SI Sorting ##########################

echo "Running postprocessing pipeline........................"
$CONDA_PREFIX/bin/python -c "
from main_pipeline import run_widgets

run_widgets(
    '${SESSION}', 
    probe_id=int('$PROBE_ID'),
    protocol='${PROTOCOL}.json',
    job_id=int(0),
    n_chunks=int(1)
)"

echo "======================================================"

#################################################################
echo "DONE"

crc-job-stats

