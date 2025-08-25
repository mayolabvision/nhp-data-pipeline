#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=0-23:59:59
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=pipe-postproc
#SBATCH --error=/ix1/pmayo/outfiles/out_%A_%a.out
#SBATCH --output=/ix1/pmayo/outfiles/out_%A_%a.out
#SBATCH --mail-type=done,fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --array=0-1

# ----- Load environment -----
module purge
module load python/ondemand-jupyter-python3.9

ENV_PATH=$(python -c "import config; print(config.ENV_PATH)")
source activate "$ENV_PATH"

# ----- Specify inputs -----
echo "======================================================"

SESSION="${1}"
PROBE_ID=$SLURM_ARRAY_TASK_ID
PROTOCOL="np_medicine.json"

echo "SESSION    =  '$SESSION'"
echo "PROBE_ID   =  $PROBE_ID"
echo "PROTOCOL   =  $PROTOCOL"
echo "ENV_PATH   =  '$ENV_PATH'"

echo "======================================================"

#################################################################
####################### RUN SI Sorting ##########################

echo "Running postprocessing pipeline........................"
$CONDA_PREFIX/bin/python -c "
from main_pipeline import run_postprocess

run_postprocess(
    '${SESSION}', 
    probe_id=int('$PROBE_ID'),
    protocol='${PROTOCOL}' 
)"

echo "======================================================"

#################################################################
echo "DONE"

crc-job-stats
