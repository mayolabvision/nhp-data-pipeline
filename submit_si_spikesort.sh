#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=0-08:00:00
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=si
#SBATCH --error=/ix1/pmayo/outfiles/si/out_%A_%a.out
#SBATCH --output=/ix1/pmayo/outfiles/si/out_%A_%a.out
#SBATCH --mail-type=done,fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --array=0-1

#export MKL_NUM_THREADS=1
#export NUMEXPR_NUM_THREADS=1
#export OMP_NUM_THREADS=1

# INPUTS = SESSION_NAME, PROBE_ID
module purge
module load python/ondemand-jupyter-python3.10
source activate /ix1/pmayo/envs/si_env 

################################################################

SESSION="${1}"
PROBE_ID=$SLURM_ARRAY_TASK_ID

echo "========================"
echo "SESSION    =  '$SESSION'"
echo "PROBE_ID   =  $PROBE_ID"

PROTOCOL_FILE="np_protocol.json"
RAW_DATA_PATH="/ix1/pmayo/lab_NHPdata"

echo "PROTOCOL   =  '$PROTOCOL_FILE'"
echo "========================"

#################################################################
##################### RUN SI (PP + Sort) ########################

echo "Running spike interface"
$CONDA_PREFIX/bin/python -c "
from run_si_pipeline import run_si_sorting

run_si_sorting(
    '${SESSION}', 
    probe_id=int('$PROBE_ID'), 
    protocol_file='${PROTOCOL_FILE}', 
    raw_data_path='${RAW_DATA_PATH}'
)"

echo "========================"

#################################################################
echo "DONE"

crc-job-stats
