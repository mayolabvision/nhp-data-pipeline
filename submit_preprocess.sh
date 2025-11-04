#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=0-16:00:00
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=NHPipe-pp
#SBATCH --error=/ix1/pmayo/outfiles/out_%A_%a.out
#SBATCH --output=/ix1/pmayo/outfiles/out_%A_%a.out
#SBATCH --mail-type=done,fail
#SBATCH --mail-user=knoneman@pitt.edu
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
################ Check that raw signal exists ###################
#module load matlab/R2023a
#NEV_PATH=$(python -c "import config; print(config.NEVUTIL_PATH)")
#
#RAW_DATA_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
#PROBES_PATH=$(python -c "import config; print(config.PROBES_PATH)")
#DATA_PATH="${RAW_DATA_PATH}/${SESSION}"
#
#echo "Checking raw signal exists........................"
#matlab -nodisplay <<EOF
#addpath(genpath('matlab'));
#addpath(genpath('${NEV_PATH}'));
#fprintf('Running rawRipple_to_binaryFile for $1\n');
#rawRipple_to_binaryFile('${DATA_PATH}', '${PROBES_PATH}');
#exit
#EOF
#
#echo "======================================================"

#################################################################
##################### RUN SI (PP + Sort) ########################

echo "Running preprocessing pipeline........................"
$CONDA_PREFIX/bin/python -c "
from main_pipeline import run_preprocess

run_preprocess(
    '${SESSION}', 
    probe_id=int('${PROBE_ID}'),
    protocol='${PROTOCOL}.json' 
)"

echo "======================================================"

#################################################################
echo "DONE"
crc-job-stats
