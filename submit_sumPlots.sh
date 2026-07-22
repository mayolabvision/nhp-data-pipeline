#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=plots
#SBATCH --error=/ix1/pmayo/outfiles/out_%A.out 
#SBATCH --output=/ix1/pmayo/outfiles/out_%A.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-02:59:59

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "My SLURM_ARRAY_TASK_COUNT is $SLURM_ARRAY_TASK_COUNT"
echo "Job started at $(date)"

module purge 
module load matlab/R2023a
module load python/ondemand-jupyter-python3.11

ENV_PATH=$(python -c "import config; print(config.ENV_PATH)")
source activate "$ENV_PATH"

RAW_DATA_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")

# ----- Specify inputs -----
echo "======================================================"

SESSION="${1}"
PROTOCOL="${3:-np-ks4}"

echo "SESSION    =  $SESSION"
echo "PROTOCOL   =  $PROTOCOL"

echo "======================================================"

#################################################################
################### RUN SI to get PROFILE #######################

echo "Retrieving protocol path........................"
RESULT=$($CONDA_PREFIX/bin/python -c "
from main_pipeline import profile_to_mat

sorter_path, short_path = profile_to_mat('${SESSION}', protocol='${PROTOCOL}.json')
print(sorter_path, short_path)
")

SORTER_HASH=$(echo "$RESULT" | awk '{print $1}')
SHORT_HASH=$(echo "$RESULT" | awk '{print $2}')

echo "SORTER_HASH = '$SORTER_HASH'"
echo "======================================================"

DATA_PATH="${RAW_DATA_PATH}/${SESSION}/tables/${SESSION}-${SHORT_HASH}.mat"
FIG_PATH="${RAW_DATA_PATH}/${SESSION}/figs/${SHORT_HASH}/"

echo "DATA_PATH: $DATA_PATH"
echo "FIG_PATH: $FIG_PATH"
echo "======================================================"

########################

echo "Plotting sorting QC........................"
matlab -nodisplay <<EOF
addpath(genpath('matlab'));
fprintf('Running plot_sortingQC for $1\n');
ia_mdirRasters('$DATA_PATH', '$FIG_PATH');
exit
EOF

####################################################################
echo "DONE"
crc-job-stats

