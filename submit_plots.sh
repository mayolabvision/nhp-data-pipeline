#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=plots
#SBATCH --error=/ix1/pmayo/outfiles/out_%A_%a.out 
#SBATCH --output=/ix1/pmayo/outfiles/out_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-02:59:59
#SBATCH --array=0-49

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
HELPERS_PATH=$(python -c "import config; print(config.HELPERS_PATH)")

# ----- Specify inputs -----
echo "======================================================"

SESSION="${1}"
PROTOCOL="${2:-np-nodrift-ks4_wr12}"
PROBE_ID="${3:-1}"

echo "SESSION    =  $SESSION"
echo "PROTOCOL   =  $PROTOCOL"
echo "PROBE_ID   =  $PROBE_ID"

echo "======================================================"

#################################################################
################### RUN SI to get PROFILE #######################

echo "Retrieving protocol path........................"
SORTER_PATH=$($CONDA_PREFIX/bin/python -c "
from main_pipeline import profile_to_mat

print(profile_to_mat(
    '${SESSION}', 
    protocol='${PROTOCOL}.json' 
))
")

echo "SORTER_PATH    =  '$SORTER_PATH'"
echo "======================================================"

#DATA_PATH="${RAW_DATA_PATH}/${SESSION}/tables/${SESSION}-${PROTOCOL}.mat"
DATA_PATH="${RAW_DATA_PATH}/${SESSION}/tables/${SESSION}-${SORTER_PATH}.mat"
FIG_PATH="${RAW_DATA_PATH}/${SESSION}/figs/${SORTER_PATH}/"

echo "DATA_PATH: $DATA_PATH"
echo "FIG_PATH: $FIG_PATH"
echo "======================================================"

########################

echo "Plotting RF Maps per cluster........................"
matlab -nodisplay <<EOF
addpath(genpath('matlab'));
addpath(genpath('$HELPERS_PATH/plotting'));
fprintf('Running ia_rfMaps for $1\n');
ia_rfMaps('$DATA_PATH', ...
    'PROBE_INDEX', $PROBE_ID, ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')),...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

########################

echo "Plotting MDIR rasters per cluster........................"
ALIGN="stim"
echo "ALIGN: $ALIGN"
# First MATLAB call: run process_NeuropixRecording_KKN
matlab -nodisplay <<EOF
addpath(genpath('matlab'));
addpath(genpath('$HELPERS_PATH/plotting'));
fprintf('Running ia_mdirRasters for $1\n');
ia_mdirRasters('$DATA_PATH', ...
    'PROBE_INDEX', $PROBE_ID, ...
    'ALIGN', '$ALIGN', ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

ALIGN="sacc"
echo "ALIGN: $ALIGN"
# First MATLAB call: run process_NeuropixRecording_KKN
matlab -nodisplay <<EOF
addpath(genpath('matlab'));
addpath(genpath('$HELPERS_PATH/plotting'));
fprintf('Running ia_mdirRasters for $1\n');
ia_mdirRasters('$DATA_PATH', ...
    'PROBE_INDEX', $PROBE_ID, ...
    'ALIGN', '$ALIGN', ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

########################

echo "Plotting PURS rasters per cluster........................"
ALIGN="targ"
PURE_ONLY="0"
echo "ALL TRIALS, ALIGNED TO TARG"
matlab -nodisplay <<EOF
addpath(genpath('matlab'));
addpath(genpath('$HELPERS_PATH/plotting'));
addpath(genpath('$HELPERS_PATH/behavior'));
fprintf('Running ia_pursRasters for $1\n');
ia_pursRasters('$DATA_PATH', ...
    'PROBE_INDEX', $PROBE_ID, ...
    'ALIGN', '$ALIGN', ...
    'PURE_ONLY', logical(str2double('$PURE_ONLY')), ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

ALIGN="purs"
PURE_ONLY="0"
echo "ALL TRIALS, ALIGNED TO PURS"
matlab -nodisplay <<EOF
addpath(genpath('matlab'));
addpath(genpath('$HELPERS_PATH/plotting'));
addpath(genpath('$HELPERS_PATH/behavior'));
fprintf('Running ia_pursRasters for $1\n');
ia_pursRasters('$DATA_PATH', ...
    'PROBE_INDEX', $PROBE_ID, ...
    'ALIGN', '$ALIGN', ...
    'PURE_ONLY', logical(str2double('$PURE_ONLY')), ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

ALIGN="targ"
PURE_ONLY="1"
echo "PURE ONLY, ALIGNED TO TARG"
matlab -nodisplay <<EOF
addpath(genpath('matlab'));
addpath(genpath('$HELPERS_PATH/plotting'));
addpath(genpath('$HELPERS_PATH/behavior'));
fprintf('Running ia_pursRasters for $1\n');
ia_pursRasters('$DATA_PATH', ...
    'PROBE_INDEX', $PROBE_ID, ...
    'ALIGN', '$ALIGN', ...
    'PURE_ONLY', logical(str2double('$PURE_ONLY')), ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

ALIGN="purs"
PURE_ONLY="1"
echo "PURE ONLY, ALIGNED TO PURS"
matlab -nodisplay <<EOF
addpath(genpath('matlab'));
addpath(genpath('$HELPERS_PATH/plotting'));
addpath(genpath('$HELPERS_PATH/behavior'));
fprintf('Running ia_pursRasters for $1\n');
ia_pursRasters('$DATA_PATH', ...
    'PROBE_INDEX', $PROBE_ID, ...
    'ALIGN', '$ALIGN', ...
    'PURE_ONLY', logical(str2double('$PURE_ONLY')), ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

####################################################################
echo "DONE"
crc-job-stats

