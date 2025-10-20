#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=0-00:59:00
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=matlab
#SBATCH --error=/ix1/pmayo/outfiles/out_%A.out
#SBATCH --output=/ix1/pmayo/outfiles/out_%A.out
#SBATCH --mail-type=done,fail
#SBATCH --mail-user=knoneman@pitt.edu

# ----- Load environment -----
module purge
module load matlab/R2023a
module load python/ondemand-jupyter-python3.11

ENV_PATH=$(python -c "import config; print(config.ENV_PATH)")
source activate "$ENV_PATH"

RAW_DATA_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
OUT_DATA_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
NEV_PATH=$(python -c "import config; print(config.NEVUTIL_PATH)")

# ----- Specify inputs -----
echo "======================================================"

SESSION="${1}"
PROTOCOL="${2:-np-nodrift-ks4_wr12}"

echo "SESSION    =  '$SESSION'"
echo "PROTOCOL   =  $PROTOCOL"

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

#################################################################
################### RUN MATLAB SCRIPT #######################

RAW_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
OUT_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
NEV_PATH=$(python -c "import config; print(config.NEVUTIL_PATH)")

#echo "Running matlab pipeline........................"
#matlab -nodisplay <<EOF
#addpath(genpath('matlab'));   % add the subdirectory to path
#fprintf('Running process_fullRecording for $1\n');
#process_fullRecording('${SESSION}', ...
#    'RAW_DATA_PATH', '$RAW_PATH', ...
#    'OUT_DATA_PATH', '$OUT_PATH', ...
#    'NEVUTIL_PATH', '$NEV_PATH', ...
#    'SORTER_PATH', '$SORTER_PATH');
#exit
#EOF

#################################################################
HELPERS_PATH=$(python -c "import config; print(config.HELPERS_PATH)")

TBL_PATH="${OUT_DATA_PATH}/${SESSION}/tables/${SESSION}-${SORTER_PATH}.mat"
echo "TBL_PATH    =  '$TBL_PATH'"
echo "======================================================"

echo "Running extra matlab fxns........................"
matlab -nodisplay <<EOF
addpath(genpath('matlab'));
addpath(genpath('${HELPERS_PATH}/behavior'));
addpath(genpath('${HELPERS_PATH}/utils'));
addpath(genpath('${HELPERS_PATH}/neurons'));
fprintf('Running addToTbl_KKN for $1\n');
addToTbl_KKN('${TBL_PATH}', ...
    'SAVE_NAME', '${PROTOCOL}');
exit
EOF

#################################################################
echo "DONE"

crc-job-stats
