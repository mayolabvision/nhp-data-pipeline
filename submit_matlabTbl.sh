#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=0-04:59:00
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=matlab
#SBATCH --error=/ix1/pmayo/outfiles/out_%A.out
#SBATCH --output=/ix1/pmayo/outfiles/out_%A.out
#SBATCH --mail-type=done,fail
#SBATCH --mail-user=mayolab@pitt.edu

# ----- Load environment -----
module purge
module load matlab/R2023a
module load python/ondemand-jupyter-python3.11

ENV_PATH=$(python -c "import config; print(config.ENV_PATH)")
source activate "$ENV_PATH"

# ----- Specify inputs -----
echo "======================================================"

SESSION="${1}"
PROTOCOL="${2:-np-ks4}"

echo "SESSION    =  '$SESSION'"
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

SORTER_PATH=$(echo "$RESULT" | awk '{print $1}')
SHORT_PATH=$(echo "$RESULT" | awk '{print $2}')

echo "SORTER_PATH = '$SORTER_PATH'"

#################################################################
################### RUN MATLAB SCRIPT #######################

RAW_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
OUT_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
NEV_PATH=$(python -c "import config; print(config.NEVUTIL_PATH)")
NET_PATH=$(python -c "import config; print(config.NASNET_PATH)")
HELP_PATH=$(python -c "import config; print(config.HELPERS_PATH)")

echo "Running matlab pipeline........................"
matlab -nodisplay <<EOF
addpath(genpath('matlab'));   % add the subdirectory to path
fprintf('Running process_recording for $1\n');
PROCESS_RECORDING('${SESSION}', ...
    'RAW_DATA_PATH', '$RAW_PATH', ...
    'OUT_DATA_PATH', '$OUT_PATH', ...
    'NEVUTIL_PATH', '$NEV_PATH', ... 
    'NASNET_PATH', '$NET_PATH', ...
    'HELPERS_PATH', '$HELP_PATH', ...
    'SORTER_PATH', '$SORTER_PATH');
exit
EOF

#################################################################
