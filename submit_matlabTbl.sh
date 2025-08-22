#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=0-00:59:00
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=pipe-mtlb
#SBATCH --error=/ix1/pmayo/outfiles/out_%A.out
#SBATCH --output=/ix1/pmayo/outfiles/out_%A.out
#SBATCH --mail-type=done,fail
#SBATCH --mail-user=knoneman@pitt.edu

# ----- Load environment -----
module purge
module load matlab/R2023a
module load python/ondemand-jupyter-python3.9

ENV_PATH=$(python -c "import config; print(config.ENV_PATH)")
source activate "$ENV_PATH"

RAW_DATA_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
OUT_DATA_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
NEV_PATH=$(python -c "import config; print(config.NEVUTIL_PATH)")

# ----- Specify inputs -----
echo "======================================================"

SESSION="${1}"
PROBE_ID=$SLURM_ARRAY_TASK_ID
PROTOCOL="np_medicine.json"

echo "SESSION    =  '$SESSION'"
echo "PROBE_ID   =  $PROBE_ID"
echo "PROTOCOL   =  $PROTOCOL"

echo "======================================================"

#################################################################
################### RUN SI to get PROFILE #######################

echo "Running matlab pipeline........................"
SORTER_PATH=$($CONDA_PREFIX/bin/python -c "
from main_pipeline import profile_to_mat

print(profile_to_mat(
    '${SESSION}', 
    protocol='${PROTOCOL}' 
))
")


echo "SORTER_PATH    =  '$SORTER_PATH'"
echo "======================================================"

#################################################################
################### RUN MATLAB SCRIPT #######################

RAW_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
OUT_PATH=$(python -c "import config; print(config.RAW_DATA_PATH)")
NEV_PATH=$(python -c "import config; print(config.NEVUTIL_PATH)")

matlab -nodisplay <<EOF
try
    addpath(genpath('matlab'));   % add the subdirectory to path
    fprintf('Running process_fullRecording for $1\n');
    process_fullRecording('${SESSION}', ...
        'RAW_DATA_PATH', '$RAW_PATH', ...
        'OUT_DATA_PATH', '$OUT_PATH', ...
        'NEVUTIL_PATH', '$NEV_PATH', ...
        'SORTER_PATH', '$SORTER_PATH');
catch err
    disp('ERROR in process_fullRecording:');
    disp(getReport(err));
    exit(1);
end
exit
EOF

#################################################################
echo "DONE"

crc-job-stats
