#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=0-08:00:00
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
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

################################################################
##################### MAKE PROBE MAP ###########################

echo "Making probe map"
$CONDA_PREFIX/bin/python -c "
from run_si_pipeline import make_kilosortChanMap

make_kilosortChanMap(
    '${RAW_DATA_PATH}/${SESSION}/metadata.json',
    probe_id=int('$PROBE_ID')
)"

echo "========================"

#################################################################
##################### GET SYNC PULSES ###########################

echo "Pulling out sync pulses"
PROBE_TYPE=$(jq -r ".probe_types[${PROBE_ID}]" "${RAW_DATA_PATH}/${SESSION}/metadata.json")

if [[ "$PROBE_TYPE" == "neuropixel" ]]; then
    run_name="${SESSION%_g0}"
    RUNIT_PATH="/ix1/pmayo/packages/CatGT-linux/runit.sh"

    if [ ! -d "/ix1/pmayo/lab_NHPdata/${SESSION}/catgt_${SESSION}" ]; then
        echo "Directory does NOT exist. Running runit.sh..."
        ${RUNIT_PATH} '-dir=/ix1/pmayo/lab_NHPdata -run='$run_name' -g=0 -t=0,0 -t_miss_ok -ni -prb=0 -bf=0,0,-1,0,9,1 -dest=/ix1/pmayo/lab_NHPdata/'$SESSION''
    else
        echo "catgt_'$SESSION' exists. Skipping runit.sh."
    fi
fi

echo "========================"

#################################################################
##################### RUN SI (PP + Sort) ########################

echo "Running spike interface"
$CONDA_PREFIX/bin/python -c "
from run_si_pipeline import run_si_preprocess

run_si_preprocess(
    '${SESSION}', 
    probe_id=int('$PROBE_ID'), 
    protocol_file='${PROTOCOL_FILE}', 
    raw_data_path='${RAW_DATA_PATH}'
)"

echo "========================"

#################################################################
echo "DONE"

crc-job-stats
