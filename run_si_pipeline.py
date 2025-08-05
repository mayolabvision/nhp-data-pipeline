# General packages
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import numpy as np
from pathlib import Path
import shutil
import json
import sys
import logging

# Spike Interface packages
import spikeinterface.full as si
from spikeinterface.sortingcomponents.motion import correct_motion_on_peaks, interpolate_motion
from spikeinterface.sorters import run_sorter

# Custom preprocessing/plotting functions
from si_utils import *

####################################################################################################


############################################################################################################################################

def run_si_preprocess(session_name, probe_id=0, drift_correct='none', remove_bad_channels=True, raw_data_path='/ix1/pmayo/lab_NHPdata/'):
    
    global_job_kwargs = dict(n_jobs=48, chunk_duration='1s', progress_bar=True)
    si.set_global_job_kwargs(**global_job_kwargs)
    
    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    data_folder, preprocess_folder, _, metadata = make_folder_paths(raw_data_path, session_name, probe_id, drift_correct)
    if preprocess_folder.exists():
        shutil.rmtree(preprocess_folder)

    preprocess_folder.mkdir(parents=True, exist_ok=True)
    figs_folder = preprocess_folder / 'figs'
    os.makedirs(figs_folder, exist_ok=True)

    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    print("Loading raw data...........................")
    if metadata["probe_types"][probe_id] == "neuropixel":
        stream_names, stream_ids = si.get_neo_streams('spikeglx', data_folder)
        RAW_RECORDING = si.read_spikeglx(data_folder, stream_name=f'imec{probe_id}.ap', load_sync_channel=False)

    print("Sampling frequency:", RAW_RECORDING.get_sampling_frequency())
    print("Number of channels:", RAW_RECORDING.get_num_channels())
    print("Number of segments:", RAW_RECORDING.get_num_segments())
    print("Number of samples:", RAW_RECORDING.get_num_samples(segment_index=0))
    print("Duration (min):", round((RAW_RECORDING.get_num_samples(segment_index=0)/RAW_RECORDING.get_sampling_frequency())/60))
    print("Data dtype:", RAW_RECORDING.get_dtype())

    #-------------------------------------------------------------------------------------------------------------------------------------- 

    print("Plotting noise distributions...........................")
    plot_noise_hists(RAW_RECORDING, save_path=os.path.join(figs_folder, 'noise_dists.png'))

    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    print("Filtering recording...........................")
    FILT_RECORDING, bad_channels = preprocess_raw_recording(RAW_RECORDING, remove_bad_channels=remove_bad_channels)
    print("Bad channels:", bad_channels)

    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    rec = FILT_RECORDING.astype(float)
    if drift_correct.lower() != 'none':
        motion_folder = preprocess_folder / 'motion'
    
        print("Estimate motion from recording...........................")
        PROC_RECORDING, motion_info = si.correct_motion(rec, preset=drift_correct, 
                                                        interpolate_motion_kwargs={'border_mode':'force_extrapolate'}, 
                                                        folder=motion_folder, output_motion_info=True)
        PROC_RECORDING = PROC_RECORDING.astype(int)
   
        print("Plot motion correction...........................")
        plot_motion_correction(FILT_RECORDING, motion_info, save_path=os.path.join(figs_folder, 'motion_info.png'))

    else:
        print("Skipping motion correction...........................")
        PROC_RECORDING = rec.astype(int) 
        
    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    print("Saving drift-corrected recording to disk...........................")
    PROC_RECORDING = PROC_RECORDING.save(folder=preprocess_folder / 'output', format='binary', dtype='int16', overwrite=True)

    #-------------------------------------------------------------------------------------------------------------------------------------- 

############################################################################################################################################

def run_si_spikesort(session_name, probe_id=0, sorter_type='kilosort4', drift_correct='none', raw_data_path='/ix1/pmayo/lab_NHPdata/'):
    
    global_job_kwargs = dict(n_jobs=12, chunk_duration='1s', progress_bar=True)
    si.set_global_job_kwargs(**global_job_kwargs)
    
    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    print("Defining sorting parameters...........................")

    default_params = si.get_default_sorter_params(sorter_name_or_class=sorter_type)
    print(f"--------------{sorter_type}--------------\n{default_params}\n------------------------------------------------------\n")

    custom_params = {
        'do_correction': False,
        'clear_cache': True
    }
    print(f"--------- custom parameters ---------\n{custom_params}\n-------------------------------------------------------------")
    
    run_suffix = format_param_suffix(default_params, custom_params)
    run_suffix = run_suffix or 'defaults'

    _, preprocess_folder, sorter_folder, metadata = make_folder_paths(raw_data_path, session_name, probe_id, drift_correct, 
                                                                      sorter_type=sorter_type, run_suffix=run_suffix)
    
    os.makedirs(sorter_folder.parent, exist_ok=True)
    #setup_logging(sorter_folder)
    
    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    print("Loading preprocessed data...........................")
    PROC_RECORDING = si.load(preprocess_folder / 'output')

    #-------------------------------------------------------------------------------------------------------------------------------------- 

    print("Running spike sorting...............................")
    SORTING = si.run_sorter(sorter_type, PROC_RECORDING, remove_existing_folder=True, folder=sorter_folder,
                            docker_image=False, verbose=True, **custom_params)

    #-------------------------------------------------------------------------------------------------------------------------------------- 

############################################################################################################################################

