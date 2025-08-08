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
from spikeinterface.core import set_global_job_kwargs
from spikeinterface.extractors import get_neo_streams, read_spikeglx
from spikeinterface.preprocessing import apply_preprocessing_pipeline
from spikeinterface.core.motion import Motion
from spikeinterface.sortingcomponents.motion.motion_interpolation import InterpolateMotionRecording
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
from spikeinterface.sortingcomponents.motion import estimate_motion
from spikeinterface.sorters import run_sorter

import medicine

# Custom preprocessing/plotting functions
from utils import *

############################################################################################################################################

def make_kilosortChanMap(meta_path, probe_id):
    meta_path = Path(meta_path)
    with open(meta_path, 'r') as f:
        metadata = json.load(f)
        
    session_name = meta_path.parent.name
    probe_type = metadata['probe_types'][probe_id]
    
    print(probe_type)
    if probe_type == 'neuropixel':
        probe_path = meta_path.parent / f"{session_name}_imec{probe_id}" / f"{session_name}_t0.imec{probe_id}.ap_kilosortChanMap.mat"
        print(probe_path.parent / f"{session_name}_t0.imec{probe_id}.ap.meta")

        if not probe_path.exists():
            MetaToCoords(probe_path.parent / f"{session_name}_t0.imec{probe_id}.ap.meta", 
                         1, badChan=np.zeros((0), dtype='int'), destFullPath='', showPlot=True)

############################################################################################################################################

def run_si_preprocess(session, probe_id=0, protocol_file='np_protocol.json', raw_data_path='/ix1/pmayo/lab_NHPdata/'): 
    
    global_job_kwargs = dict(n_jobs=int(int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))-1), chunk_duration='1s', progress_bar=True)
    set_global_job_kwargs(**global_job_kwargs)
   
    #-------------------------------------------------------------------------------------------------------------------------------------- 
     
    # Load in protocol that contains info about current run of pipeline
    with open(os.path.join(os.getcwd(), "protocols", protocol_file), 'r') as f:
        protocol = json.load(f)
    
    # Initial steps, setting up directories based on protocol/metadata
    metadata, data_folder, motion_folder, preprocess_folder, _, _, _ = make_folder_paths(session, protocol, probe_id=probe_id, raw_data_path=raw_data_path)

    #-------------------------------------------------------------------------------------------------------------------------------------- 
    if os.path.exists(preprocess_folder) and any(f.endswith('.raw') for f in os.listdir(preprocess_folder)): 
        print("preprocess_folder exists and contains a .raw file.")
    else:
        print("-------------------LOADING RAW RECORDING-------------------")
        if metadata["probe_types"][probe_id] == "neuropixel":
            stream_names, stream_ids = get_neo_streams('spikeglx', data_folder)
            RAW_RECORDING = read_spikeglx(data_folder, stream_name=f'imec{probe_id}.ap', load_sync_channel=False)

        print_recording_details(RAW_RECORDING) 
        #-------------------------------------------------------------------------------------------------------------------------------------- 
        
        print("-------------------PREPROCESSING RECORDING-------------------")
        PP_RECORDING = apply_preprocessing_pipeline(RAW_RECORDING, protocol['preprocessing'])
        
        #PP_RECORDING = PP_RECORDING.frame_slice(start_frame=0,
        #                              end_frame=int(900*30000))
        
        print_recording_details(PP_RECORDING) 
        #-------------------------------------------------------------------------------------------------------------------------------------- 

        if protocol["motion_correction"]["drift_preset"] != "none":
            # Check that motion estimation hasn't already been for this session/protocol
            if not (Path(motion_folder).is_dir() and all((Path(motion_folder) / f).exists() for f in ['motion.npy'])):
                print("-------------------ESTIMATING MOTION-------------------")
                motion_params = {
                k: v for k, v in protocol["motion_correction"].items() if k != "drift_preset"}
                    
                peaks, peak_locations = estimate_peaks(motion_folder, PP_RECORDING)
                
                if protocol["motion_correction"]["drift_preset"] == "medicine":
                    # Run MEDiCINe to estimate motion
                    medicine.run_medicine(
                        peak_amplitudes=peaks['amplitude'],
                        peak_depths=peak_locations['y'],
                        peak_times=peaks['sample_index'] / PP_RECORDING.get_sampling_frequency(),
                        output_dir=motion_folder,
                        **motion_params
                    )

                else:
                    motion = estimate_motion(recording=PP_RECORDING,
                             peaks=peaks,
                             peak_locations=peak_locations)

                    np.save(motion_folder / 'motion.npy', motion['displacement'])
                    np.save(motion_folder / 'time_bins.npy', motion['temporal_bins_s'])
                    np.save(motion_folder / 'depth_bins.npy', motion['spatial_bins_um'])

            print("-------------------LOADING MOTION-------------------")
            # Load motion estimated by MEDiCINe
            motion = np.load(motion_folder / 'motion.npy')
            time_bins = np.load(motion_folder / 'time_bins.npy')
            depth_bins = np.load(motion_folder / 'depth_bins.npy')

            print("-------------------CORRECTING MOTION-------------------")
            motion_object = Motion(
                displacement=motion,
                temporal_bins_s=time_bins,
                spatial_bins_um=depth_bins,
            )
            
            MC_RECORDING = InterpolateMotionRecording(
                PP_RECORDING.astype(float),
                motion_object,
                border_mode='force_extrapolate',
            )
        else:
            MC_RECORDING = PP_RECORDING.astype(float)

        MC_RECORDING = MC_RECORDING.astype(int)
        print_recording_details(MC_RECORDING) 

        print("-------------------SAVING RECORDING-------------------")
        MC_RECORDING = MC_RECORDING.save(folder=preprocess_folder, format='binary', dtype='int16', overwrite=True)


#-------------------------------------------------------------------------------------------------------------------------------------- 

def estimate_peaks(motion_folder,recording):
    
    peaks_path = Path(motion_folder).parent / "peaks.npy"
    peak_locs_path = Path(motion_folder).parent / "peak_locations.npy"
    
    if peaks_path.exists() and peak_locs_path.exists():
        peaks = np.load(peaks_path, allow_pickle=True)
        peak_locations = np.load(peak_locs_path, allow_pickle=True)
        print("Loaded existing peaks and peak locations.")
    else: 
        print("Calculating peaks and peak locations...")

        peaks = detect_peaks(recording=recording, method='locally_exclusive')
        peak_locations = localize_peaks(recording, peaks, method="monopolar_triangulation")
        
        np.save(peaks_path, peaks)
        np.save(peak_locs_path, peak_locations)
    
    return peaks, peak_locations

############################################################################################################################################

def run_si_sorting(session, probe_id=0, protocol_file='np_protocol.json', raw_data_path='/ix1/pmayo/lab_NHPdata/'): 

    #-------------------------------------------------------------------------------------------------------------------------------------- 
     
    # Load in protocol that contains info about current run of pipeline
    with open(os.path.join(os.getcwd(), "protocols", protocol_file), 'r') as f:
        protocol = json.load(f)
    
    # Initial steps, setting up directories based on protocol/metadata
    metadata, data_folder, motion_folder, sorter_folder, _, _ = make_folder_paths(session, protocol, probe_id=probe_id, raw_data_path=raw_data_path)

    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    print("-------------------LOADING RAW RECORDING-------------------")
    if metadata["probe_types"][probe_id] == "neuropixel":
        stream_names, stream_ids = get_neo_streams('spikeglx', data_folder)
        RAW_RECORDING = read_spikeglx(data_folder, stream_name=f'imec{probe_id}.ap', load_sync_channel=False)

    print_recording_details(RAW_RECORDING) 
    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    print("-------------------PREPROCESSING RECORDING-------------------")
    PP_RECORDING = apply_preprocessing_pipeline(RAW_RECORDING, protocol['preprocessing'])
    
    #PP_RECORDING = PP_RECORDING.frame_slice(start_frame=0,
    #                              end_frame=int(900*30000))
    
    print_recording_details(PP_RECORDING) 
    #-------------------------------------------------------------------------------------------------------------------------------------- 

    if protocol["motion_correction"]["drift_preset"] != "none":
        print("-------------------LOADING MOTION-------------------")
        # Load motion estimated by MEDiCINe
        motion = np.load(motion_folder / 'motion.npy')
        time_bins = np.load(motion_folder / 'time_bins.npy')
        depth_bins = np.load(motion_folder / 'depth_bins.npy')

        print("-------------------CORRECTING MOTION-------------------")
        motion_object = Motion(
            displacement=motion,
            temporal_bins_s=time_bins,
            spatial_bins_um=depth_bins,
        )
        
        MC_RECORDING = InterpolateMotionRecording(
            PP_RECORDING.astype(float),
            motion_object,
            border_mode='force_extrapolate',
        )

    else:
        MC_RECORDING = PP_RECORDING.astype(float)

    print_recording_details(MC_RECORDING) 
    
    #-------------------------------------------------------------------------------------------------------------------------------------- 
    
    print("-------------------SPIKE SORTING-------------------")
    sorting_params = protocol['sorting']
    SORTING = run_sorter(recording=MC_RECORDING, **sorting_params)    

############################################################################################################################################

