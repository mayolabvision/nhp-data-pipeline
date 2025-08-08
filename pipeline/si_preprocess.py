# si_utils.py
import os
from pathlib import Path
import numpy as np

from spikeinterface.core import set_global_job_kwargs
from spikeinterface.extractors import get_neo_streams, read_spikeglx
from spikeinterface.preprocessing import apply_preprocessing_pipeline
from spikeinterface.core.motion import Motion
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
from spikeinterface.sortingcomponents.motion import estimate_motion, correct_motion_on_peaks, InterpolateMotionRecording

global_job_kwargs = dict(n_jobs=int(int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))-1), chunk_duration='1s', progress_bar=True)
set_global_job_kwargs(**global_job_kwargs)

def load_or_compute_peaks(preprocess_path, recording, protocol):
    peaks_path = Path(preprocess_path).parent / "peaks.npy"
    peak_locs_path = Path(preprocess_path).parent / "peak_locations.npy"

    if peaks_path.exists() and peak_locs_path.exists():
        print("Loading existing peaks and peak locations.")
        peaks = np.load(peaks_path, allow_pickle=True)
        peak_locations = np.load(peak_locs_path, allow_pickle=True)
    else:
        print("Calculating peaks and peak locations...")
        peaks = detect_peaks(recording=recording, **protocol['detect_kwargs'])
        peak_locations = localize_peaks(recording, peaks, **protocol['localize_peaks_kwargs'])
        
        np.save(peaks_path, peaks)
        np.save(peak_locs_path, peak_locations)

    return peaks, peak_locations

def run_preprocessing_with_motion_correction(data_path, probe_id, protocol, preprocess_path):
    os.makedirs(preprocess_path, exist_ok=True)

    if not (preprocess_path / 'motion.npy').exists():
        # This gets stream info if needed (unused here but might be needed)
        stream_names, stream_ids = get_neo_streams('spikeglx', data_path)

        raw_recording = read_spikeglx(data_path, stream_name=f'imec{probe_id}.ap', load_sync_channel=False)
        pp_recording1 = apply_preprocessing_pipeline(raw_recording, protocol['motion_correction']['preprocessing'])

        peaks, peak_locations = load_or_compute_peaks(preprocess_path, pp_recording1, protocol['motion_correction'])

        motion_object = estimate_motion(
                            recording=pp_recording1,
                            peaks=peaks,
                            peak_locations=peak_locations,
                            **protocol['motion_correction']['estimate_motion_kwargs'])

        np.save(preprocess_path / 'motion.npy', motion_object.displacement)
        np.save(preprocess_path / 'time_bins.npy', motion_object.temporal_bins_s)
        np.save(preprocess_path / 'depth_bins.npy', motion_object.spatial_bins_um)

        mc_peak_locations = correct_motion_on_peaks(
                                peaks=peaks,
                                peak_locations=peak_locations,
                                motion=motion_object,
                                recording=pp_recording1)
        np.save(preprocess_path / 'mc_peak_locations.npy', mc_peak_locations)

    else:
        # Load in estimated motion
        motion = np.load(preprocess_path / 'motion.npy')
        time_bins = np.load(preprocess_path / 'time_bins.npy')
        depth_bins = np.load(preprocess_path / 'depth_bins.npy')
        
        motion_object = Motion(
            displacement=motion,
            temporal_bins_s=time_bins,
            spatial_bins_um=depth_bins,
        )

    # Preprocessing recording
    pp_recording = apply_preprocessing_pipeline(raw_recording, protocol['preprocessing'])

    # Use interpolation to correct for estimated motion
    mc_recording = InterpolateMotionRecording(
        pp_recording.astype(float), 
        motion_object,
        **protocol['motion_correction']['interpolate_motion_kwargs']
    )
    
    return mc_recording

def run_preprocessing_without_motion_correction(data_path, probe_id, protocol, preprocess_path):
    stream_names, stream_ids = get_neo_streams('spikeglx', data_path)
    raw_recording = read_spikeglx(data_path, stream_name=f'imec{probe_id}.ap', load_sync_channel=False)

    pp_recording = apply_preprocessing_pipeline(raw_recording, protocol['preprocessing'])

    return pp_recording.astype(float)

def save_processed_recording(recording, preprocess_path):
    recording = recording.astype(int)
    recording = recording.save(folder=preprocess_path, format='binary', dtype='int16', overwrite=True)
    return recording
