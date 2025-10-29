# si_utils.py
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

from spikeinterface.core import set_global_job_kwargs, BaseRecording
from spikeinterface.extractors import get_neo_streams, read_spikeglx
from spikeinterface.preprocessing import apply_preprocessing_pipeline, compute_motion
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

def detect_excessive_motion(raw_recording, preprocess_path, threshold_um=1000, min_duration_sec=30):
    
    if not (preprocess_path.parent.parent / 'rigid_fast_motion.npy').is_file():
        preprocessing_dict = {
        'bandpass_filter': {'freq_min': 300, 'freq_max': 5000, 'dtype': 'int16'},
        'phase_shift': {},
        'common_reference': {'operator': 'median', 'reference': 'global'}
        }

        pp_recording = apply_preprocessing_pipeline(raw_recording, preprocessing_dict)
       
        motion = compute_motion(pp_recording, preset="rigid_fast") 
        
        np.save(preprocess_path.parent.parent / 'rigid_fast_motion.npy', motion.displacement)
        np.save(preprocess_path.parent.parent / 'rigid_fast_time_bins.npy', motion.temporal_bins_s)
        np.save(preprocess_path.parent.parent / 'rigid_fast_depth_bins.npy', motion.spatial_bins_um)
    
    motion = np.load(preprocess_path.parent.parent / 'rigid_fast_motion.npy')
    time_bins = np.load(preprocess_path.parent.parent / 'rigid_fast_time_bins.npy')
    depth_bins = np.load(preprocess_path.parent.parent / 'rigid_fast_depth_bins.npy') 
    
    time_bins = time_bins[0] - time_bins[0][0]       # start time (in sec) at 0
    motion = (motion[0] - motion[0][0]).squeeze()    # initial motion at 0 µm

    motion_smooth = savgol_filter(motion, window_length=11, polyorder=3)

    # --- calculate onset time ---
    onset_time_sec = None
    onset_time_frame = None

    # mask: don't worry about initial noise issues
    mask = time_bins >= 1200

    above_thresh = (np.abs(motion_smooth) > threshold_um) & mask

    if np.any(above_thresh):
        idx = np.where(above_thresh)[0]

        # Find breaks between contiguous stretches
        breaks = np.where(np.diff(idx) > 1)[0] + 1
        segments = np.split(idx, breaks)

        for segment in segments:
            duration = time_bins[segment[-1]] - time_bins[segment[0]]
            if duration >= min_duration_sec:
                onset_time_sec = time_bins[segment[0]]
                break

    return onset_time_sec 

def run_preprocessing_with_motion_correction(raw_recording, protocol, preprocess_path):
    os.makedirs(preprocess_path, exist_ok=True)
    
    if not (preprocess_path / 'motion.npy').exists():
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
            displacement=[motion[0]],
            temporal_bins_s=[time_bins[0]],
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

def run_preprocessing_without_motion_correction(raw_recording, protocol, preprocess_path):
    pp_recording = apply_preprocessing_pipeline(raw_recording, protocol['preprocessing'])

    return pp_recording.astype(float)

def save_processed_recording(recording, preprocess_path):
    recording = recording.astype(int)
    recording = recording.save(folder=preprocess_path, format='binary', dtype='int16', overwrite=True)
    return recording

######################################################################################################

def find_missing_extensions(extensions_dir, requested_extensions):
    """
    Compare the requested extensions against the existing ones in `extensions_dir`
    and return a dictionary of the missing extensions.
    """
    # Get all folder names in the extensions directory
    existing = {name for name in os.listdir(extensions_dir)
                if os.path.isdir(os.path.join(extensions_dir, name))}
    
    # Keep only requested ones that are NOT in existing
    missing = {k: v for k, v in requested_extensions.items() if k not in existing}
    
    return missing

######################################################################################################
def add_extension_arrays_to_metrics(extensions_path: Path, metrics: pd.DataFrame) -> pd.DataFrame:
    """
    Loop through extension directories in `extensions_path`, find .npy files
    whose shape[0] matches the number of rows in `metrics`, and add them as new columns.

    Each .npy filename (without extension) becomes the column name.
    If the file has extra dimensions beyond the first, the row stores that subarray.
    """
    n_rows = metrics.shape[0]

    skip_dirs = {"templates", "correlograms", "principal_components"}

    for ext_dir in extensions_path.iterdir():
        if ext_dir.is_dir() and ext_dir.name not in skip_dirs:
            for npy_file in ext_dir.glob("*.npy"):
                data = np.load(npy_file, allow_pickle=True)

                # Only add if first dimension matches
                if data.shape[0] == n_rows:
                    col_name = npy_file.stem  # filename without ".npy"

                    # If data is more than 1D, keep row-wise slices as arrays
                    if data.ndim > 1:
                        metrics.loc[:, col_name] = [row for row in data]
                    else:
                        metrics[col_name] = data

    return metrics

