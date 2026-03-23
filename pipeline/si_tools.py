# si_utils.py
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d

from spikeinterface.core import set_global_job_kwargs, BaseRecording
from spikeinterface.extractors import get_neo_streams, read_spikeglx
from spikeinterface.preprocessing import apply_preprocessing_pipeline, bandpass_filter, compute_motion
from spikeinterface.core.motion import Motion
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
from spikeinterface.sortingcomponents.motion import estimate_motion, correct_motion_on_peaks, InterpolateMotionRecording

cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
global_job_kwargs = dict(n_jobs=int(cpus - 1), chunk_duration='5s', progress_bar=True)
set_global_job_kwargs(**global_job_kwargs)

def load_or_compute_peaks(preprocess_path, recording, protocol):
    job_kwargs = dict(chunk_duration="2s", n_jobs=8, progress_bar=True)

    peaks_path = Path(preprocess_path).parent / "peaks.npy"
    peak_locs_path = Path(preprocess_path).parent / "peak_locations.npy"

    if peaks_path.exists() and peak_locs_path.exists():
        print("Loading existing peaks and peak locations.")
        peaks = np.load(peaks_path, allow_pickle=True)
        peak_locations = np.load(peak_locs_path, allow_pickle=True)
    else:
        print("Calculating peaks and peak locations...")
        peaks = detect_peaks(recording=recording, **protocol['detect_kwargs'], **job_kwargs)
        peak_locations = localize_peaks(recording, peaks, **protocol['localize_peaks_kwargs'], **job_kwargs)
        
        np.save(peaks_path, peaks)
        np.save(peak_locs_path, peak_locations)

    return peaks, peak_locations

def detect_motion_cutoffs(recording, self):
    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
    job_kwargs = dict(chunk_duration="2s", n_jobs=int(cpus - 2), progress_bar=True)
    
    save_path = Path(self.shake_path)
    if not (save_path / 'motion.npy').is_file():
        pp_recording = apply_preprocessing_pipeline(recording, self.protocol["shake_trimming"]["preprocessing"])
        
        depths = recording.get_channel_locations()
        max_depth = int(depths.max(axis=0)[1]) 

        est_kwargs = self.protocol["shake_trimming"].get("estimate_motion_kwargs", {})
        est_kwargs["win_scale_um"] = max_depth

        motion = compute_motion(pp_recording,
                    preset=self.protocol["shake_trimming"].get("estimate_motion_kwargs", {}).pop("method", "medicine"),
                    detect_kwargs=self.protocol["shake_trimming"].get("detect_kwargs", {}),
                    localize_peaks_kwargs=self.protocol["shake_trimming"].get("localize_peaks_kwargs", {}),
                    estimate_motion_kwargs=est_kwargs, **job_kwargs
                    )
         
        np.save(save_path / 'motion.npy', motion.displacement[0])
        np.save(save_path / 'time_bins.npy', motion.temporal_bins_s[0])
        np.save(save_path / 'depth_bins.npy', motion.spatial_bins_um)
    
    motion = np.load(save_path / 'motion.npy')
    time_bins = np.load(save_path / 'time_bins.npy')
    depth_bins = np.load(save_path / 'depth_bins.npy')

    time_bins = time_bins - time_bins[0]       # start time (in sec) at 0
    motion = (motion - motion[0]).squeeze()    # initial motion at 0 µm

    crop_endSec = find_motion_window(time_bins, motion, **self.protocol["shake_trimming"]["window_params"])

    return crop_endSec

def find_motion_window(time_bins, motion,
                       min_duration=1800,
                       ratio_threshold=6,
                       smooth_sigma=50,
                       padding_sec=5,
                       min_jump_um=10):
    """
    Find the motion window based on large jumps in the motion trace.
    - end_time is determined by the first jump that exceeds the threshold
      while respecting min_duration.
    """
    time_bins = np.asarray(time_bins)
    motion = np.asarray(motion).squeeze()

    motion_s = gaussian_filter1d(motion, sigma=smooth_sigma)
    dm = np.abs(np.diff(motion_s))

    typical = np.median(dm)
    if typical == 0:
        typical = np.mean(dm)

    jump_idx = np.where(dm > ratio_threshold * typical)[0]

    # always start at the first time bin
    start_time = time_bins[0]
    end_time = time_bins[-1]

    if len(jump_idx) == 0:
        return time_bins[-1]

    for idx in jump_idx:
        if idx < 10:
            continue

        jump_time = time_bins[idx]

        # compute mean motion before and after the jump
        mean_before = np.mean(motion_s[max(0, idx-1000):idx])
        mean_after  = np.mean(motion_s[idx:min(len(motion_s), idx+1000)])

        if abs(mean_after - mean_before) < min_jump_um:
            continue  # ignore small jumps

        # enforce minimum duration
        if jump_time - start_time >= min_duration:
            end_time = jump_time
            break

    # final safeguard
    if end_time - start_time < min_duration:
        return time_bins[-1]

    return end_time - padding_sec
 
def run_motion_correction(raw_recording, protocol, preprocess_path):
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


def apply_pp_pipeline(raw_recording, protocol):
    pp_recording = apply_preprocessing_pipeline(raw_recording, protocol["preprocessing"])
    pp_recording = pp_recording.astype(float)

    return pp_recording

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

