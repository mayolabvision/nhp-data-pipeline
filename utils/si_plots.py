import os
import numpy as np
import matplotlib.pyplot as plt
import spikeinterface.full as si
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
from spikeinterface.sortingcomponents.motion import correct_motion_on_peaks

job_kwargs = dict(n_jobs=-1, chunk_duration='1s', progress_bar=True)

def plot_noise_hists(rec, save_path=None):
    # we can estimate the noise on the scaled traces (microV) or on the raw one (which is in our case int16).
    noise_levels_microV = si.get_noise_levels(rec, return_scaled=True, **job_kwargs)
    noise_levels_int16 = si.get_noise_levels(rec, return_scaled=False, **job_kwargs)
    
    # Create subplots side by side
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    
    # Plot scaled traces (microV)
    axes[0].hist(noise_levels_microV, bins=np.arange(5, 30, 2.5))
    axes[0].set_title('scaled traces (microV)')
    axes[0].set_xlabel('noise [microV]')
    
    # Plot raw traces (int16)
    axes[1].hist(noise_levels_int16, bins=np.arange(5, 30, 2.5))
    axes[1].set_title('raw traces (int16)')
    axes[1].set_xlabel('noise [int16 units]')

    if save_path is not None:
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close(fig) 

def plot_peaks_from_recording(rec,save_path=None):
    noise_levels_int16 = si.get_noise_levels(rec, return_scaled=False)

    peaks = detect_peaks(rec,  method='locally_exclusive', noise_levels=noise_levels_int16,
                         detect_threshold=5, radius_um=50., **job_kwargs)
    peak_locations = localize_peaks(rec, peaks, method='center_of_mass', radius_um=50., **job_kwargs)

    # check for drifts
    fs = rec.sampling_frequency
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(peaks['sample_index'] / fs, peak_locations['y'], color='k', marker='.',  alpha=0.002)

    if save_path is not None:
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close(fig) 

    return peaks, peak_locations

def plot_motion_correction(rec,motion_info,save_path=None):
    
    fig = plt.figure(figsize=(14, 8))
    si.plot_motion_info(motion_info, rec, figure=fig, depth_lim=(400, 1000), color_amplitude=True, amplitude_cmap="inferno", scatter_decimate=10,)

    if save_path is not None:
        plt.savefig(save_path)
        plt.close(fig) 

def plot_peaks_with_drift_correction(peaks,peak_locations,motion,rec,save_path=None):
    
    corrected_peak_locations = correct_motion_on_peaks(peaks, peak_locations, motion, rec)

    # check for drifts
    fs = rec.sampling_frequency
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(peaks['sample_index'] / fs, corrected_peak_locations['y'], color='k', marker='.', alpha=0.002)

    if save_path is not None:
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close(fig) 


