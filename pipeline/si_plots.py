# si_plots.py
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d

from spikeinterface.core import set_global_job_kwargs, BaseRecording
from spikeinterface.extractors import get_neo_streams, read_spikeglx
from spikeinterface.preprocessing import apply_preprocessing_pipeline
from spikeinterface.core.motion import Motion
from spikeinterface.core import get_noise_levels
import spikeinterface.widgets as sw
from spikeinterface.widgets import UnitSummaryWidget

from medicine.plotting import plot_motion_correction

from kilosort.io import load_ops
from kilosort.data_tools import (
    mean_waveform, cluster_templates, get_good_cluster, get_cluster_spikes,
    get_spike_waveforms, get_best_channels
    )

global_job_kwargs = dict(n_jobs=int(int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))-1), chunk_duration='1s', progress_bar=True)
set_global_job_kwargs(**global_job_kwargs)

def plot_probe_motion(profile):
    motion = np.load(Path(profile.shake_path) / 'motion.npy')
    time_bins = np.load(Path(profile.shake_path) / 'time_bins.npy')
    
    motion_s = gaussian_filter1d(motion, sigma=profile.protocol['shake_trimming']['window_params']['smooth_sigma'])

    fig, (ax1) = plt.subplots(1, 1, figsize=(15, 4))

    ax1.plot(time_bins, motion, color="black", linewidth=2)
    ax1.plot(time_bins, motion_s, color="red", linewidth=2)
 
    ax1.axvline(0, linestyle='--')
    if profile.crop_endSec is not None:
        ax1.axvline(profile.crop_endSec, linestyle='--')

    ax1.set_xlabel("Time [sec]")
    ax1.set_ylabel("Motion [$\\mu$m]")

    fig.suptitle(f"{profile.session} --- {profile.metadata['probe_label'][profile.probe_id]}", fontsize=16, y=0.98)
    plt.tight_layout()
    fig.text(0.5, 0.88, f"{profile.metadata['hardware_config'][profile.probe_id]}, {profile.metadata['probe_config'][profile.probe_id]}: depth = {profile.metadata['probe_depth_mm'][profile.probe_id]}mm, grid hole = {profile.metadata['probe_gridHole'][profile.probe_id]}", ha='center', fontsize=12)
    
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.97])  # leave space for labels
    fig.savefig(Path(profile.figs_path) / "probe_motion_cutoff.png", dpi=300, bbox_inches="tight")

def plot_preprocessing_steps(raw_recording, profile):
    # collect preprocessing steps
    profile.protocol["preprocessing"].pop("motion_crop", None)
    pp_steps = list(profile.protocol['preprocessing'].items())

    # make subplots: one for raw + one for each preprocessing step
    fig1, axs = plt.subplots(
        ncols=len(pp_steps) + 1, 
        figsize=(30, 10),
        squeeze=False,
        sharex=True,
        sharey=True)

    # plot raw
    im = sw.plot_traces(
        raw_recording, backend='matplotlib', clim=(-50, 50), 
        ax=axs[0, 0], mode="map")
    axs[0, 0].set_title("raw")

    # loop through steps progressively
    recording = raw_recording
    for i, (step_name, step_params) in enumerate(pp_steps, start=1):
        print(f"\n=== Step {i}/{len(pp_steps)}: '{step_name}' ===")
        print(f"Parameters: {step_params}")
        
        # apply pipeline with only the first i steps
        sub_pipeline = dict(pp_steps[:i])
        recording = apply_preprocessing_pipeline(recording, sub_pipeline)

        # plot
        im = sw.plot_traces(
            recording, backend='matplotlib', clim=(-50, 50), 
            ax=axs[0, i], mode="map")
        axs[0, i].set_title(step_name)

    # axis labels (shared)
    fig1.text(0.5, 0.04, 'Time (seconds)', ha='center', fontsize=14)
    fig1.text(0.04, 0.5, 'Channels', va='center', rotation='vertical', fontsize=14)
    
    fig1.suptitle(f"{profile.session} --- {profile.metadata['probe_label'][profile.probe_id]}", fontsize=16, y=0.98)
    fig1.text(0.5, 0.94, f"{profile.metadata['hardware_config'][profile.probe_id]}, {profile.metadata['probe_config'][profile.probe_id]}: depth = {profile.metadata['probe_depth_mm'][profile.probe_id]}mm, grid hole = {profile.metadata['probe_gridHole'][profile.probe_id]}", ha='center', fontsize=12)
    
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])  # leave space for labels
    fig1.savefig(Path(profile.figs_path) / "preprocessing_steps.png", dpi=300, bbox_inches="tight")

def plot_probe_peaks(raw_recording, profile):
    peaks = np.load(Path(profile.preprocess_path).parent / "peaks.npy", allow_pickle=True)
    peak_locations = np.load(Path(profile.preprocess_path).parent / "peak_locations.npy", allow_pickle=True)

    fig2, axs = plt.subplots(
        ncols=4, 
        figsize=(15,10),
        squeeze=False,
        sharex=True,
        sharey=False,
        gridspec_kw={'width_ratios': [1, 1, 1, 0.5]})

    sw.plot_probe_map(raw_recording, ax=axs[0,0], with_channel_ids=True)
    axs[0,0].scatter(peak_locations['x'], peak_locations['y'], color='purple', alpha=0.002)
    axs[0,0].set_ylim(-100,500)

    sw.plot_probe_map(raw_recording, ax=axs[0,1], with_channel_ids=True)
    axs[0,1].scatter(peak_locations['x'], peak_locations['y'], color='purple', alpha=0.002)
    axs[0,1].set_ylim(3500,4100)

    sw.plot_probe_map(raw_recording, ax=axs[0,2], with_channel_ids=True)
    axs[0,2].scatter(peak_locations['x'], peak_locations['y'], color='purple', alpha=0.002)
    axs[0,2].set_ylim(7200,7800)

    sw.plot_peak_activity(raw_recording, peaks, ax=axs[0,3])
    axs[0,3].set_ylim(-100,7800)
    axs[0,3].set_title('')

    fig2.suptitle(f"{profile.session} --- {profile.metadata['probe_label'][profile.probe_id]}", fontsize=16, y=0.98)
    fig2.text(0.5, 0.94, f"{profile.metadata['hardware_config'][profile.probe_id]}, {profile.metadata['probe_config'][profile.probe_id]}: depth = {profile.metadata['probe_depth_mm'][profile.probe_id]}mm, grid hole = {profile.metadata['probe_gridHole'][profile.probe_id]}", ha='center', fontsize=12)

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.97])  # leave space for labels
    fig2.savefig(Path(profile.figs_path) / "probe_peaks.png", dpi=300, bbox_inches="tight")

def plot_noise_levels(raw_recording, profile):
    fig3, axs = plt.subplots(
        ncols=2, 
        figsize=(10,5),
        squeeze=False,
        sharey=True
    )

    # we can estimate the noise on the scaled traces (microV) or on the raw one (which is in our case int16).
    noise_levels_microV = get_noise_levels(raw_recording, return_in_uV=True)
    noise_levels_int16 = get_noise_levels(raw_recording, return_in_uV=False)

    _ = axs[0,0].hist(noise_levels_microV, bins=np.arange(5, 30, 2.5), color='gray', edgecolor='black')
    axs[0,0].set_xlabel('noise [microV]')

    _ = axs[0,1].hist(noise_levels_int16, bins=np.arange(0, 15, 2.5), color='gray', edgecolor='black')
    axs[0,1].set_xlabel('noise [int16]')

    fig3.text(0.04, 0.5, '# of channels', va='center', rotation='vertical', fontsize=10)
    fig3.suptitle(f"{profile.session} --- {profile.metadata['probe_label'][profile.probe_id]}", fontsize=14, y=0.98)
    fig3.text(0.5, 0.9, f"{profile.metadata['hardware_config'][profile.probe_id]}, {profile.metadata['probe_config'][profile.probe_id]}: depth = {profile.metadata['probe_depth_mm'][profile.probe_id]}mm, grid hole = {profile.metadata['probe_gridHole'][profile.probe_id]}", ha='center', fontsize=10)

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.97])  # leave space for labels
    fig3.savefig(Path(profile.figs_path) / "noise_levels.png", dpi=300, bbox_inches="tight")

def plot_motion_correction_traces(raw_recording, profile):
    motion = np.load(profile.preprocess_path / 'motion.npy')
    time_bins = np.load(profile.preprocess_path / 'time_bins.npy')
    depth_bins = np.load(profile.preprocess_path / 'depth_bins.npy')

    peaks = np.load(Path(profile.preprocess_path).parent / "peaks.npy", allow_pickle=True)
    peak_locations = np.load(Path(profile.preprocess_path).parent / "peak_locations.npy", allow_pickle=True)

    #motion_object = Motion(
    #    displacement=[motion[0]],
    #    temporal_bins_s=[time_bins[0]],
    #    spatial_bins_um=depth_bins,
    #)

    #sw.plot_motion(motion_object)

    fig4 = plot_motion_correction(peak_times = peaks['sample_index'] / raw_recording.sampling_frequency, 
                           peak_depths = peak_locations['y'], 
                           peak_amplitudes = peaks['amplitude'], 
                           time_bins = time_bins[0], 
                           depth_bins = depth_bins, 
                           motion = motion[0])

    fig4.suptitle(f"{profile.session} --- {profile.metadata['probe_label'][profile.probe_id]}", fontsize=18, y=0.98)
    fig4.text(0.5, 0.93, f"{profile.metadata['hardware_config'][profile.probe_id]}, {profile.metadata['probe_config'][profile.probe_id]}: depth = {profile.metadata['probe_depth_mm'][profile.probe_id]}mm, grid hole = {profile.metadata['probe_gridHole'][profile.probe_id]}", ha='center', fontsize=12)

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.96])  # leave space for labels
    fig4.savefig(Path(profile.figs_path) / "motion_correction.png", dpi=300, bbox_inches="tight")


def plot_analyzer_sess(analyzer, metrics, profile):
    fig = plt.figure(figsize=(12, 12))  # adjust as needed
    gs = GridSpec(10, 5, figure=fig)

    # Create axes with your specified layout
    ax0 = fig.add_subplot(gs[0:9, 0]) 
    ax1 = fig.add_subplot(gs[0:2, 1:4]) 
    ax2 = fig.add_subplot(gs[2:4, 1:4]) 
    ax3 = fig.add_subplot(gs[4:7, 1:3]) 
    ax4 = fig.add_subplot(gs[4:6, 3:4])
    ax5 = fig.add_subplot(gs[6:7, 3:4])
    ax6 = fig.add_subplot(gs[7:9, 1:2])
    ax7 = fig.add_subplot(gs[7:9, 2:3])
    ax8 = fig.add_subplot(gs[7:9, 3:4])

    # Now call your plotting functions on each axis
    sw.plot_unit_depths(analyzer, ax=ax0)

    sw.plot_all_amplitudes_distributions(analyzer, unit_ids=[0,1,2,3,4,5], ax=ax1)
    ax1.set_xlabel("clusters")
    ax1.set_xticks([])
    ax1.set_ylabel("amplitude")

    sw.plot_unit_presence(analyzer, time_range=[0, 60], ax=ax2)
    ax2.set_ylabel("clusters")
    cbar = ax2.figure.axes[-1]  # colorbar axis is usually last
    cbar.set_ylabel("normalized activity")  # label it

    sw.plot_template_similarity(analyzer, ax=ax3)
    ax3.set_xlabel("cluster 1")
    ax3.set_ylabel("cluster 2")
    cbar = ax3.figure.axes[-1]  # colorbar axis is usually last
    cbar.set_ylabel("template similarity")  # label it

    ax4.hist(metrics["firing_rate"], bins=20)
    ax4.set_xlabel("firing rate (Hz)")
    ax4.set_ylabel("number of clusters")

    ax5.hist(metrics["presence_ratio"], bins=20)
    ax5.set_xlabel("presence_ratio")
    ax5.set_ylabel("number of clusters")

    ax6.hist(metrics["snr"], bins=20)
    ax6.set_xlabel("SNR")
    ax6.set_ylabel("number of clusters")

    ax7.hist(metrics["rp_contamination"], bins=20)
    ax7.set_xlabel("rp_contamination")

    ax8.hist(metrics["amplitude_cutoff"], bins=20)
    ax8.set_xlabel("amplitude_cutoff")

    fig.suptitle(f"{profile.session} --- {profile.metadata['probe_label'][profile.probe_id]}", fontsize=16, y=0.98)
    fig.text(0.5, 0.94, f"{profile.metadata['hardware_config'][profile.probe_id]}, {profile.metadata['probe_config'][profile.probe_id]}: depth = {profile.metadata['probe_depth_mm'][profile.probe_id]}mm, grid hole = {profile.metadata['probe_gridHole'][profile.probe_id]}", ha='center', fontsize=12)

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])  # leave space for labels
    fig.savefig(Path(profile.figs_path) / "analyzer_summary.png", dpi=300, bbox_inches="tight")

def plot_si_units(analyzer, profile, cluster_id): 
 
    fig = plt.figure(figsize=(12, 10))
    widget = UnitSummaryWidget(analyzer, unit_id=cluster_id, backend="matplotlib", figure=fig)

    fig.suptitle(f"{profile.session} --- {profile.metadata['probe_label'][profile.probe_id]}", fontsize=16, y=0.98)
    fig.text(0.5, 0.94, f"{profile.metadata['hardware_config'][profile.probe_id]}, {profile.metadata['probe_config'][profile.probe_id]}: depth = {profile.metadata['probe_depth_mm'][profile.probe_id]}mm, grid hole = {profile.metadata['probe_gridHole'][profile.probe_id]}", ha='center', fontsize=12)

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])  # leave space for labels
    fig.savefig(Path(profile.figs_path) / 'clusters' / f"clust{cluster_id:04d}_si.png", dpi=300, bbox_inches="tight")

def plot_ks_units(profile, cluster_id): 

    ops = np.load(profile.sorter_path / 'sorter_output' / 'ops.npy', allow_pickle=True).item()
    fname = ops['filename']
    if isinstance(fname, str) and 'PosixPath' in fname:
        import re
        fname = re.search(r"PosixPath\('(.+)'\)", fname).group(1)

        ops['filename'] = str(fname)
        np.save(profile.sorter_path / 'sorter_output' / 'ops.npy', ops)

    # Get the mean spike waveform and mean templates for the cluster
    mean_wv, spike_subset = mean_waveform(cluster_id, profile.sorter_path / 'sorter_output', n_spikes=500,
                                          bfile=None, best=True)
    mean_temp = cluster_templates(cluster_id, profile.sorter_path / 'sorter_output', mean=True,
                                  best=True, spike_subset=spike_subset)
    chan = get_best_channels(profile.sorter_path / 'sorter_output')[cluster_id]
    # Get n spike times for this cluster
    spike_times, _ = get_cluster_spikes(cluster_id, profile.sorter_path / 'sorter_output', n_spikes=500)
    # Time in s for spike time axis
    t2 = spike_times / ops['fs']
    t = (np.arange(ops['nt']) / ops['fs']) * 1000
    # Get single-channel waveform for each spike
    waves = get_spike_waveforms(spike_times, profile.sorter_path / 'sorter_output', chan=chan)


    fig = plt.figure(figsize=(12, 10))  # adjust as needed
    gs = GridSpec(10, 7, figure=fig)

    # Create axes with your specified layout
    ax0 = fig.add_subplot(gs[0:9, 0:1]) 
    ax1 = fig.add_subplot(gs[0:4, 1:6])
    ax2 = fig.add_subplot(gs[4:9, 1:6])


    im = ax0.imshow(waves.T, aspect='auto', extent=[t[0], t[-1], t2[0], t2[-1]])
    ax0.set_xlabel('Time (ms)')
    ax0.set_ylabel('Spike time (s)')
    cbar = plt.colorbar(im, ax=ax0, orientation='horizontal', pad=0.1)
    cbar.set_label('Amplitude (µV)')

    ax1.plot(t, mean_wv, c='black', linestyle='dashed', linewidth=2, label='waveform')
    ax1.plot(t, mean_temp, linewidth=1, label='template')
    ax1.set_title(f'unit_id: {cluster_id}')
    ax1.set_xlabel('Time (ms)')
    ax1.legend()

    n_waveforms = waves.shape[1]  # 500
    cmap = plt.cm.bwr  # blue-white-red
    colors = [cmap(i / n_waveforms) for i in range(n_waveforms)]

    # Plot each waveform with color gradient and transparency
    cmap = plt.cm.bwr  # blue → red
    norm = Normalize(vmin=0, vmax=n_waveforms-1)  # normalize waveform index

    # Plot each waveform with color gradient and transparency
    for i in range(n_waveforms):
        ax2.plot(t, waves[:, i], color=cmap(norm(i)), alpha=0.2)

    ax2.set_xlabel('Time (ms)')
    ax2.set_ylabel('Amplitude (µV)')
    ax2.set_title('Spike waveforms (blue → red, transparent)')

    # Add colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # needed for colorbar
    cbar = plt.colorbar(sm, ax=ax2, orientation='vertical', pad=0.01, shrink=0.85)
    cbar.set_label('Time in session (blue = start, red = end)')

    fig.suptitle(f"{profile.session} --- {profile.metadata['probe_label'][profile.probe_id]}", fontsize=16, y=0.98)
    fig.text(0.5, 0.94, f"{profile.metadata['hardware_config'][profile.probe_id]}, {profile.metadata['probe_config'][profile.probe_id]}: depth = {profile.metadata['probe_depth_mm'][profile.probe_id]}mm, grid hole = {profile.metadata['probe_gridHole'][profile.probe_id]}", ha='center', fontsize=12)

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])  # leave space for labels
    fig.savefig(Path(profile.figs_path) / 'clusters' / f"clust{cluster_id:04d}_ks.png", dpi=300, bbox_inches="tight")



