import numpy as np
from pathlib import Path
import pandas as pd
from scipy.io import loadmat
import os
import json

import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

from .base import RecordingProfile
from config import RAW_DATA_PATH
from .path_utils import get_preprocess_hash, get_motion_hash, save_params, get_sorter_hash
from .si_tools import run_preprocessing_with_motion_correction, run_preprocessing_without_motion_correction, save_processed_recording, detect_excessive_motion
from .si_tools import find_missing_extensions, add_extension_arrays_to_metrics
from .si_plots import plot_motion_screening, plot_preprocessing_steps, plot_probe_peaks, plot_noise_levels, plot_motion_correction_traces
from .ks_tools import convert_npy_to_mat, get_best_channels

from spikeinterface.sorters import run_sorter
from spikeinterface.core import load, BaseRecording
from spikeinterface import create_sorting_analyzer, load_sorting_analyzer
from spikeinterface.extractors import read_binary
from spikeinterface.qualitymetrics import compute_quality_metrics

from probeinterface.generator import generate_linear_probe

class PlexonProfile(RecordingProfile):
    def prep_session_data(self):
        # Pull out raw data and save to sub-folder for each probe
        self.data_path = Path(RAW_DATA_PATH) / self.session / f"{self.metadata['hardware_config'][self.probe_id]}_{self.metadata['probe_label'][self.probe_id]}"
        self.probe_path = self.data_path / "prb.mat"
       
        prb = loadmat(self.probe_path)
        self.num_channels = int(len(prb['chanMap']))
 
        self.preprocess_hash = get_preprocess_hash(self.protocol["motion_screening"] | self.protocol["preprocessing"])    
        self.pp_hash, self.motion_hash, self.pp_params, self.motion_params = get_motion_hash(self.protocol['motion_correction'])
        self.preprocess_path = self.data_path / "preprocess" / self.preprocess_hash / self.pp_hash / self.motion_hash
        save_params(self.preprocess_path.parent.parent / "params.json", self.protocol['preprocessing'])
        save_params(self.preprocess_path.parent / "params.json", self.pp_params)
        
        self.sorter_hash, self.sorter_params, _ = get_sorter_hash(self.protocol['sorting'])
        self.full_hash = "-".join([self.preprocess_hash, self.motion_hash, self.sorter_hash])
        self.sorter_path = self.data_path / "sorting" / self.full_hash
 
        self.analyzer_path = self.sorter_path / 'analyzer'
        self.metrics_path = self.sorter_path / 'quality_metrics'
        
        self.tbl_path = self.data_path.parent / "tables" / f"{self.session}-{self.full_hash}.mat"
        self.figs_path = self.data_path.parent / "figs" / self.full_hash / f"{self.metadata['hardware_config'][self.probe_id]}_{self.metadata['probe_label'][self.probe_id]}"
        save_params(self.figs_path / "params.json", self.protocol)
    
    def motion_screening(self):
        if self.protocol["motion_screening"]["enabled"]:
            raw_recording = read_binary(file_paths=self.data_path / "raw.bin", num_channels=self.num_channels, 
                                        sampling_frequency=30000, dtype="float64", time_axis=1)
            
            cutoff_time_sec = detect_excessive_motion(raw_recording, self.preprocess_path, 
                                                      threshold_um=self.protocol['motion_screening']['motion_thresh_um'],
                                                      min_duration_sec=self.protocol['motion_screening']['min_duration_sec'])
            self.cutoff_frame = None if cutoff_time_sec is None else int(cutoff_time_sec * raw_recording.get_sampling_frequency())

            if not (self.figs_path / "motion_screening.png").is_file():             
                plot_motion_screening(self, cutoff_time_sec)
                print(f"===== motion screening plotted =====")
        else:
            self.cutoff_frame = None

    def preprocessing(self):
        if not (self.preprocess_path / 'params.json').is_file(): 
            raw_recording = read_binary(file_paths=self.data_path / "raw.bin", num_channels=self.num_channels, 
                                        sampling_frequency=30000, dtype="float64", time_axis=1)
           
            print("Duration of raw recording (min):", round((raw_recording.get_num_samples(segment_index=0)/raw_recording.get_sampling_frequency())/60))
            trim_recording = raw_recording.frame_slice(start_frame=0, end_frame=self.cutoff_frame)
            print("Duration of trimmed recording (min):", round((trim_recording.get_num_samples(segment_index=0)/trim_recording.get_sampling_frequency())/60))
     
            plot_noise_levels(trim_recording,self)
            print(f"===== distributions of noise_levels plotted =====")
            plot_preprocessing_steps(trim_recording,self)
            print(f"===== traces for preprocessing_steps plotted =====")
            
            if self.protocol.get('motion_correction'):
                mc_recording = run_preprocessing_with_motion_correction(trim_recording, self.protocol, self.preprocess_path)
                print(f"===== preprocessing with motion correction complete =====")
                
                plot_probe_peaks(trim_recording,self)
                print(f"===== activity peaks on probe plotted =====")
                plot_motion_correction_traces(trim_recording,self)    
                print(f"===== motion-corrected traces plotted =====")
                
            else:
                mc_recording = run_preprocessing_without_motion_correction(trim_recording, self.protocol, self.preprocess_path)
                
            mc_recording = save_processed_recording(mc_recording, self.preprocess_path / 'output')
            save_params(self.preprocess_path / "params.json", self.motion_params)
            
            print(f"✓ Preprocessed data saved: {self.preprocess_path}") 
        else:
            print("Preprocessed data already exists, skipping save_preprocessed_recording")

    def spike_sorting(self):
        if not (self.sorter_path / 'params.json').is_file():
            recording = load(self.preprocess_path / 'output')

            _, _, custom_sorter_params = get_sorter_hash(self.protocol['sorting'])
            sorting = run_sorter(recording=recording, folder=self.sorter_path, 
                                 verbose=True, save_preprocessed_copy=False, docker_image=False,
                                 clear_cache=True, **custom_sorter_params)
            
            save_params(self.data_path / "sorting" / self.full_hash / "params.json", 
                            self.protocol)
            print(f"✓ Spike sorting outputs saved: {self.sorter_path}")
        else:
            print("Sorted outputs already exist, skipping spike sorting")

        convert_npy_to_mat(self.sorter_path / 'sorter_output')

    def postprocessing(self):
        recording = load(self.preprocess_path / 'output')
        sorting = load(self.sorter_path)

        if not self.analyzer_path.is_dir():
            analyzer = create_sorting_analyzer(sorting=sorting, recording=recording, format="binary_folder", sparse=True,
                                               return_in_uV=True, folder=self.analyzer_path)
            print(f"✓ Postprocessing sorting_analyzer saved: {self.analyzer_path}")
        else:
            analyzer = load_sorting_analyzer(self.analyzer_path)
            print("sorting_analyzer already exists, loading in extensions")

        extensions_dir = self.analyzer_path / "extensions"
        extensions_dir.mkdir(parents=True, exist_ok=True)
        extensions_to_run = find_missing_extensions(extensions_dir, self.protocol["postprocessing"]) 
        if extensions_to_run:
            analyzer.compute(extensions_to_run)
            
            save_params(self.analyzer_path / "params.json", self.protocol['postprocessing'])
            print("✓ Extensions ran:", extensions_to_run)
        else:
            print("all postprocessing extensions have already been ran")

    def quality_metrics(self):
        analyzer = load_sorting_analyzer(self.analyzer_path)

        #if not self.metrics_path.is_file():
        metrics = compute_quality_metrics(
                analyzer,
                metric_names=list(self.protocol['quality_metrics'].keys()),
                metric_params=self.protocol['quality_metrics'])
       
        metrics["sess_name"] = self.metadata["sess_name"]
        metrics["probe_id"] = self.probe_id
        metrics["cluster_id"] = metrics.index
        for key in ["probe_label", "probe_type", "probe_config",
                    "hardware_config", "probe_depth_mm", "probe_gridHole"]:
            value = self.metadata[key][self.probe_id]
            metrics[key] = [value] * len(metrics)        

        # Specify the order you want
        if self.protocol['sorting']['sorter_name'] == 'kilosort4':
            metrics["best_channel"] = get_best_channels(self.sorter_path / 'sorter_output')
            new_cols = ["sess_name","probe_id","cluster_id","best_channel","probe_label","probe_type","probe_config",
                        "hardware_config","probe_depth_mm","probe_gridHole"]
        else:
            new_cols = ["sess_name","probe_id","cluster_id","probe_label","probe_type","probe_config",
                        "hardware_config","probe_depth_mm","probe_gridHole"]

        metrics = metrics[new_cols + [col for col in metrics.columns if col not in new_cols]]
        metrics = add_extension_arrays_to_metrics(self.analyzer_path / "extensions", metrics)

        save_params(self.metrics_path / "params.json", self.protocol['quality_metrics'])
        metrics.to_csv(self.metrics_path / 'cluster_metrics.csv', index=False)
        print("✓ Quality metrics calculated.")

