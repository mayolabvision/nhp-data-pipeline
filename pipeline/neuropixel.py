import numpy as np
from pathlib import Path
import pandas as pd
import os
import shutil
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import json
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

from .base import RecordingProfile
from config import RAW_DATA_PATH
from .path_utils import get_preprocess_hash, get_motion_hash, save_params, get_sorter_hash
from .si_tools import run_preprocessing_with_motion_correction, run_preprocessing_without_motion_correction, save_processed_recording, detect_probe_motion
from .si_tools import find_missing_extensions, add_extension_arrays_to_metrics
from .si_plots import plot_probe_motion, plot_preprocessing_steps, plot_probe_peaks, plot_noise_levels, plot_motion_correction_traces
from .ks_tools import convert_npy_to_mat, get_best_channels

from spikeinterface.sorters import run_sorter
from spikeinterface.core import load, BaseRecording
from spikeinterface import create_sorting_analyzer, load_sorting_analyzer
from spikeinterface.extractors import get_neo_streams, read_spikeglx
from spikeinterface.qualitymetrics import compute_quality_metrics

class NeuropixelProfile(RecordingProfile):
    def prep_session_data(self):
        self.data_path = Path(RAW_DATA_PATH) / self.session / f"{self.session}_{self.metadata['hardware_config'][self.probe_id]}"
        
        self.preprocess_hash = get_preprocess_hash(self.protocol["preprocessing"])    
        self.pp_hash, self.motion_hash, self.pp_params, self.motion_params = get_motion_hash(self.protocol['motion_correction'])
        self.preprocess_path = self.data_path / "preprocess" / self.preprocess_hash / self.pp_hash / self.motion_hash
        save_params(self.preprocess_path.parent.parent / "params.json", self.protocol['preprocessing'])
        save_params(self.preprocess_path.parent / "params.json", self.pp_params)
        
        self.sorter_hash, self.sorter_params, _ = get_sorter_hash(self.protocol['sorting'])
        self.full_hash = "-".join([self.preprocess_hash, self.pp_hash, self.motion_hash, self.sorter_hash])
        self.sorter_path = self.data_path / "sorting" / self.full_hash
 
        self.analyzer_path = self.sorter_path / 'analyzer'
        self.metrics_path = self.sorter_path / 'quality_metrics'
        
        self.tbl_path = self.data_path.parent / "tables" / f"{self.session}-{self.full_hash}.mat"

        self.figs_path = self.data_path.parent / "figs" / self.full_hash / f"{self.metadata['hardware_config'][self.probe_id]}_{self.metadata['probe_label'][self.probe_id]}"
        save_params(self.figs_path / "params.json", self.protocol)

    def preprocessing(self):
        # Check if data has already been preprocessed
        if not (self.preprocess_path / 'params.json').is_file():
            print(f"--- Loading in raw data for applying preprocessing ---")
            stream_names, stream_ids = get_neo_streams('spikeglx', self.data_path)
            raw_recording = read_spikeglx(self.data_path, stream_name=f'imec{self.probe_id}.ap', load_sync_channel=False)
            print("--------------------------------------------------")
            print("Sampling frequency:", raw_recording.get_sampling_frequency())
            print("Number of channels:", raw_recording.get_num_channels())
            print("Number of segments:", raw_recording.get_num_segments())
            print("Number of samples:", raw_recording.get_num_samples(segment_index=0))
            print("Duration of recording (min):", round((raw_recording.get_num_samples(segment_index=0)/raw_recording.get_sampling_frequency())/60))
            print("Data dtype:", raw_recording.get_dtype())
            print("--------------------------------------------------")
            
            plot_noise_levels(raw_recording,self)
            print(f"===== distributions of noise_levels plotted =====")
            plot_preprocessing_steps(raw_recording,self)
            print(f"===== traces for preprocessing_steps plotted =====")

            if self.protocol.get('motion_correction'):
                print(f"--- Preprocessing data with motion correction ---")
                mc_recording = run_preprocessing_with_motion_correction(raw_recording, self.protocol, self.preprocess_path)
                save_processed_recording(mc_recording, self.preprocess_path / 'output')
                print(f"✓ Preprocessed data saved: {self.preprocess_path}")
                
                plot_probe_peaks(raw_recording,self)
                print(f"===== activity peaks on probe plotted =====")
                plot_motion_correction_traces(raw_recording,self)    
                print(f"===== motion-corrected traces plotted =====")
            
            else:
                print(f"--- Preprocessing data without motion correction ---")
                if (self.preprocess_path / 'output').exists(): #and (self.preprocess_path / 'motion.npy').is_file():
                    mc_recording = load(self.preprocess_path / 'output')
                else:
                    mc_recording = run_preprocessing_without_motion_correction(raw_recording, self.protocol, self.preprocess_path)
                    save_processed_recording(mc_recording, self.preprocess_path / 'output')
                    print(f"✓ Preprocessed data saved: {self.preprocess_path}")
                                    
                print(f"Estimating probe motion...........................")
                if not (self.preprocess_path / 'motion.npy').is_file():
                    detect_probe_motion(mc_recording, self.preprocess_path)
                
                plot_probe_motion(self)
                print(f"===== probe motion estimated and plotted =====")
                
            save_params(self.preprocess_path / "params.json", self.motion_params)
            
            print("Preprocessing of data compete!!!")
        else:
            print("Preprocessed data already exists, skipping this step")
        
        convert_npy_to_mat(self.preprocess_path)
    
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

