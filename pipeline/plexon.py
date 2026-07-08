import numpy as np
from pathlib import Path
import pandas as pd
from scipy.io import loadmat
import os
import json
import hashlib

import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

from .base import RecordingProfile
from config import RAW_DATA_PATH, PROBES_PATH
from .path_utils import get_preprocess_hash, save_params, get_sorter_hash
from .si_tools import save_processed_recording, get_mean_waveforms 
from .si_tools import find_missing_extensions, add_extension_arrays_to_metrics
from .si_plots import plot_preprocessing_steps, plot_probe_peaks, plot_noise_levels, plot_analyzer_sess, plot_si_units, plot_ks_units
from .ks_tools import convert_npy_to_mat, get_best_channels

from spikeinterface.sorters import run_sorter
from spikeinterface.core import set_global_job_kwargs, load, BaseRecording
from spikeinterface.preprocessing import apply_preprocessing_pipeline
from spikeinterface import create_sorting_analyzer, load_sorting_analyzer
from spikeinterface.extractors import read_binary
from spikeinterface.qualitymetrics import compute_quality_metrics

from .ripple_probe_maker import get_probe
from probeinterface.io import read_probeinterface

cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
global_job_kwargs = dict(n_jobs=int(cpus - 1), chunk_duration='5s', progress_bar=True)
set_global_job_kwargs(**global_job_kwargs)

class PlexonProfile(RecordingProfile):
    def prep_session_data(self):
        # Pull out raw data and save to sub-folder for each probe
        self.data_path = Path(RAW_DATA_PATH) / self.session / f"{self.session}_{self.metadata['hardware_config'][self.probe_id]}" 
        
        self.probe_path = self.data_path / "prbMap.json"
        if not (self.probe_path).is_file():
            get_probe(self.data_path, PROBES_PATH, probe_id=self.probe_id)
 
        self.preprocess_hash = get_preprocess_hash(self.protocol["preprocessing"])    
        self.preprocess_path = self.data_path / "preprocess" / self.preprocess_hash 
        
        self.sorter_hash, self.sorter_params, _ = get_sorter_hash(self.protocol['sorting'])
        self.full_hash = "-".join([self.preprocess_hash, self.sorter_hash])
        self.sorter_path = self.data_path / "sorting" / self.full_hash

        # Get only the part after "sorting"
        short_sorter_path = str(self.sorter_path).split("sorting", 1)[1].lstrip("/\\")
        self.short_hash = hashlib.sha256(short_sorter_path.encode()).hexdigest()[:16]

        self.analyzer_path = self.sorter_path / 'analyzer'
        self.metrics_path = self.sorter_path / 'quality_metrics'

        self.tbl_path = self.data_path.parent / "tables" / f"{self.session}-{self.short_hash}.mat"
        self.figs_path = self.data_path.parent / "figs" / self.short_hash / f"{self.metadata['hardware_config'][self.probe_id]}_{self.metadata['probe_label'][self.probe_id]}"
        save_params(self.figs_path / "params.json", self.protocol)
    
    def preprocessing(self):
        print(f"===================================================================")

        print(f"Loading in raw data for applying preprocessing......................")
        with open(self.data_path / "ripple_info.json", "r") as f:
            ripple_info = json.load(f)

        raw_recording = read_binary(file_paths=self.data_path / "raw_signal.bin",
                                    sampling_frequency=ripple_info["Fs"],
                                    num_channels=ripple_info["num_channels"],
                                    dtype=ripple_info["dtype_python"],
                                    gain_to_uV=ripple_info.get("gain_to_uV"),
                                    offset_to_uV=ripple_info.get("offset_to_uV"))

        print("--------------------- RAW RECORDING -----------------------------")
        print("Sampling frequency:", raw_recording.get_sampling_frequency())
        print("Number of channels:", raw_recording.get_num_channels())
        print("Number of segments:", raw_recording.get_num_segments())
        print("Number of samples:", raw_recording.get_num_samples(segment_index=0))
        print("Duration of recording (min):", round((raw_recording.get_num_samples(segment_index=0)/raw_recording.get_sampling_frequency())/60))
        print("Data dtype:", raw_recording.get_dtype())
        print("--------------------------------------------------")

        prb = read_probeinterface(self.probe_path)
        raw_recording = raw_recording._set_probes(prb)

        if not (self.preprocess_path / 'params.json').is_file():
            os.makedirs(self.preprocess_path, exist_ok=True)

            print(f"Applying preprocessing pipeline to raw data..........................")
            pp_recording = apply_preprocessing_pipeline(raw_recording, self.protocol["preprocessing"]) 

            print(f"Preprocessing data without motion correction..........................")
            save_processed_recording(pp_recording, self.preprocess_path / 'output')
            print(f"✓ Preprocessed data saved: {self.preprocess_path}")

            save_params(self.preprocess_path / "params.json", self.protocol["preprocessing"])
            convert_npy_to_mat(self.preprocess_path)

            print(f"===================================================================")
            print("Preprocessing of data compete!!!")
            print(f"===================================================================")

        else:
            print(f"===================================================================")
            print("Preprocessed data already exists, skipping this step")
            print(f"===================================================================")


        print(f"Plotting preprocessing steps..........................")
        if not (Path(self.figs_path) / "noise_levels.png").is_file():
            plot_noise_levels(raw_recording, self)
            print(f"===== distributions of noise_levels plotted =====")

        if not (Path(self.figs_path) / "preprocessing_steps.png").is_file():
            plot_preprocessing_steps(raw_recording, self)
            print(f"===== traces for preprocessing_steps plotted =====")

        print(f"===================================================================")
        print("~~~~~~~~~~~~~PIPELINE STAGE 1: PREPROCESSING COMPLETE~~~~~~~~~~~~~~~")
        print(f"===================================================================")
    
    def shake_trimming(self):
        pass

    def spike_sorting(self):
    
        print(f"Loading in preprocessed recording......................")
        recording = load(self.preprocess_path / 'output')
        _, _, custom_sorter_params = get_sorter_hash(self.protocol['sorting'])

        print("--------------------------------------------------")
        print(custom_sorter_params)
        print("--------------------------------------------------")
    
        if not (self.sorter_path / 'params.json').is_file():

            sorting = run_sorter(recording=recording, folder=self.sorter_path,
                                 verbose=True, save_preprocessed_copy=False, docker_image=False,
                                 clear_cache=True, remove_existing_folder=True, **custom_sorter_params)

            save_params(self.data_path / "sorting" / self.full_hash / "params.json",
                            self.protocol)
            convert_npy_to_mat(self.sorter_path / 'sorter_output')
            print(f"✓ Spike sorting outputs saved: {self.sorter_path}")

        else:
            print("Sorted outputs already exist, skipping spike sorting")

        print(f"===================================================================")
        print("~~~~~~~~~~~~~~PIPELINE STAGE 3: SPIKE SORTING COMPLETE~~~~~~~~~~~~~~")
        print(f"===================================================================")


    def postprocessing(self):
        print(f"Loading in preprocessed recording......................")
        recording = load(self.preprocess_path / 'output')

        print(f"Loading in sorting outputs......................")
        sorting = load(self.sorter_path)

        print(f"Postprocessing sorted recording......................")
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
            cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
            job_kwargs = dict(n_jobs=int(cpus - 1), chunk_duration='1s', progress_bar=True)

            analyzer.compute(extensions_to_run, **job_kwargs)

            save_params(self.analyzer_path / "params.json", self.protocol['postprocessing'])
            print("✓ Extensions ran:", extensions_to_run)
        else:
            print("all postprocessing extensions have already been ran")

        print(f"===================================================================")
        print("~~~~~~~~~~~~~~PIPELINE STAGE 4: POSTPROCESSING COMPLETE~~~~~~~~~~~~~")
        print(f"===================================================================")


    def quality_metrics(self):
        print(f"Loading in sorting analyzer......................")
        analyzer = load_sorting_analyzer(self.analyzer_path)

        if not (Path(self.metrics_path) / "cluster_metrics666.csv").is_file():
            print(f"Computing quality metrics......................")
            metrics = compute_quality_metrics(
                    analyzer,
                    metric_names=list(self.protocol['quality_metrics'].keys()),
                    metric_params=self.protocol['quality_metrics'])

            metrics["sess_name"] = self.metadata["sess_name"]
            metrics["monkey"] = self.metadata["monkey"]
            metrics["experimenter"] = self.metadata["experimenter"]
            metrics["probe_id"] = self.probe_id
            metrics["cluster_id"] = metrics.index
            for key in ["probe_label", "probe_type", "probe_config",
                        "hardware_config", "probe_depth_mm", "probe_gridHole"]:
                value = self.metadata[key][self.probe_id]
                metrics[key] = [value] * len(metrics)

            # Specify the order you want
            if self.protocol['sorting']['sorter_name'] == 'kilosort4':
                metrics["best_channel"] = get_best_channels(self.sorter_path / 'sorter_output')
                new_cols = ["sess_name","monkey","experimenter","probe_id","cluster_id","best_channel","probe_label","probe_type","probe_config",
                            "hardware_config","probe_depth_mm","probe_gridHole"]
            else:
                new_cols = ["sess_name","monkey","experimenter","probe_id","cluster_id","probe_label","probe_type","probe_config",
                            "hardware_config","probe_depth_mm","probe_gridHole"]
            metrics = metrics[new_cols + [col for col in metrics.columns if col not in new_cols]]
            metrics = add_extension_arrays_to_metrics(self.analyzer_path / "extensions", metrics)

            mean_wfs = get_mean_waveforms(analyzer)
            metrics['mean_waveform'] = [
                mean_wfs.get(unit_id, np.array([]))
                for unit_id in metrics['cluster_id']
            ]
            print("✓ Waveforms for each unit pulled.")

            save_params(self.metrics_path / "params.json", self.protocol['quality_metrics'])
            metrics.to_csv(self.metrics_path / 'cluster_metrics.csv', index=False)
            print("✓ Quality metrics calculated.")
        else:
            metrics = pd.read_csv(Path(self.metrics_path) / "cluster_metrics.csv")
            print("cluster_metrics already exists, loading in df")

        print(f"Plotting analyzer summary and metrics......................")
        if not (Path(self.figs_path) / "analyzer_summary666.png").is_file():
            plot_analyzer_sess(analyzer, metrics, self)
            print(f"===== analyzer session summary plotted =====")

        print(f"===================================================================")
        print("~~~~~~~~~~~~~PIPELINE STAGE 5: QUALITY METRICS CALCULATED~~~~~~~~~~~")
        print(f"===================================================================")


    def post_widgets(self, job_id=0, n_chunks=1):
        print("Loading cluster metrics.....................")
        metrics = pd.read_csv(Path(self.metrics_path) / "cluster_metrics.csv")

        print(f"Loading in sorting analyzer......................")
        analyzer = load_sorting_analyzer(self.analyzer_path)

        os.makedirs(Path(self.figs_path) / "clusters", exist_ok=True)

        # ---- Chunking logic ----
        n_rows = len(metrics)
        chunk_size = int(np.ceil(n_rows / n_chunks))

        start_idx = job_id * chunk_size
        end_idx = min((job_id + 1) * chunk_size, n_rows)

        print(f"Job {job_id}/{n_chunks} processing rows {start_idx}:{end_idx}")
        metrics_chunk = metrics.iloc[start_idx:end_idx]

        # ---- Loop only over this chunk ----
        for i, (idx, row) in enumerate(metrics_chunk.iterrows(), start=1):
            cluster_id = row['cluster_id']
            print('.', end='', flush=True)

            out_file = Path(self.figs_path) / "clusters" / f"clust{cluster_id:04d}_si.png"

            if not out_file.is_file():
                plot_si_units(analyzer, self, cluster_id)
                plot_ks_units(self, cluster_id)

        print(f"\n===== analyzer per units plotted (job {job_id}) =====")
