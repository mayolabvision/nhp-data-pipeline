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
from .path_utils import get_preprocess_hash, get_motion_hash, get_trim_hash, save_params, get_sorter_hash, get_sparse_hash
from .si_tools import run_motion_correction, save_processed_recording, detect_motion_cutoffs, get_mean_waveforms
from .si_tools import find_missing_extensions, add_extension_arrays_to_metrics, select_sparse_contacts
from .si_plots import plot_probe_motion, plot_preprocessing_steps, plot_probe_peaks, plot_noise_levels, plot_motion_correction_traces, plot_analyzer_sess, plot_si_units, plot_ks_units
from .ks_tools import convert_npy_to_mat, get_best_channels

from spikeinterface import aggregate_channels
from spikeinterface.sorters import run_sorter
from spikeinterface.core import set_global_job_kwargs, load, BaseRecording
from spikeinterface.preprocessing import apply_preprocessing_pipeline
from spikeinterface import create_sorting_analyzer, load_sorting_analyzer
from spikeinterface.extractors import get_neo_streams, read_spikeglx
from spikeinterface.qualitymetrics import compute_quality_metrics

cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
global_job_kwargs = dict(n_jobs=int(cpus - 1), chunk_duration='5s', progress_bar=True)
set_global_job_kwargs(**global_job_kwargs)

class NeuropixelProfile(RecordingProfile):
    def prep_session_data(self):
        self.data_path = Path(RAW_DATA_PATH) / self.session / f"{self.session}_{self.metadata['hardware_config'][self.probe_id]}"
        
        self.sparse_hash = get_sparse_hash(self.protocol["sparse"])    
      
        # preprocess/ preprocess hash / MC preprocess hash / MC hash
        if self.sparse_hash == "notsparse":
            self.preprocess_hash = get_preprocess_hash(self.protocol["preprocessing"])    
        else:
            preprocess_hash = get_preprocess_hash(self.protocol["preprocessing"])    
            self.preprocess_hash = preprocess_hash + self.sparse_hash
 
        self.pp_hash, self.motion_hash, self.pp_params, self.motion_params = get_motion_hash(self.protocol['motion_correction'])
        self.preprocess_path = self.data_path / "preprocess" / self.preprocess_hash / self.pp_hash / self.motion_hash
        save_params(self.preprocess_path.parent.parent / "params.json", self.protocol['preprocessing'])
        save_params(self.preprocess_path.parent / "params.json", self.pp_params)
       
        # shake_trimming / preprocess + shake hashes
        if self.protocol.get('shake_trimming'):
            ppShake_hash, moShake_hash, self.trim_hash, self.shake_params = get_trim_hash(self.protocol['shake_trimming'])    
            self.shake_hash = (ppShake_hash + moShake_hash)    
            self.shake_path = self.data_path / "shake_trimming" / self.shake_hash
            save_params(self.shake_path / "params.json", self.shake_params)
        else:
            self.trim_hash = "nocrop"
            self.shake_hash = "nocrop"

        # sorting / preprocess + MC preprocess + MC + shake + trim + sorter hashes
        self.sorter_hash, self.sorter_params, _ = get_sorter_hash(self.protocol['sorting'])

        self.full_hash = "-".join([self.preprocess_hash, self.pp_hash, self.motion_hash, self.shake_hash, self.trim_hash, self.sorter_hash])
        self.sorter_path = self.data_path / "sorting" / self.full_hash
 
        self.analyzer_path = self.sorter_path / 'analyzer'
        self.metrics_path = self.sorter_path / 'quality_metrics'
        
        self.tbl_path = self.data_path.parent / "tables" / f"{self.session}-{self.full_hash}.mat"

        self.figs_path = self.data_path.parent / "figs" / self.full_hash / f"{self.metadata['hardware_config'][self.probe_id]}_{self.metadata['probe_label'][self.probe_id]}"
        save_params(self.figs_path / "params.json", self.protocol)

    def preprocessing(self):
        print(f"===================================================================")
        # Check if data has already been preprocessed
        
        print(f"Loading in raw data for applying preprocessing......................")
        stream_names, stream_ids = get_neo_streams('spikeglx', self.data_path)
        recording = read_spikeglx(self.data_path, stream_name=f'imec{self.probe_id}.ap', load_sync_channel=False)

        if self.protocol.get('sparse'):
            keep_ids, remove_ids = select_sparse_contacts(recording, **self.protocol["sparse"])
            raw_recording = recording.remove_channels(remove_ids)
        else:
            raw_recording = recording
            
        print("--------------------- RAW RECORDING -----------------------------")
        print("Sampling frequency:", raw_recording.get_sampling_frequency())
        print("Number of channels:", raw_recording.get_num_channels())
        print("Number of segments:", raw_recording.get_num_segments())
        print("Number of samples:", raw_recording.get_num_samples(segment_index=0))
        print("Duration of recording (min):", round((raw_recording.get_num_samples(segment_index=0)/raw_recording.get_sampling_frequency())/60))
        print("Data dtype:", raw_recording.get_dtype())
        print("--------------------------------------------------")
        
        if not (self.preprocess_path / 'params.json').is_file():
            os.makedirs(self.preprocess_path, exist_ok=True)

            print(f"Applying preprocessing pipeline to raw data..........................")
            #split_recording = raw_recording.split_by("group")
           
            pp_recording = apply_preprocessing_pipeline(raw_recording, self.protocol["preprocessing"])
            #pp_recording = aggregate_channels(pp_recording)

            if self.protocol.get('motion_correction'):
                print(f"Preprocessing data with motion correction.............................")
                mc_recording = run_motion_correction(raw_recording, self.protocol, self.preprocess_path)
                save_processed_recording(mc_recording, self.preprocess_path / 'output')
                print(f"✓ Preprocessed data saved: {self.preprocess_path}")
                
                plot_probe_peaks(raw_recording, self)
                print(f"===== activity peaks on probe plotted =====")
                plot_motion_correction_traces(raw_recording, self)    
                print(f"===== motion-corrected traces plotted =====")
            
            else:
                print(f"Preprocessing data without motion correction..........................")
                save_processed_recording(pp_recording, self.preprocess_path / 'output')
                print(f"✓ Preprocessed data saved: {self.preprocess_path}")
                                    
            save_params(self.preprocess_path / "params.json", self.motion_params)
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
        if self.protocol.get('shake_trimming'):
            os.makedirs(self.shake_path, exist_ok=True)
            
            print(f"Loading in raw data for estimating probe motion......................")
            stream_names, stream_ids = get_neo_streams('spikeglx', self.data_path)
            raw_recording = read_spikeglx(self.data_path, stream_name=f'imec{self.probe_id}.ap', load_sync_channel=False)
            
            print(f"Estimating probe motion to crop out shaking...........................")
            crop_endSec = detect_motion_cutoffs(raw_recording, self)
            
            Fs = raw_recording.get_sampling_frequency()
            crop_recording = raw_recording.frame_slice(start_frame=0, end_frame=crop_endSec*Fs) 
            print("--------------------------------------------------")
            print("Duration of full recording (min):", round((raw_recording.get_num_samples(segment_index=0)/raw_recording.get_sampling_frequency())/60))
            print("Duration of crop recording (min):", round((crop_recording.get_num_samples(segment_index=0)/crop_recording.get_sampling_frequency())/60))
            
            print("--------------------------------------------------")
            
            print(f"===================================================================")
            print("Shake trimming of data compete!!!")
            print(f"===================================================================")
        
            convert_npy_to_mat(self.shake_path)
        
        else:
            print(f"===================================================================")
            print("Shake trimming is disabled, skipping this step")
            print(f"===================================================================")
        
        print(f"===================================================================")
        print("~~~~~~~~~~~~~PIPELINE STAGE 2: SHAKE TRIMMING COMPLETE~~~~~~~~~~~~~~")
        print(f"===================================================================")
        
    def spike_sorting(self):
        # Temporary deleting of sorting results to make changes
        #if self.sorter_path.exists():
        #    shutil.rmtree(self.sorter_path)

        print(f"Loading in preprocessed recording......................")
        recording = load(self.preprocess_path / 'output')
        _, _, custom_sorter_params = get_sorter_hash(self.protocol['sorting'])
        
        if self.protocol.get('shake_trimming'):
            print(f"Finding time to trim recording due to excessive motion.......................")
            self.crop_endSec = detect_motion_cutoffs(recording, self)
            custom_sorter_params['tmax'] = self.crop_endSec
            
            if not (Path(self.figs_path) / "probe_motion_cutoff.png").is_file(): 
                plot_probe_motion(self)
                print(f"===== probe motion w/ cutoff plotted =====")
     
        if custom_sorter_params.get('whitening_range') == 666 or custom_sorter_params.get('nearest_chans') == 666:
            prb = recording.get_probe().to_dataframe()
            prb = prb.sort_values('y')
            yp = np.diff(prb.y.values)
            ypitch = yp[yp>0].min()

            print(f"ypitch = {ypitch}")

            # set params so 110um spatial search 
            if ypitch==20:
                custom_sorter_params['nearest_chans'] = 16
                custom_sorter_params['whitening_range'] = 30
            elif ypitch==40:
                custom_sorter_params['nearest_chans'] = 8
                custom_sorter_params['whitening_range'] = 16
            else: # set to defaults
                custom_sorter_params['nearest_chans'] = 16
                custom_sorter_params['whitening_range'] = 16
   
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
        if self.protocol.get('shake_trimming'):
            self.crop_endSec = detect_motion_cutoffs(recording, self)
            Fs = recording.get_sampling_frequency()
            recording = recording.frame_slice(start_frame=0, end_frame=self.crop_endSec*Fs) 

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

        if not (Path(self.metrics_path) / "cluster_metrics.csv").is_file(): 
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
        if not (Path(self.figs_path) / "analyzer_summary.png").is_file(): 
            plot_analyzer_sess(analyzer, metrics, self) 
            print(f"===== analyzer session summary plotted =====")
        
        print(f"===================================================================")
        print("~~~~~~~~~~~~~PIPELINE STAGE 5: QUALITY METRICS CALCULATED~~~~~~~~~~~")
        print(f"===================================================================")


    def post_widgets(self, job_id=0, n_chunks=1):
        print("Loading cluster metrics.....................")
        metrics = pd.read_csv(Path(self.metrics_path) / "cluster_metrics.csv")      

        print(f"Loading in sorting analyzer......................")
        analyzer = load_sorting_analyzer(self.analyzer_path, 
                                         load_extensions=True, format="binary_folder")

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

