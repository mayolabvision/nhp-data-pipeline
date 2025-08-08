import numpy as np
from pathlib import Path
import os
import json

from .base import RecordingProfile
from config import RAW_DATA_PATH
from .catgt_utils import run_catgt
from .SGLXMetaToCoords import MetaToCoords
from .path_utils import get_preprocess_hash, get_motion_hash, save_params
from .si_preprocess import run_preprocessing_with_motion_correction, run_preprocessing_without_motion_correction, save_processed_recording

class NeuropixelProfile(RecordingProfile):
    def prep_session_data(self):
        # 1. Set raw data path
        self.data_path = Path(RAW_DATA_PATH) / self.session / f"{self.session}_imec{self.probe_id}"
        if not self.data_path.is_dir():
            raise FileNotFoundError(f"Raw data directory not found: {self.data_path}")
        print(f"✓ Raw data found: {self.data_path}")

        # 2. Run CatGT if needed, pulls out imec sync pulses
        run_catgt(self.session, Path(RAW_DATA_PATH))

    def make_probe_map(self):
        self.probe_path = self.data_path / f"{self.session}_t0.imec{self.probe_id}.ap_kilosortChanMap.mat"

        if not self.probe_path.exists():
            MetaToCoords(self.probe_path.parent / f"{self.session}_t0.imec{self.probe_id}.ap.meta",
                         1, badChan=np.zeros((0), dtype='int'), destFullPath='', showPlot=True)
            print(f"✓ Probe map made: {self.probe_path}")
        else:
            print(f"Probe map already exists, skipping make_probe_map")

    def preprocessing(self):
        self.preprocess_hash = get_preprocess_hash(self.protocol['preprocessing'])    
        save_params(self.data_path / "preprocess" / self.preprocess_hash / "params.json",
            self.protocol['preprocessing'])
 
        if 'motion_correction' in self.protocol and self.protocol['motion_correction']:
            self.pp_hash, self.motion_hash, self.motion_params = get_motion_hash(self.protocol['motion_correction'])
            save_params(self.data_path / "preprocess" / self.preprocess_hash / self.pp_hash / "params.json",
                self.protocol['motion_correction']['preprocessing'])
            
            save_params(self.data_path / "preprocess" / self.preprocess_hash / self.pp_hash / self.motion_hash / "params.json",
                self.motion_params)
            
            self.preprocess_path = self.data_path / "preprocess" / self.preprocess_hash / self.pp_hash / self.motion_hash
            mc_recording = run_preprocessing_with_motion_correction(self.data_path, self.probe_id, self.protocol, self.preprocess_path)
            print(f"✓ Preprocessing w/ drift correction complete...")

        else:
            self.preprocess_path = self.data_path / "preprocess" / self.preprocess_hash / "nodrift"
            mc_recording = run_preprocessing_without_motion_correction(self.data_path, self.probe_id, self.protocol, self.preprocess_path)
            print(f"✓ Preprocessing w/out drift correction complete...")

        if not (self.preprocess_path / 'output' / 'traces_cached_seg0.raw').is_file():
            mc_recording = save_processed_recording(mc_recording, self.preprocess_path / 'output')
            print(f"✓ Preprocessed data saved: {self.preprocess_path}")
        else:
            print(f"Preprocessed data already exists, skipping save_preprocessed_recording")
            

