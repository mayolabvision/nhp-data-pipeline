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

class BehaviorProfile(RecordingProfile):
    def prep_session_data(self):
        # Pull out raw data and save to sub-folder for each probe
        self.data_path = Path(RAW_DATA_PATH) / self.session
        
        self.tbl_path = self.data_path / "tables" / f"{self.session}.mat"
        self.figs_path = self.data_path / "figs" 
    
    def preprocessing(self):
        pass

    def spike_sorting(self):
        pass

    def postprocessing(self):
        pass

    def quality_metrics(self):
        pass
