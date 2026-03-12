# recording_profiles/base.py
from abc import ABC, abstractmethod
from pathlib import Path
import json
import hashlib

from config import RAW_DATA_PATH

class RecordingProfile(ABC):
    def __init__(self, session: str, probe_id: int, protocol_path: str):
        self.session = session
        self.probe_id = probe_id
        self.protocol_path = protocol_path
        self.num_channels = None
        self.protocol = None
        self.data_path = None
        self.metadata = None
        self.probe_path = None
        self.preprocess_hash = None
        self.preprocess_path = None
        self.motion_hash = None
        self.motion_params = None
        self.sorter_path = None
        self.sorter_hash = None
        self.sorter_params = None
        self.full_hash = None
        self.analyzer_path = None
        self.metrics_path = None
        self.tbl_path = None
        self.figs_path = None
        self.crop_startSec = None
        self.crop_endSec = None

    def load_metadata(self):
        """Load metadata.json for this session."""
        metadata_path = Path(RAW_DATA_PATH) / self.session / "metadata.json"
        with open(metadata_path, "r") as f:
            self.metadata = json.load(f)

    def load_protocol(self):
        """Load protocol parameters."""
        with open(self.protocol_path, 'r') as f:
            self.protocol = json.load(f)

        if 'motion_crop' in self.protocol.get('preprocessing', {}) and self.protocol['preprocessing']['motion_crop'] is not True:
            del self.protocol['preprocessing']['motion_crop']

    @abstractmethod
    def prep_session_data(self):
        """Check and prepare raw data for processing."""
        return self

    @abstractmethod
    def preprocessing(self):
        """Applying preprocessing to raw recording."""
        pass
    
    @abstractmethod
    def spike_sorting(self):
        """Applying spike sorting to preprocessed recording."""
        pass
    
    @abstractmethod
    def postprocessing(self):
        """Applying postprocessing."""
        pass
    
    @abstractmethod
    def quality_metrics(self):
        """Calculate quality metrics from sorting_analyzer."""
        pass
    
