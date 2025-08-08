# recording_profiles/base.py
from abc import ABC, abstractmethod
from pathlib import Path
import json
import hashlib

from config import RAW_DATA_PATH

class RecordingProfile(ABC):
    def __init__(self, session: str, probe_id: str, protocol_path: str):
        self.session = session
        self.probe_id = probe_id
        self.protocol_path = protocol_path
        self.protocol = None
        self.data_path = None
        self.metadata_path = None
        self.metadata = None
        self.probe_path = None
        self.preprocess_hash = None
        self.preprocess_path = None
        self.motion_hash = None
        self.motion_params = None

    def load_metadata(self):
        """Load metadata.json for this session."""
        self.metadata_path = Path(RAW_DATA_PATH) / self.session / "metadata.json"
        with open(self.metadata_path, "r") as f:
            self.metadata = json.load(f)

    def load_protocol(self):
        """Load protocol parameters."""
        with open(self.protocol_path, 'r') as f:
            self.protocol = json.load(f)

    @abstractmethod
    def prep_session_data(self):
        """Check and prepare raw data for processing."""
        pass

    @abstractmethod
    def make_probe_map(self):
        """Check and prepare raw data for processing."""
        pass

    @abstractmethod
    def preprocessing(self):
        """Applying preprocessing to raw recording."""
        pass
 
