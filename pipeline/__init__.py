import json
from pathlib import Path
from .neuropixel import NeuropixelProfile
from config import RAW_DATA_PATH

# Mapping from probe_type to class
RECORDING_PROFILE_MAP = {
    "neuropixel": NeuropixelProfile,
    # add others as needed
}

def get_recording_profile(session, probe_id):
    """Load metadata.json for this session."""
    metadata_path = Path(RAW_DATA_PATH) / session / "metadata.json"
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found at {metadata_path}")

    with open(metadata_path, "r") as f:
        metadata = json.load(f)

    probe_type = metadata["probe_type"][probe_id]
    #print(f"Detected probe type: {probe_type}")

    try:
        return RECORDING_PROFILE_MAP[probe_type.lower()]
    except KeyError:
        raise ValueError(f"No recording profile found for probe type '{probe_type}'")


