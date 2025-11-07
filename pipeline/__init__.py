import json
from pathlib import Path
from datetime import date
from .behavior import BehaviorProfile
from .fhc import FHCProfile
from .neuropixel import NeuropixelProfile
from .plexon import PlexonProfile
from config import RAW_DATA_PATH, PROTOCOLS_PATH

# Mapping from probe_type to class
RECORDING_PROFILE_MAP = {
    "neuropixel": NeuropixelProfile,
    "plexon": PlexonProfile,
    "fhc": FHCProfile,
    "behavior": BehaviorProfile
}

def get_recording_profile(session, probe_id):
    """Load metadata.json for this session."""
    metadata_path = Path(RAW_DATA_PATH) / session / "metadata.json"
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found at {metadata_path}")

    with open(metadata_path, "r") as f:
        metadata = json.load(f)

    if not metadata.get("probe_type"):   # empty list, None, or missing key
        probe_type = "behavior"
    else:
        probe_type = metadata["probe_type"][probe_id]

    try:
        return RECORDING_PROFILE_MAP[probe_type.lower()]
    except KeyError:
        raise ValueError(f"No recording profile found for probe type '{probe_type}'")

def save_protocol_to_dict(protocol_name, full_hash, session):
    map_path = RAW_DATA_PATH / "protocol_map.json"
    today = date.today().isoformat()  # YYYY-MM-DD format

    if map_path.exists():
        with open(map_path, "r") as f:
            protocol_map = json.load(f)
    else:
        protocol_map = {}

    # If hash is already present, update sessions only
    if full_hash in protocol_map:
        # Check if session already exists
        existing_sessions = [s["session"] for s in protocol_map[full_hash]["sessions"]]
        if session not in existing_sessions:
            protocol_map[full_hash]["sessions"].append({"session": session, "date": today})

    else:
        # Load the protocol params from protocol file
        with open(Path(PROTOCOLS_PATH) / protocol_name, "r") as f:
            protocol_params = json.load(f)

        # Add new entry
        protocol_map[full_hash] = {
            "protocol_name": protocol_name,
            "sessions": [
                {"session": session, "date": today}
            ],
            "protocol_params": protocol_params
        }

    # Save back to file
    with open(map_path, "w") as f:
        json.dump(protocol_map, f, indent=4)
