# main_pipeline.py
from pathlib import Path

from pipeline import get_recording_profile
from config import PROTOCOLS_PATH

def run_preprocess(session, probe_id=0, protocol=None):
    # Get appropriate recording profile class, based on probe type
    profile_cls = get_recording_profile(session, probe_id)
    profile = profile_cls(session, probe_id, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()

    # Check data exists and is formatted correctly
    profile.prep_session_data()

    # Make probe map
    profile.make_probe_map()

    profile.preprocessing()

def run_sorting(session, probe_id=0, protocol=None):
    # Get appropriate recording profile class, based on probe type
    profile_cls = get_recording_profile(session, probe_id)
    profile = profile_cls(session, probe_id, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()

    # Check data exists and is formatted correctly
    profile.prep_session_data()

    # Make probe map
    profile.make_probe_map()

    profile.preprocessing()

