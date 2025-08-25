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
    profile_cls = get_recording_profile(session, probe_id)
    profile = profile_cls(session, probe_id, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()
    profile.prep_session_data()
    
    # Run sorting
    profile.spike_sorting()

def run_postprocess(session, probe_id=0, protocol=None):
    profile_cls = get_recording_profile(session, probe_id)
    profile = profile_cls(session, probe_id, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()
    profile.prep_session_data()
    
    profile.postprocessing()
    profile.quality_metrics()

def profile_to_mat(session, protocol=None):
    profile_cls = get_recording_profile(session, 0)
    profile = profile_cls(session, 0, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()
    profile.prep_session_data()

    # Get only the part after "sorting"
    sorter_path = str(profile.sorter_path).split("sorting", 1)[1].lstrip("/\\")

    return sorter_path
