# main_pipeline.py
from pathlib import Path
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"

from pipeline import get_recording_profile, save_protocol_to_dict
from config import PROTOCOLS_PATH

def run_preprocess(session, probe_id=0, protocol=None):
    profile_cls = get_recording_profile(session, probe_id)
    profile = profile_cls(session, probe_id, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()
    print(f"---------- ✓ Protocol and metadata loaded ----------")
    print(f"metadata         =  {profile.metadata}")
    
    profile.prep_session_data()
    print(f"---------- ✓ Paths established with hash ----------")
    print(f"data_path        =  {profile.data_path}")
    print(f"figs_path        =  {profile.figs_path}")
    print(f"preprocess_hash  =  {profile.preprocess_hash}")
    print(f"motion_hash      =  {profile.motion_hash}")
    
    profile.preprocessing()
    print(f"------------ ✓ Preprocessing complete -------------")
    
    return profile

def run_shakeTrimming(session, probe_id=0, protocol=None):
    profile_cls = get_recording_profile(session, probe_id)
    profile = profile_cls(session, probe_id, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()
    profile.prep_session_data()
    
    # Estimate probe drift to detect high motion
    profile.shake_trimming()
    print(f"------------ ✓ Shaking motion estimation complete -------------")
    
    return profile

def run_sorting(session, probe_id=0, protocol=None):
    profile_cls = get_recording_profile(session, probe_id)
    profile = profile_cls(session, probe_id, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()
    profile.prep_session_data()
    
    # Run sorting
    profile.spike_sorting()
    print(f"-------------- ✓ Spike sorting complete --------------")

    return profile

def run_postprocess(session, probe_id=0, protocol=None):
    profile_cls = get_recording_profile(session, probe_id)
    profile = profile_cls(session, probe_id, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()
    profile.prep_session_data()
    profile.spike_sorting()
    
    profile.postprocessing()
    print(f"-------------- ✓ Postprocessing complete --------------")
    
    profile.quality_metrics()
    print(f"-------------- ✓ Quality metrics calculated --------------")

    save_protocol_to_dict(protocol, profile.full_hash, session)

    return profile

def run_widgets(session, probe_id=0, protocol=None, job_id=0, n_chunks=1):
    
    profile_cls = get_recording_profile(session, probe_id)
    profile = profile_cls(session, probe_id, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()
    profile.prep_session_data()
    
    profile.post_widgets(job_id=job_id, n_chunks=n_chunks)
    print(f"-------------- ✓ Plotting widgets --------------")

    return profile


def profile_to_mat(session, protocol=None):
    profile_cls = get_recording_profile(session, 0)
    profile = profile_cls(session, 0, Path(PROTOCOLS_PATH) / protocol)

    profile.load_metadata()
    profile.load_protocol()
    profile.prep_session_data()

    # Get only the part after "sorting"
    if profile.sorter_path is not None:
        sorter_path = str(profile.sorter_path).split("sorting", 1)[1].lstrip("/\\")
    else:
        sorter_path = None

    return sorter_path
