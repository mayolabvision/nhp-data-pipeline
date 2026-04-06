import json
import hashlib
from pathlib import Path
import copy

from spikeinterface.preprocessing import get_motion_parameters_preset
from spikeinterface.sorters import get_default_sorter_params

def print_protocol(profile):
    print(f"===================================================================")
    print("Preprocess hash: ", profile.preprocess_hash)
    print("Motion PP hash:  ", profile.pp_hash)
    print("Motion hash:     ", profile.motion_hash)
    print("Shake hash:      ", profile.shake_hash)
    print("Trim hash:       ", profile.trim_hash)
    print("Sorter hash:     ", profile.sorter_hash)
    print("FULL HASH:       ", profile.full_hash)
    print(f"===================================================================") 

def save_params(params_path: Path, params: dict):
    """Save params dict to JSON file at params_path if the file does not already exist."""
    if not params_path.is_file():
        params_path.parent.mkdir(parents=True, exist_ok=True)  # ensure directory exists
        with open(params_path, "w") as f:
            json.dump(params, f, indent=4)

def get_preprocess_hash(preprocessing_params: dict):
    # 1. Serialize the params with sorted keys for consistency
    params_str = json.dumps(preprocessing_params, sort_keys=True)
    
    # 2. Hash the serialized string
    params_hash = hashlib.md5(params_str.encode('utf-8')).hexdigest()
    
    #print("preprocess_hash:", params_hash)
    
    return params_hash

def get_sparse_hash(sparse_params: dict):
    if not sparse_params:
        return "notsparse"

    # 1. Serialize the params with sorted keys for consistency
    params_str = json.dumps(sparse_params, sort_keys=True)
    
    # 2. Hash the serialized string
    params_hash = hashlib.md5(params_str.encode('utf-8')).hexdigest()
    
    return params_hash


def get_motion_hash(motion_params: dict):
    # Handle empty motion_params
    if not motion_params:
        return "nodrift", "nodrift", {}, {}

    drift_preset = motion_params['estimate_motion_kwargs']['method']
    drift_preset = drift_preset.split("_")[0] if "dredge" in drift_preset else drift_preset

    default_params = get_motion_parameters_preset(drift_preset)

    # 1. Get custom params (only differences + 'method')
    custom_params = prune_params(motion_params, default_params)
    #print("Custom params:", custom_params)
    
    pp_str = json.dumps(custom_params['preprocessing'], sort_keys=True)
    pp_hash = hashlib.md5(pp_str.encode('utf-8')).hexdigest()
    #print("pp_hash:", pp_hash)
   
    custom_params_motion = {k: v for k, v in custom_params.items() if k != 'preprocessing'} 
    params_str = json.dumps(custom_params_motion, sort_keys=True)
    params_hash = hashlib.md5(params_str.encode('utf-8')).hexdigest()
    #print("motion_hash:", params_hash)
    
    # 2. Get full params by merging custom into default
    full_params = merge_params(default_params, custom_params)

    return pp_hash, params_hash, motion_params['preprocessing'], full_params

def get_trim_hash(motion_params: dict):
    drift_preset = motion_params['estimate_motion_kwargs']['method']
    drift_preset = drift_preset.split("_")[0] if "dredge" in drift_preset else drift_preset

    default_params = get_motion_parameters_preset(drift_preset)

    # 1. Get custom params (only differences + 'method')
    custom_params = prune_params(motion_params, default_params)
    #print("Custom params:", custom_params)
    
    pp_str = json.dumps(custom_params['preprocessing'], sort_keys=True)
    pp_hash = hashlib.md5(pp_str.encode('utf-8')).hexdigest()
    
    trim_str = json.dumps(custom_params['window_params'], sort_keys=True)
    trim_hash = hashlib.md5(trim_str.encode('utf-8')).hexdigest()

    custom_params_motion = {k: v for k, v in custom_params.items() 
                        if k not in ['window_params', 'preprocessing']}
    params_str = json.dumps(custom_params_motion, sort_keys=True)
    params_hash = hashlib.md5(params_str.encode('utf-8')).hexdigest()
    #print("motion_hash:", params_hash)
    
    # 2. Get full params by merging custom into default
    full_params = merge_params(default_params, custom_params)

    return pp_hash, params_hash, trim_hash, full_params

def get_sorter_hash(sorter_params: dict):
    default_params = get_default_sorter_params(sorter_params['sorter_name'])

    # 1. Get custom params (only differences + 'method')
    custom_params = prune_params(sorter_params, default_params)
    #print("Custom params:", custom_params)
    
    params_str = json.dumps(custom_params, sort_keys=True)
    params_hash = hashlib.md5(params_str.encode('utf-8')).hexdigest()
    #print("sorter_hash:", params_hash)
    
    # 2. Get full params by merging custom into default
    full_params = merge_params(default_params, custom_params)
    
    return params_hash, full_params, custom_params

def prune_params(custom, default):
    """
    Recursively keep keys in 'custom' only if:
    - key == 'method' (always keep)
    - OR value differs from default's value
    - 'preprocessing' is always included as-is, without comparison
    """
    pruned = {}
    for key, val in custom.items():
        # Always include 'preprocessing' without comparing
        if key == 'preprocessing':
            pruned[key] = val
            continue

        # If key not in default, keep it
        if key not in default:
            pruned[key] = val
            continue

        # If values are dicts, recurse
        if isinstance(val, dict) and isinstance(default[key], dict):
            nested = prune_params(val, default[key])
            if nested:  # only keep if nested dict not empty
                pruned[key] = nested
        else:
            # Keep if different
            if val != default[key]:
                pruned[key] = val
    return pruned

def merge_params(default, custom):
    """Recursively merge custom into default, overwriting default values."""
    merged = copy.deepcopy(default)
    for key, val in custom.items():
        if isinstance(val, dict) and key in merged and isinstance(merged[key], dict):
            merged[key] = merge_params(merged[key], val)
        else:
            merged[key] = val
    return merged

