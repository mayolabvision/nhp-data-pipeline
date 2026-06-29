from pathlib import Path
import json
import os
import shutil
import hashlib
from spikeinterface.sorters import get_default_sorter_params

def shorten_filename(filepath, key_file=None):
    """
    Shorten a filepath whose name exceeds OS limits by hashing the long suffix
    after a fixed prefix (everything up to and including the first '-').

    The prefix is preserved verbatim; the remainder is replaced with the first
    16 hex chars of its SHA-256 digest. A JSON key file records the mapping so
    the original suffix can always be recovered.

    Args:
        filepath (str | Path): Full path to the file (need not exist yet).
        key_file (str | Path | None): Path to the JSON key file. Defaults to
            a file named 'filename_key.json' in the same directory.

    Returns:
        Path: The shortened path.
    """
    filepath = Path(filepath)
    parent = filepath.parent
    name = filepath.name

    # Split on the first '-' to isolate the fixed prefix
    split = name.split('-', 1)
    if len(split) == 1:
        return filepath  # nothing to shorten

    prefix, suffix_with_ext = split[0], split[1]

    # Separate the file extension from the suffix
    dot_idx = suffix_with_ext.rfind('.')
    if dot_idx != -1:
        long_suffix = suffix_with_ext[:dot_idx]
        ext = suffix_with_ext[dot_idx:]
    else:
        long_suffix = suffix_with_ext
        ext = ''

    short_hash = hashlib.sha256(long_suffix.encode()).hexdigest()[:16]
    short_name = f"{prefix}-{short_hash}{ext}"
    short_path = parent / short_name

    if key_file is None:
        key_file = parent / 'filename_key.json'
    key_file = Path(key_file)

    # Load existing key file or start fresh
    if key_file.exists():
        with open(key_file, 'r') as f:
            key_map = json.load(f)
    else:
        key_map = {}

    if short_hash not in key_map:
        key_map[short_hash] = long_suffix
        with open(key_file, 'w') as f:
            json.dump(key_map, f, indent=2)

    return short_path


def restore_filename(short_path, key_file=None):
    """
    Reverse shorten_filename: look up the short hash in the key file and
    return the original full path.

    Args:
        short_path (str | Path): The shortened path returned by shorten_filename.
        key_file (str | Path | None): Path to the JSON key file. Defaults to
            'filename_key.json' in the same directory as short_path.

    Returns:
        Path: The original path, or short_path unchanged if not found.
    """
    short_path = Path(short_path)
    parent = short_path.parent
    name = short_path.name

    split = name.split('-', 1)
    if len(split) == 1:
        return short_path

    prefix, hash_with_ext = split[0], split[1]
    dot_idx = hash_with_ext.rfind('.')
    if dot_idx != -1:
        short_hash = hash_with_ext[:dot_idx]
        ext = hash_with_ext[dot_idx:]
    else:
        short_hash = hash_with_ext
        ext = ''

    if key_file is None:
        key_file = parent / 'filename_key.json'
    key_file = Path(key_file)

    if not key_file.exists():
        return short_path

    with open(key_file, 'r') as f:
        key_map = json.load(f)

    if short_hash not in key_map:
        return short_path

    original_name = f"{prefix}-{key_map[short_hash]}{ext}"
    return parent / original_name

def make_output_structure(session, protocol, probe_id=0, raw_data_path='/ix1/pmayo/lab_NHPdata'):
    
    # Where raw data is stored
    data_folder = Path(raw_data_path) / session
    
    # Load in metadata that contains info about recording session
    with open(os.path.join(raw_data_path, session, 'metadata.json'), 'r') as f:
        metadata = json.load(f)

    if metadata["probe_types"][probe_id] == "neuropixel":
        # Path of session data, for loading SpikeGLX into SI
        output_folder  =  Path(data_folder) / f"{session}_imec{probe_id}"

        # Path to put motion data
        os.makedirs(output_folder / 'motion', exist_ok=True)
        
        motion = protocol['motion_correction']
        motion_suffix = motion['drift_preset'] + '_' + '_'.join(
            f"{(''.join(word[0] for word in key.split('_')) if '_' in key else key[:3])}{value}"
            for key, value in motion.items()
            if key != 'drift_preset'
        )

        motion_folder = output_folder / 'motion' / motion_suffix
        motion_folder.mkdir(parents=True, exist_ok=True)
        
        preprocess_folder = output_folder / 'preprocess' / motion_suffix
        preprocess_folder.mkdir(parents=True, exist_ok=True)

        # Path to put sorting outputs
        default_params = get_default_sorter_params(sorter_name_or_class=protocol['sorting']['sorter_name'])

        os.makedirs(output_folder / protocol['sorting']['sorter_name'], exist_ok=True)

        sorter_suffix = format_sorter_suffix(default_params, protocol['sorting'])
        sorter_folder = output_folder / 'sorting' / protocol['sorting']['sorter_name'] / f'{motion_suffix}_{sorter_suffix}'
       
        # Path to save MATLAB table
        table_path = data_folder / 'tables' / f"{session}-{motion_suffix}-{protocol['sorting']['sorter_name']}_{sorter_suffix}.mat"
        os.makedirs(data_folder / 'tables', exist_ok=True)

        # Folder to save figures
        fig_folder = data_folder / 'figs' / f"{motion_suffix}-{protocol['sorting']['sorter_name']}_{sorter_suffix}"
        os.makedirs(fig_folder, exist_ok=True)

    else:
        output_folder = data_folder

    return metadata, data_folder, motion_folder, preprocess_folder, sorter_folder, table_path, fig_folder

def format_sorter_suffix(default_params, custom_params):
    """
    Generate a suffix based on differing parameters, using automatic nicknames.

    Parameters:
        default_params (dict): Default parameter dictionary.
        custom_params (dict): Custom parameter dictionary.

    Returns:
        str: Suffix string for folder naming (e.g., 'wr12_mcd24').
    """

    # Optional manual overrides
    nickname_overrides = {
        'batch_size': 'bs',
        'highpass_cutoff': 'hpc',
        'do_correction': 'cor',
        # Add others if you really want custom ones
    }

    # Exclude keys that shouldn't appear in suffix
    exclude_keys = {
        'sorter_name', 'clear_cache', 'save_extra_vars', 'save_preprocessed_copy',
        'torch_device', 'delete_recording_dat', 'bad_channels', 'pool_engine',
        'n_jobs', 'chunk_duration', 'progress_bar', 'mp_context',
        'max_threads_per_worker', 'drift_smoothing', 'use_binary_file'
    }

    def auto_nickname(key):
        if key in nickname_overrides:
            return nickname_overrides[key]
        elif '_' in key:
            return ''.join(part[0] for part in key.split('_'))
        else:
            return key[:3]

    suffix_parts = []

    for key in default_params:
        if key in custom_params and key not in exclude_keys:
            default_val = default_params[key]
            custom_val = custom_params[key]

            if custom_val != default_val:
                nickname = auto_nickname(key)

                # Format the value
                if isinstance(custom_val, float) and 0 < custom_val < 1:
                    val_str = f"{int(custom_val * 100)}"
                elif isinstance(custom_val, bool):
                    val_str = "1" if custom_val else "0"
                else:
                    val_str = str(custom_val)

                suffix_parts.append(f"{nickname}{val_str}")

    return "_".join(suffix_parts)

def print_recording_details(recording):
    print("Sampling frequency:", recording.get_sampling_frequency())
    print("Number of channels:", recording.get_num_channels())
    print("Number of segments:", recording.get_num_segments())
    print("Number of samples:", recording.get_num_samples(segment_index=0))
    print("Duration of recording (min):", round((recording.get_num_samples(segment_index=0)/recording.get_sampling_frequency())/60))
    print("Data dtype:", recording.get_dtype())

def write_recording_details(rec, save_path):
    num_channels = rec.get_num_channels()
    sampling_frequency = rec.get_sampling_frequency()
    num_segments = rec.get_num_segments()
    num_samples = rec.get_num_samples(segment_index=0)
    total_time = num_samples / sampling_frequency
    dtype = rec.get_dtype()

    # opening a .txt file to output the values of raw_rec
    with open(save_path, 'w') as file:
        # Step 4: Write the metadata to the file
        file.write(f'Number of channels: {num_channels}\n')
        file.write(f'Sampling frequency: {sampling_frequency} Hz\n')
        file.write(f'Number of segments: {num_segments}\n')
        file.write(f'Number of samples: {num_samples}\n')
        file.write(f'Total time: {total_time} seconds\n')
        file.write(f'Data type: {dtype}\n')

def shorten_sorter_path(sorter_path: str, key_dir: str | Path | None = None) -> str:
    if key_dir is None:
        key_dir = Path.cwd()
    key_dir = Path(key_dir)

    short_id = hashlib.sha256(sorter_path.encode()).hexdigest()[:16]

    key_file = key_dir / "sorter_path_key.json"
    key_map = json.loads(key_file.read_text()) if key_file.exists() else {}

    field = f"id_{short_id}"
    if field not in key_map:
        key_map[field] = sorter_path
        key_file.write_text(json.dumps(key_map, indent=2))

    return short_id

