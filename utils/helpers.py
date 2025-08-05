from pathlib import Path
import json
import os
import shutil
import spikeinterface.full as si

def make_folder_paths(raw_data_path, session_name, probe_id, drift_correct, sorter_type=None, run_suffix=None):

    # Set up directories, check if data has already been preprocessed and sorted 
    data_folder  =  Path(raw_data_path) / session_name

    # Load in metadata that contains info about recording session
    with open(Path(data_folder) / "metadata.json", 'r') as f:
        metadata = json.load(f)

    # Define output folder path based on probe type used
    if metadata["probe_types"][probe_id] == "neuropixel":
        output_folder  =  Path(data_folder) / f"{session_name}_imec{probe_id}" 
    else:
        output_folder  =  Path(data_folder) / f"{session_name}" 

    # Define where to output 'preprocessed' data, named based on which type of drift correction applied
    preprocess_folder = Path(output_folder) / "preprocess" / f"{drift_correct}"
    
    if sorter_type is not None:
        if run_suffix is not None:
            sorter_folder   =  Path(output_folder) / sorter_type / f"si_{drift_correct}_{run_suffix}"
        else:
            sorter_folder   =  Path(output_folder) / sorter_type / f"si_{drift_correct}"
    else:
        sorter_folder = None

    return data_folder, preprocess_folder, sorter_folder, metadata

def format_param_suffix(default_params, custom_params):
    """
    Compare custom_params to default_params and return a formatted suffix string
    for folder naming, including only parameters that differ.

    Parameters:
        default_params (dict): Dictionary of default parameters.
        custom_params (dict): Dictionary of custom parameters.

    Returns:
        str: Formatted suffix string (e.g., 'wr12_mcd24_ccg20').
    """

    # Define nickname mapping
    nicknames = {
        'batch_size': 'bs',
        'th_universal': 'thu',
        'th_learned': 'thl',
        'artifact_threshold': 'arth',
        'whitening_range': 'wr',
        'highpass_cutoff': 'hpc',
        'binning_depth': 'bdp',
        'sig_interp': 'sip',
        'min_template_size': 'mts',
        'template_sizes': 'tsz',
        'nearest_channels': 'nch',
        'nearest_templates': 'ntp',
        'max_channel_distance': 'mcd',
        'max_peels': 'mpl',
        'templates_from_data': 'tfd',
        'n_templates': 'ntm',
        'n_pcs': 'npc',
        'th_single_ch': 'ths',
        'acg_theshold': 'acg',
        'ccg_threshold': 'ccg',
        'cluster_neighbors': 'cln',
        'cluster_downsampling': 'cld',
        'max_cluster_subset': 'mcs',
        'x_centers': 'xct',
        'duplicate_spike_ms': 'dsm',
        'position_limit': 'plm',
        'do_car': 'car',
        'invert_sign': 'inv',
        'do_correction': 'cor',
        'skip_kilosort_preprocessing': 'skp',
        'keep_good_only': 'kgo'
    }

    # Params to ignore
    exclude_keys = {'clear_cache','save_extra_vars','save_preprocessed_copy','torch_device', 'delete_recording_dat', 'bad_channels',
                    'pool_engine','n_jobs','chunk_duration','progress_bar','mp_context','max_threads_per_worker', 'drift_smoothing', 'use_binary_file'}

    suffix_parts = []

    # Ensure custom_params are in the same order as default_params
    for key in default_params:
        if key in custom_params and key not in exclude_keys:
            default_val = default_params[key]
            custom_val = custom_params[key]

            if custom_val != default_val:
                nickname = nicknames.get(key, key)

                # Format the value
                if isinstance(custom_val, float) and 0 < custom_val < 1:
                    val_str = f"{int(custom_val * 100)}"
                elif isinstance(custom_val, bool):
                    val_str = "1" if custom_val else "0"
                else:
                    val_str = str(custom_val)

                suffix_parts.append(f"{nickname}{val_str}")

    return "_".join(suffix_parts)

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
