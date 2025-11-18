function rawRipple_to_binaryFile(data_path,probes_path,probe_index)
% Convert Ripple NS5 raw signals to a single binary file, for spike sorting
%
% This function reads raw neural signals from Ripple NS5 files in a session,
% selects channels corresponding to Plexon probes, concatenates them across 
% multiple tasks, and saves the combined data as a binary file (.bin) in 
% channels x samples format.
%
% Additionally, it creates a JSON file ('ripple_info.json') storing essential
% metadata such as:
%   - Sampling frequency (Fs)
%   - Number of channels and samples
%   - Data type information (MATLAB and Python)
%   - Gain and offset to convert to microvolts (uV)
%
% INPUTS:
%   data_path   : Path to the session folder containing NS5 files and metadata.json
%   probes_path : Path to folder containing Plexon probe configurations (*.mat)
%
% OUTPUT:
%   Creates a binary file [sess_name].bin in data_path containing the concatenated
%   raw signals (samples x channels) and a ripple_info.json file with metadata.
%
% NOTES:
%   - Only Plexon probes are processed.
%   - Tasks with 'fstm' in the name are skipped (microstim data)
%   - Raw NS5 data is assumed to already be in physical units (µV), so gain_to_uV = 1.
%   - If the output files already exist, the function returns without overwriting.

json_text = fileread(fullfile(data_path,'metadata.json'));
metadata = jsondecode(json_text);

if ~isequal(metadata.probe_type{probe_index},'plexon')
    return
end

hw_config = metadata.hardware_config{probe_index};

if exist(fullfile(data_path, [metadata.sess_name, '_', metadata.hardware_config{probe_index}], 'ripple_info.json'), 'file') == 2
    return
end

filePattern = fullfile(data_path, '*.ns5');
raw_files = dir(filePattern);
raw_filenames = {raw_files.name}.';
nevnames = cellfun(@(q) q(1:end-4), raw_filenames, 'uni', 0);
raw_filepaths = arrayfun(@(x) fullfile(x.folder, x.name), raw_files, 'UniformOutput', false);

recording_times = cellfun(@(l) l.hdr.timeOrigin, cellfun(@(q) read_nsx(q,'readdata',false), raw_filepaths, 'uni', 0), 'uni', 0);
[~,idx] = sort(recording_times);
nevnames = nevnames(idx);
nevpaths = raw_filepaths(idx);

% Define possible task keywords
task_keywords = {'rfmp', 'rfMapping', 'purs', 'pursuit', 'mdir', 'dirmem', 'fstm', 'cfix'};

% Initialize cell array for tasks
tasks = cell(size(nevnames));

% Loop through each nevname
for i = 1:numel(nevnames)
    name = nevnames{i};
    found = false;
    for j = 1:numel(task_keywords)
        pattern = [task_keywords{j}, '\w*'];  % keyword followed by letters/numbers
        match = regexp(name, pattern, 'match', 'once');
        if ~isempty(match)
            tasks{i} = match;
            found = true;
            break;
        end
    end
    if ~found
        tasks{i} = 'unknown';
    end
end

ns5_tasks = cell(sum(cellfun(@(q) ~contains(q,'fstm'), nevnames, 'uni', 1)),1);
for nevnum = 1:length(nevnames)
    nevpath = nevpaths{nevnum};
    this_task = tasks{nevnum};

    if ~contains(this_task,'fstm')
        fprintf('\n---- loading raw signal for %s ----\n', this_task);
        %[~, out_ns5, ~] = extract_nevout(nevpath, 'KEEP_INT', true);
        [~, out_ns5, ~] = extract_nevout(nevpath, 'KEEP_INT', false);

        if ~isempty(out_ns5.data(ismember(out_ns5.hdr.label, string(1:512)),1))
            this_prb = load(fullfile(probes_path,[metadata.probe_config{probe_index},'.mat']));

            if isequal(hw_config,"elecA")
                these_chans = this_prb.chanMap;
            elseif isequal(hw_config,"elecB")
                these_chans = this_prb.chanMap+128;
            elseif isequal(hw_config,"elecC")
                these_chans = this_prb.chanMap+256;
            elseif isequal(hw_config,"elecD")
                these_chans = this_prb.chanMap+384;
            end

            % Channels x samples
            this_ns5 = out_ns5.data(ismember(out_ns5.hdr.label, string(these_chans)),:);

            % Save each task as samples x channels 
            ns5_tasks{nevnum} = this_ns5';
        else
            fprintf('\n---- no raw signal for %s ----\n', this_task);
            return
        end
    end
end

% Samples x Channels (across all tasks)
ns5_data = vertcat(ns5_tasks{:});

% Saving data to binary file
full_bin_path = fullfile(data_path, [metadata.sess_name, '_', hw_config], 'raw_signal.bin'); 
if ~exist(fileparts(full_bin_path),'dir'), mkdir(fileparts(full_bin_path)); end
fid_data = fopen(full_bin_path, 'wb'); % Open file in append mode ('a')
%fwrite(fid_data, ns5_data, 'int16');
fwrite(fid_data, ns5_data, 'double');
fclose(fid_data);

% Storing important parameters to json
ripple_info = struct();
ripple_info.Fs = double(out_ns5.hdr.Fs);
ripple_info.num_samples = size(ns5_data,1);
ripple_info.num_channels = size(ns5_data,2);
%ripple_info.dtype_matlab = 'int16';
ripple_info.dtype_matlab = 'double';
%ripple_info.dtype_python = 'int16';
ripple_info.dtype_python = 'float64';
%ripple_info.gain_to_uV = out_ns5.hdr.scale(these_chans(1)); % e.g., 0.25 μV per int16 unit
%ripple_info.offset_to_uV = 0;

json_text = jsonencode(ripple_info, 'PrettyPrint', true);
fid_info = fopen(fullfile(data_path, [metadata.sess_name, '_', hw_config], 'ripple_info.json'), 'w');
fwrite(fid_info, json_text, 'char');
fclose(fid_info);

end
