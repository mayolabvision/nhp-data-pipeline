function rawRipple_to_binaryFile(data_path,probes_path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

json_text = fileread(fullfile(data_path,'metadata.json'));
metadata = jsondecode(json_text);

pr_configs = metadata.probe_config(cellfun(@(q) isequal(q,'plexon'), metadata.probe_type));
hw_configs = metadata.hardware_config(cellfun(@(q) isequal(q,'plexon'), metadata.probe_type));

num_probes = numel(pr_configs);
if num_probes == 0
    return
end

full_bin_path = fullfile(data_path,[metadata.sess_name,'.bin']); 
if exist(full_bin_path, 'file') == 2
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

for nevnum = 1:length(nevnames)
    nevpath = nevpaths{nevnum};
    this_task = tasks{nevnum};

    if ~contains(this_task,'fstm')
        fprintf('\n---- loading raw signal for %s ----\n', this_task);
        [~, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false, 'READ_LFP', false);

        if ~isempty(out_ns5.data(ismember(out_ns5.hdr.label, string(1:512)),1))
            this_ns5 = [];
            for p = 1:num_probes
                this_prb = load(fullfile(probes_path,[pr_configs{p},'.mat']));
                this_hw = hw_configs{p};

                if isequal(this_hw,"elecA")
                    these_chans = this_prb.chanMap;
                elseif isequal(this_hw,"elecB")
                    these_chans = this_prb.chanMap+128;
                elseif isequal(this_hw,"elecC")
                    these_chans = this_prb.chanMap+256;
                elseif isequal(this_hw,"elecD")
                    these_chans = this_prb.chanMap+384;
                end

                % Channels x samples
                this_ns5 = [this_ns5; out_ns5.data(ismember(out_ns5.hdr.label, string(these_chans)),:)];
            end

            fid_write = fopen(full_bin_path, 'a'); % Open file in append mode ('a')
            fwrite(fid_write, this_ns5', 'double');

            fclose(fid_write);  
        else
            fprintf('\n---- no raw signal for %s ----\n', this_task);
            delete(full_bin_path);
            return
        end

    end
end

end
