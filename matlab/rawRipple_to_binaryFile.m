function rawRipple_to_binaryFile(data_path,probe_id,probes_path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

json_text = fileread(fullfile(data_path,'metadata.json'));
metadata = jsondecode(json_text);

if ~isequal(metadata.probe_type{probe_id},'plexon')
    return
end

prb_folder = fullfile(data_path, sprintf('%s_%s',metadata.hardware_config{probe_id},metadata.probe_label{probe_id}));
if exist(fullfile(prb_folder,'prb.mat'), 'file')
    return
end

if ~exist(prb_folder,'dir')
    mkdir(prb_folder);
end

full_bin_path = fullfile(prb_folder, 'raw.bin'); 
if exist(full_bin_path, 'file') == 2
    delete(full_bin_path);
end

prb = load(fullfile(probes_path,[metadata.probe_config{probe_id} '.mat']));

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
%disp(tasks)

for nevnum = 1:length(nevnames)
    nevpath = nevpaths{nevnum};
    this_task = tasks{nevnum};

    if ~contains(this_task,'fstm')
        [~, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false, 'READ_LFP', false);
        if ~isempty(out_ns5.data(ismember(out_ns5.hdr.label, string(1:512)),1))
            if isequal(metadata.hardware_config{probe_id},"elecA")
                these_chans = prb.chanMap;
            elseif isequal(metadata.hardware_config{probe_id},"elecB")
                these_chans = (prb.chanMap+129)-1;
            elseif isequal(metadata.hardware_config{probe_id},"elecC")
                these_chans = (prb.chanMap+257)-1;
            elseif isequal(metadata.hardware_config{probe_id},"elecD")
                these_chans = (prb.chanMap+385)-1;
            end

            % Channels x samples
            this_ns5 = out_ns5.data(ismember(out_ns5.hdr.label, string(these_chans)),:);

            fid_write = fopen(full_bin_path, 'a'); % Open file in append mode ('a')
            fwrite(fid_write, this_ns5, 'double');

            fclose(fid_write);  
        else
            fprintf('\n---- no raw signal for %s ----\n', this_task);
        end

    end
end

% Save new probe map to data_path, with correct channels
prb.chanMap = (1:size(this_ns5,1))';
prb.chanMap0ind = prb.chanMap - 1;

chanMap = prb.chanMap; chanMap0ind = prb.chanMap0ind; connected = prb.connected; 
kcoords = prb.kcoords; name = prb.name; xcoords = prb.xcoords; ycoords = prb.ycoords; 

new_prb_path = fullfile(prb_folder, 'prb.mat');
save(new_prb_path,"chanMap","chanMap0ind","connected","kcoords","name","xcoords","ycoords");

end
