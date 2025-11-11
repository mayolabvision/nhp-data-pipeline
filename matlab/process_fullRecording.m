function process_fullRecording(session_name,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% session_name = 'kendra_scrappy_0136a_g0' 
% session_name is the name of a datafolder, which contains ripple data
% (.ns5, .nev, etc...) and folders of spikeglx imec probe data, catgt align pulses 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default paths to add to the MATLAB path    
defaultRAW_PATH  =  '/Volumes/lab_NHPdata';
defaultOUT_PATH  =  '/Volumes/home/DATA';
defaultNEV_PATH  =  '/Users/kendranoneman/Packages/nevutils';

p = inputParser;
addRequired(p, 'session_name', @ischar);
addParameter(p, 'RAW_DATA_PATH', defaultRAW_PATH, @ischar); 
addParameter(p, 'OUT_DATA_PATH', defaultOUT_PATH, @ischar); 
addParameter(p, 'NEVUTIL_PATH', defaultNEV_PATH, @ischar);
addParameter(p, 'NASNET_PATH', [], @ischar); % only used for plex
addParameter(p, 'SORTER_PATH', [], @ischar);

% Parse inputs
parse(p, session_name, varargin{:});
RAW_PATH       =  p.Results.RAW_DATA_PATH;
OUT_PATH       =  p.Results.OUT_DATA_PATH;
NEV_PATH       =  p.Results.NEVUTIL_PATH;
NET_PATH       =  p.Results.NASNET_PATH;
SORTER_PATH    =  p.Results.SORTER_PATH;

addpath(genpath(NEV_PATH));

if isequal(SORTER_PATH, 'None')
    SORTER_PATH = [];
end

% Get the directory of the current script or function
currentDir = fileparts(mfilename('fullpath'));
parentDirOneLevelUp = fileparts(currentDir);
addpath(genpath(parentDirOneLevelUp));
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(fullfile(OUT_PATH, session_name), 'dir'), mkdir(fullfile(OUT_PATH, session_name)); end

session_path = fullfile(RAW_PATH, session_name);
filePattern = fullfile(RAW_PATH, session_name, '*.ns5');

S1 = struct();

if isfile(fullfile(session_path,'metadata.json'))
    metadata = loadMetadataJSON(fullfile(session_path,'metadata.json'));

    eye_chan_labels = metadata.HEeye_VEeye_diode_pupil(1:2)';
    diode_chan_label = metadata.HEeye_VEeye_diode_pupil{3};
    pupil_chan_label = metadata.HEeye_VEeye_diode_pupil{4};

    S1.sess_name = metadata.sess_name;
    S1.metadata = metadata;
else
    metadata.probe_type = 'behavior';

    eye_chan_labels = {'10241','10242'};
    diode_chan_label = '10243';
    pupil_chan_label = '10244';

    S1.sess_name = session_name;
end

% Create the search pattern to find files that start with 'filename' and end with '.ns5'
raw_files = dir(filePattern);
raw_filenames = {raw_files.name}.';
nevnames = cellfun(@(q) q(1:end-4), raw_filenames, 'uni', 0);
raw_filepaths = arrayfun(@(x) fullfile(x.folder, x.name), raw_files, 'UniformOutput', false);

recording_times = cellfun(@(l) l.hdr.timeOrigin, cellfun(@(q) read_nsx(q,'readdata',false), raw_filepaths, 'uni', 0), 'uni', 0);
[~,idx] = sort(recording_times);
nevnames = nevnames(idx);
nevpaths = raw_filepaths(idx);

% Define possible task keywords
task_keywords = {'rfmp', 'rfMapping', 'purs', 'pursuit', 'mdir', 'dirmem', 'fstm', 'cfix', 'frvw', 'oflu'};
tasks = cell(size(nevnames));

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

taskTypes = unique(cellfun(@(q) regexp(q, '[a-zA-Z]+', 'match', 'once'), tasks, 'uni', 0));
disp(tasks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extracting Ripple data from nev/ns5 datafiles 
tic

goodFlag = true;
prev_tempdata = struct([]);
rip_time_start = 0;
for nevnum = 1:length(nevnames) % loop through nev files, in chronological
    nevpath = nevpaths{nevnum};
    this_task = tasks{nevnum};

    fprintf('\n---- generating nev_out for %s ----\n', this_task);

    %------------------------------------- NEUROPIXELS -------------------------------------%
    if ismember('neuropixel',metadata.probe_type)
        % Load in align/sync pulses 
        if nevnum==1
            alignCodes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name],[session_name,'_tcat.nidq.bfv_8_0_9.txt']));
            alignTimes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name],[session_name,'_tcat.nidq.bft_8_0_9.txt']));
            alignTimes = alignTimes(alignCodes>0);
        end

        [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false, 'READ_LFP', false, 'alignPulseEnabled', true);
        startAcquisition = datetime(out_ns5.hdr.timeOrigin, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS');

        [dat, ~] = format_datTrials(nev, out_ns5, 'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);

        firstSyncPulse = startAcquisition + seconds(dat(1).trialcodes(2,3));
        ripple_pulse_timeStamps = cellfun(@(w) (firstSyncPulse + seconds(w)) - seconds(dat(1).trialcodes(2,3)), cellfun(@(q) q(2,3), {dat.trialcodes}.', 'uni', 0), 'uni', 1);     

        if nevnum == 1
            np_pulse_timeStamps = cellfun(@(w) firstSyncPulse + seconds(w), num2cell(alignTimes - alignTimes(1)), 'uni', 1);
        end

        [np_mask, ripple_mask] = match_syncPulses_RipToNP(np_pulse_timeStamps, ripple_pulse_timeStamps);
        fprintf('\n dat has %d rows\n', numel(dat))
        fprintf('np_mask = %d/%d, ripple_mask = %d/%d \n', sum(np_mask), length(np_mask), sum(ripple_mask), length(ripple_mask))  

        if contains(session_name, 'kendra_scrappy_0136a') 
            [dat,these_alignTimes,goodFlag] = fix_specificSessions(session_name,np_mask,ripple_mask,alignTimes,dat,goodFlag);
        else
            these_alignTimes = alignTimes(np_mask);
            if sum(np_mask) < length(ripple_mask)
                dat = dat(ripple_mask);
                fprintf('\n dat NOW has %d rows', numel(dat))
            end
        end 

        tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);

    %------------------------------------- PLEXON -------------------------------------%
    elseif ismember('plexon',metadata.probe_type)
        if contains(this_task, 'fstm')
            continue
        end
            
        if ~isempty(NET_PATH)
            addpath(genpath(NET_PATH));

            if exist([nevpath,'.ns2'], 'file') == 2
                [nev, out_ns5, out_ns2] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'), 'READ_LFP', true);
                lfp = extract_lfpData(nev,out_ns2,mappings.ripChan_num); 

                [dat, ~] = format_datTrials(nev, out_ns5, 'NEURAL_CHANNELS', mappings.ripChan_num, 'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);
                tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task, 'LFP', lfp);
            else
                [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'), 'READ_LFP', false);
                [dat, ~] = format_datTrials(nev, out_ns5, 'NEURAL_CHANNELS', mappings.ripChan_num, 'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);
                tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);
            end
        else
            [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false);

            [dat, ~, tempdata] = format_datTrials(nev, out_ns5, 'PREV_TEMPDATA', prev_tempdata, ...
                                                 'EYE_CHAN_LABELS', eye_chan_labels, ...
                                                 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);

            prev_tempdata = tempdata;
            tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);
            tbl.time_sec = tbl.time_sec+rip_time_start;

            rip_time_start = rip_time_start + double(out_ns5.hdr.nSamples/out_ns5.hdr.Fs);

        end

    %------------------------------------- FHC SINGLE ELECTRODE -------------------------------------%
    elseif ismember('fhc',metadata.probe_type)
        addpath(genpath(NET_PATH));

        [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'));

        hw_config = metadata.hardware_config{1};
        parts = strsplit(hw_config, '_');
        this_hw = parts{1};
        this_chan = str2double(parts{2});

        if isequal(this_hw,"elecA")
            neural_chan = this_chan;
        elseif isequal(this_hw,"elecB")
            neural_chan = this_chan+128;
        elseif isequal(this_hw,"elecC")
            neural_chan = this_chan+256;
        elseif isequal(this_hw,"elecD")
            neural_chan = this_chan+384;
        end

        [dat, ~, tempdata] = format_datTrials(nev, out_ns5, 'PREV_TEMPDATA', prev_tempdata, ...
                                              'NEURAL_CHANNELS', neural_chan, 'EYE_CHAN_LABELS', eye_chan_labels, ...
                                              'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);
        
        if ~isempty(dat)
            prev_tempdata = tempdata;
            tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);
        else
            tbl = [];
        end

        [~, fname, ~] = fileparts(nevpath);   % 'kendra_scrappy_0066a_mdir1'
        this_task = erase(fname, session_name); % 'a_mdir1'

    %------------------------------------- BEHAVIOR ONLY -------------------------------------%
    else 
        [nev, out_ns5, ~] = extract_nevout(nevpath);
        if ~isempty(nev)
            [dat, ~] = format_datTrials(nev, out_ns5, 'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);
            tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);
        else
            dat = []; tbl = [];
        end
    end

    % Convert structures to a cell array of string representations
    if ~isempty(tbl)
        merged_struct = merge_taskParams(tbl);
        S1.(this_task).params = merged_struct;
    end
    S1.(this_task).hdr = out_ns5.hdr;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Extracting sorting data from SpikeInterface outputs (e.g., kilosort)
    if SORTER_PATH
        hashes = split(SORTER_PATH, '-');
        preprocess_hash = hashes{1}; pp_hash = hashes{2}; motion_hash = hashes{3}; %sorter_hash = hashes{4};

        sorting_all = []; motion_info = struct([]);
        for probe = 1:numel(metadata.probe_type)

            si_path = fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], 'sorting', SORTER_PATH);

            if isequal(metadata.probe_type{probe},'plexon')
                ripple_info = loadMetadataJSON(fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], 'ripple_info.json'));

                % Pull out SpikeInterface sorting outputs
                [spikes_perTrial, sorting, ~] = parse_SortingToTbl(tbl, fullfile(si_path, 'sorter_output'), 'Fs', ripple_info.Fs);

            elseif isequal(metadata.probe_type{probe},'neuropixel')
                lfp_meta = readMetaFile(fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], [session_name, '_t0.', metadata.hardware_config{probe}, '.lf.meta']));
                ap_meta = readMetaFile(fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], [session_name, '_t0.', metadata.hardware_config{probe}, '.ap.meta']));

                lfp_fs = lfp_meta.imSampRate;
                trial_starts_sec = cellfun(@(q,v) q-(v./1000), num2cell(these_alignTimes), tbl.ALIGN_PULSE(:,1), 'uni', 1);
                trial_ends_sec = trial_starts_sec + (tbl.END_TRIAL./1000);
    
                trial_starts_lfp_samp = floor(trial_starts_sec * lfp_fs);
                trial_ends_lfp_samp = ceil(trial_ends_sec * lfp_fs);
    
                tbl.imec_sec = [trial_starts_sec trial_ends_sec];
                tbl.imecLFP_samp = [trial_starts_lfp_samp trial_ends_lfp_samp];

                % Pull out SpikeInterface sorting outputs
                [spikes_perTrial,sorting,~] = parse_SortingToTbl(tbl,fullfile(si_path,'sorter_output'), 'NP_ALIGN_PULSES', these_alignTimes, 'Fs', ap_meta.imSampRate);

            end

            tbl.(sprintf('spiketimes_%d',probe)) = spikes_perTrial;  

            if nevnum==1
                sorting.probe_index = probe;
                fields = fieldnames(sorting);
                fields(strcmp(fields, 'probe_index')) = [];
                sorting = orderfields(sorting, ['probe_index'; fields]);

                if isfile(fullfile(si_path,'quality_metrics','cluster_metrics.csv'))
                    metrics = parse_clusterMetrics(si_path);
                    sorting.clusters = [metrics removevars(sorting.clusters, 'cluster_id')];
                else
                    sorting.clusters.sess_name = repmat(metadata.sess_name,height(sorting.clusters),1);
                    sorting.clusters.probe_id = repmat(probe-1,height(sorting.clusters),1);
                    sorting.clusters.probe_index = repmat(probe,height(sorting.clusters),1);
                    sorting.clusters.probe_label = repmat(metadata.probe_label{probe},height(sorting.clusters),1);
                    sorting.clusters.probe_type = repmat(metadata.probe_type{probe},height(sorting.clusters),1);
                    sorting.clusters.probe_config = repmat(metadata.probe_config{probe},height(sorting.clusters),1);
                    sorting.clusters.hardware_config = repmat(metadata.hardware_config{probe},height(sorting.clusters),1);
                    sorting.clusters.probe_depth_mm = repmat(metadata.probe_depth_mm(probe),height(sorting.clusters),1);
                    sorting.clusters.probe_gridHole = repmat({S1.metadata.probe_gridHole{probe}},height(sorting.clusters),1);
                end
                sorting.clusters.probe_index = repmat(probe,height(sorting.clusters),1);
                sorting.clusters = movevars(sorting.clusters,{'probe_index'},'After','sess_name');

                motion_path = fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], 'preprocess', preprocess_hash, pp_hash, motion_hash);
                load(fullfile(motion_path,'motion.mat'));
                load(fullfile(motion_path,'depth_bins.mat'));
                load(fullfile(motion_path,'time_bins.mat'));

                motion_info(probe).probe_index = probe;
                motion_info(probe).probe_label = metadata.probe_label{probe};
                motion_info(probe).probe_type = metadata.probe_type{probe};
                motion_info(probe).probe_config = metadata.probe_config{probe};
                motion_info(probe).hardware_config = metadata.hardware_config{probe};
                motion_info(probe).probe_depth_mm = metadata.probe_depth_mm(probe);
                motion_info(probe).motion = double(motion);
                motion_info(probe).depth_bins = double(depth_bins);
                motion_info(probe).time_bins = double(time_bins);

                sorting_all = [sorting_all; sorting];

                if probe==numel(metadata.probe_type)
                    jsonStr = fileread(fullfile(si_path,'params.json'));
                    protocolStruct = jsondecode(jsonStr);
                    S1.protocol = protocolStruct;
                    S1.sorting = sorting_all;
                end
            end
        end
    end

    S1.(this_task).dat = dat;
    S1.(this_task).tbl = tbl;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Final steps and saving table 

% Re-order fields of struct 
ff = fieldnames(S1);
newOrder = ff(ismember(ff, {'sess_name', 'metadata', 'protocol'}));
if ismember('sorting', ff)
    newOrder = [newOrder; 'sorting', 'motion_info'];
end
newOrder = [newOrder; ff(~ismember(ff, newOrder))];
S1 = orderfields(S1, newOrder);
S = unify_taskTables(S1,taskTypes);

if ~exist(fullfile(OUT_PATH, session_name, 'tables'), 'dir'), mkdir(fullfile(OUT_PATH, session_name, 'tables')); end 

if SORTER_PATH
    save(fullfile(OUT_PATH,session_name,'tables',sprintf('%s-%s.mat',session_name,SORTER_PATH)), 'S', '-v7.3');
else
    if NET_PATH
        save(fullfile(OUT_PATH,session_name,'tables',sprintf('%s-nasnet.mat',session_name)), 'S', '-v7.3');
    else
        save(fullfile(OUT_PATH,session_name,'tables',sprintf('%s.mat',session_name)), 'S', '-v7.3');
    end
end

tc = toc;
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf(sprintf('Total elapsed time was %2.2f minutes',tc/60))
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
end
