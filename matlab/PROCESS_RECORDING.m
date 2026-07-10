function PROCESS_RECORDING(session_name,varargin)
%PROCESS_RECORDING  Process a single recording session into a unified data table.
%
%   PROCESS_RECORDING(session_name) reads the raw electrophysiology and/or
%   behavioral data for one session, aligns/extracts trials task-by-task,
%   optionally attaches spike-sorted units (nasnet or SpikeInterface), and
%   saves a single struct (S) of unified trial tables to disk.
%
%   Examples:
%       PROCESS_RECORDING('kendra_scrappy_0136a_g0')
%
%       PROCESS_RECORDING('kendra_scrappy_0136a_g0', 'RAW_DATA_PATH',
%       '/path/to/data', 'SORTER_HASH', '149823789-ajdslkj2u')
%
% ------------------------------------------------------------------------
% INPUT
%   session_name : name of a datafolder (e.g. 'kendra_scrappy_0136a_g0'),
%                  containing ripple data (.ns5, .nev, ...) and/or folders
%                  of spikeglx imec probe data + catgt align pulses.
%
% NAME-VALUE PARAMETERS (all optional; see defaults below)
%   RAW_DATA_PATH  : root folder containing raw session data
%   OUT_DATA_PATH  : root folder where processed output gets saved
%   NEVUTIL_PATH   : path to the nevutils toolbox (ripple nev/ns5 reading)
%   NASNET_PATH    : path to nasnet spike-sorting networks (plexon/FHC only)
%   HELPERS_PATH   : path to shared helper functions (behavior parsing, etc.)
%   SORTER_HASH    : hyphen-delimited hash string identifying a
%                    preprocess/motion/sorter run whose SpikeInterface
%                    outputs should be pulled in (empty/'None' to skip)
%   BEHAVIOR_ONLY  : true to skip all neural processing and only extract
%                    behavior (default false)
%   INCLUDE_ITI    : true to retain inter-trial-interval (ITI) rows in the
%                    output tables (default true)
%
% OUTPUT
%   Nothing is returned; a struct S (one field per task type, containing
%   unified trial tables/metadata) is saved as a .mat file under
%   OUT_DATA_PATH/session_name/tables/.
%
% WORKFLOW OVERVIEW
%   1. Parse inputs and set up paths.
%   2. Load session metadata (probe type, channel labels, etc.).
%   3. Discover and chronologically order raw .ns5 files, tagging each
%      with a task name.
%   4. Loop over files: extract nev/ns5 data and build a trial table,
%      branching on probe type (neuropixel / plexon / FHC / behavior-only).
%   5. If SORTER_HASH is set, attach SpikeInterface sorting output
%      (spike times + cluster metadata) for each probe.
%   6. Merge per-task results into one struct and save to disk.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default paths to add to the MATLAB path, setup for H2P cluster
defaultRAW_PATH   =  '/Volumes/lab_NHPdata';
defaultOUT_PATH   =  '/Volumes/lab_NHPdata-processed';
defaultNEV_PATH   =   '/Users/kendranoneman/Packages/nevutils';
defaultHELP_PATH  =  '/Users/kendranoneman/Projects/mayo/helperfunctions';

%% Parse required/optional input arguments
p = inputParser;
addRequired(p, 'session_name', @ischar);
addParameter(p, 'RAW_DATA_PATH', defaultRAW_PATH, @ischar); 
addParameter(p, 'OUT_DATA_PATH', defaultOUT_PATH, @ischar); 
addParameter(p, 'NEVUTIL_PATH', defaultNEV_PATH, @ischar);
addParameter(p, 'NASNET_PATH', [], @ischar); % only used for plex or FHC
addParameter(p, 'HELPERS_PATH', defaultHELP_PATH, @ischar);
addParameter(p, 'SORTER_HASH', [], @ischar);
addParameter(p, 'BEHAVIOR_ONLY', false, @islogical);
addParameter(p, 'INCLUDE_ITI', true, @islogical);

% Parse inputs
parse(p, session_name, varargin{:});
RAW_PATH       =  p.Results.RAW_DATA_PATH;
OUT_PATH       =  p.Results.OUT_DATA_PATH;
NEV_PATH       =  p.Results.NEVUTIL_PATH;
NET_PATH       =  p.Results.NASNET_PATH;
HELPERS_PATH   =  p.Results.HELPERS_PATH;
SORTER_HASH    =  p.Results.SORTER_HASH;
BEHAV_ONLY     =  p.Results.BEHAVIOR_ONLY;
INCLUDE_ITI    =  p.Results.INCLUDE_ITI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup: add required paths to MATLAB path, create output dir, load metadata
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

addpath(genpath(NEV_PATH));
addpath(fullfile(HELPERS_PATH,'behavior'));
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); % repo root, one level up from this file

if ~exist(fullfile(OUT_PATH, session_name), 'dir'), mkdir(fullfile(OUT_PATH, session_name)); end

session_path = fullfile(RAW_PATH, session_name);
filePattern = fullfile(RAW_PATH, session_name, '*.ns5');

% Load session-level metadata (probe type/config, eye/diode/pupil channel labels, etc.)
S1 = struct();
[metadata, eye_chan_labels, diode_chan_label, pupil_chan_label] = parse_sessionMetadata(session_path);
S1.sess_name = metadata.sess_name;
S1.metadata = metadata;

% Behavior-only sessions (or sessions with no detected probe) skip all
% spike-sorting/neural-network paths
if BEHAV_ONLY | isempty(metadata.probe_type)
    NET_PATH = [];
    SORTER_HASH = [];
    metadata.probe_type = 'behavior';
    disp('BEHAVIOR ONLY');
end

% Normalize the "no sorter" sentinel value to an empty path
if isequal(SORTER_HASH, 'None')
    SORTER_HASH = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parsing and ordering tasks in chronological order
[nevnames, nevpaths, tasks, taskTypes] = parse_and_order_tasks(filePattern);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extracting Ripple data from nev/ns5 datafiles
% Main loop: for each chronologically-ordered .ns5 file, extract nev/ns5
% data and build a per-task trial table (tbl). Behavior differs by probe
% type (neuropixel / plexon / FHC single electrode / behavior-only), so
% each branch below handles its own alignment, spike-sorting, and
% channel-metadata logic before converting to a common table format.
tic

rip_time_start = 0; goodFlag = true;
prev_tempdata = struct([]);
for nevnum = 1:length(nevnames) % loop through nev files, in chronological order
    nevpath = nevpaths{nevnum};
    this_task = tasks{nevnum};

    fprintf('\n---- generating nev_out for %s ----\n', this_task);

    %------------------------------------- NEUROPIXELS -------------------------------------%
    % Aligns ripple (nev/ns5) trial timestamps to imec/nidq sync pulses so
    % that trials can later be mapped onto neuropixel sample indices.
    if ismember('neuropixel',metadata.probe_type)
        % Load in align/sync pulses
        if nevnum==1
            alignCodes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name], [session_name,'_tcat.nidq.bfv_8_0_9.txt']));
            alignTimes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name], [session_name,'_tcat.nidq.bft_8_0_9.txt']));
            alignTimes = alignTimes(alignCodes>0);
        end

        [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false, 'READ_LFP', false, 'alignPulseEnabled', true);
        startAcquisition = datetime(out_ns5.hdr.timeOrigin, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS');

        [dat, dat_iti, ~, ~, ~] = format_datTrials(nev, out_ns5, 'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);

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
            INCLUDE_ITI = false;
        else
            these_alignTimes = alignTimes(np_mask);
            if sum(np_mask) < length(ripple_mask)
                dat = dat(ripple_mask);
                dat_iti = dat_iti(ripple_mask(1:end-1));
            end
        end 

        % Pulling out continuus stream of data
        [eyedata, eye_times, trial_events] = extract_eyeStream(nev, out_ns5, 'EYE_CHAN_LABELS', eye_chan_labels);
        np_win_sec = [these_alignTimes(1), these_alignTimes(1)+(numel(eye_times)/1000)] - 1;
        eye_times = eye_times + (these_alignTimes(1) - alignTimes(1));

        tbl = convert_smithDat_mayoTbl(dat, dat_iti, 'TASK_NAME', this_task, 'HELPERS_PATH', HELPERS_PATH, 'INCLUDE_ITI', INCLUDE_ITI);
    
    %------------------------------------- PLEXON -------------------------------------%
    % Handles multi-electrode plexon arrays; runs nasnet spike-sorting
    % (per hardware config / electrode group) when NET_PATH is set and no
    % SORTER_HASH override is given, otherwise extracts unsorted data.
    elseif ismember('plexon',metadata.probe_type)
        % 'fstm'/'fast' tasks are handled without spike-sorting
        if any(contains(this_task, {'fstm', 'fast'}))
            [nev, out_ns5, ~] = extract_nevout(nevpath);
            if ~isempty(nev)
                [dat, dat_iti, ~, tempdata, ~] = format_datTrials(nev, out_ns5, 'PREV_TEMPDATA', prev_tempdata, ... 
                                        'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);
                prev_tempdata = tempdata;
                tbl = convert_smithDat_mayoTbl(dat, dat_iti, 'TASK_NAME', this_task, 'HELPERS_PATH', HELPERS_PATH, 'INCLUDE_ITI', INCLUDE_ITI);
                tbl = removevars(tbl, 'time_sec');
            end
        else
            % Run nasnet spike-sorting per hardware config (elecA-D map to
            % contiguous 128-channel blocks) and build cluster metadata tables
            if ~isempty(NET_PATH) & isempty(SORTER_HASH)
                addpath(genpath(NET_PATH));

                [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'), 'READ_LFP', false);

                [dat, dat_iti, ~, tempdata, ~] = format_datTrials(nev, out_ns5, 'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);

                sorting_all = [];
                for probe = 1:numel(metadata.hardware_config)
                    sorting.probe_index = probe;
                    % Map hardware config label to its 128-channel block
                    neural_chans = electrode_group_channels(metadata.hardware_config{probe});

                    [dat2, dat_iti, ~, ~, chans] = format_datTrials(nev, out_ns5, 'PREV_TEMPDATA', prev_tempdata, ...
                                    'NEURAL_CHANNELS', neural_chans, 'PROBE_INDEX', probe, ...
                                    'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);

                    % Build per-cluster metadata table (best channel, sort code, probe geometry, etc.)
                    chans_tbl = build_nasnet_clusters_table(chans, metadata, probe);
                    sorting.clusters = categoricalize_columns(chans_tbl);
                    sorting_all = [sorting_all; sorting];

                    fieldsToCopy = {sprintf('spiketimes_%d', probe), sprintf('netlabels_%d', probe)};
                    for f = fieldsToCopy
                        [dat.(f{1})] = deal(dat2.(f{1}));
                    end
                end

                if nevnum==1
                    S1.sorting = sorting_all;
                end
                
                prev_tempdata = tempdata;
                tbl = convert_smithDat_mayoTbl(dat, dat_iti, 'TASK_NAME', this_task, 'HELPERS_PATH', HELPERS_PATH, 'INCLUDE_ITI', INCLUDE_ITI);
                tbl = removevars(tbl, 'time_sec');
            else
                % No nasnet sorting requested (or a SORTER_HASH override will
                % supply sorting later) - just extract unsorted trial data,
                % accumulating a running time offset across files
                [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false);

                [dat, dat_iti, ~, tempdata, ~] = format_datTrials(nev, out_ns5, 'PREV_TEMPDATA', prev_tempdata, ...
                                'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);

                prev_tempdata = tempdata;
                tbl = convert_smithDat_mayoTbl(dat, dat_iti, 'TASK_NAME', this_task, 'HELPERS_PATH', HELPERS_PATH, 'INCLUDE_ITI', INCLUDE_ITI);
                tbl.time_sec = tbl.time_sec+rip_time_start;
                rip_time_start = rip_time_start + double(out_ns5.hdr.nSamples/out_ns5.hdr.Fs);
           end
        end
    %------------------------------------- FHC SINGLE ELECTRODE -------------------------------------%
    % Single-electrode FHC recordings; each file corresponds to one unit,
    % so channel/cluster metadata accumulates across files and is
    % finalized once the last file has been processed.
    elseif ismember('fhc',metadata.probe_type)
        if any(contains(this_task, {'fstm', 'fast'}))
            [nev, out_ns5, ~] = extract_nevout(nevpath);
            if ~isempty(nev)
                [dat, dat_iti, ~, tempdata, ~] = format_datTrials(nev, out_ns5, 'PREV_TEMPDATA', prev_tempdata, ... 
                                        'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);
            end
        else
            addpath(genpath(NET_PATH));

            [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'));

            % Hardware config encodes electrode group + channel (e.g. 'elecA_12');
            % map that to an absolute neural channel index
            hw_config = metadata.hardware_config{1};
            parts = strsplit(hw_config, '_');
            this_hw = parts{1};
            this_chan = str2double(parts{2});
            neural_chan = electrode_group_offset(this_hw) + this_chan;

            [dat, dat_iti, ~, tempdata, chans] = format_datTrials(nev, out_ns5, 'PREV_TEMPDATA', prev_tempdata, ...
                                                  'NEURAL_CHANNELS', neural_chan, 'EYE_CHAN_LABELS', eye_chan_labels, ...
                                                  'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);
            probe = 1;
            if nevnum==1
                sorting.probe_index = probe;
                S1.sorting = sorting;
                chans_tbl_all = [];
            end

            chans_tbl = build_fhc_clusters_row(this_task, chans, metadata, probe);
            chans_tbl_all = [chans_tbl_all; chans_tbl];

        end   
 
        % On the last file, deduplicate and finalize the accumulated
        % per-unit cluster metadata table into S1.sorting.clusters
        if nevnum==length(nevnames)
            S1.sorting.clusters = finalize_fhc_clusters_table(chans_tbl_all);
        end

        prev_tempdata = tempdata;
        [~, fname, ~] = fileparts(nevpath);   % 'kendra_scrappy_0066a_mdir1'
        this_task = erase(fname, session_name); % 'a_mdir1'

        tbl = convert_smithDat_mayoTbl(dat, dat_iti, 'TASK_NAME', this_task, 'HELPERS_PATH', HELPERS_PATH, 'INCLUDE_ITI', INCLUDE_ITI);
        tbl = removevars(tbl, 'time_sec');

    %------------------------------------- BEHAVIOR ONLY -------------------------------------%
    % No probe/neural data available (or BEHAVIOR_ONLY was requested);
    % extract behavior-only trial data with no spike-sorting.
    else
        [nev, out_ns5, ~] = extract_nevout(nevpath);
        if ~isempty(nev)
            [dat, dat_iti, ~, ~, ~] = format_datTrials(nev, out_ns5, 'EYE_CHAN_LABELS', eye_chan_labels, 'DIODE_CHAN_LABEL', diode_chan_label, 'PUPIL_CHAN_LABEL', pupil_chan_label);
            tbl = convert_smithDat_mayoTbl(dat, dat_iti, 'TASK_NAME', this_task, 'HELPERS_PATH', HELPERS_PATH, 'INCLUDE_ITI', INCLUDE_ITI);
            tbl = removevars(tbl, 'time_sec');
        else
            dat = []; tbl = [];
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Extracting sorting data from SpikeInterface outputs (e.g., kilosort)
    % When SORTER_HASH is provided, attach externally-computed spike sorting
    % (per probe) to the current task's trial table, converting sorted spike
    % times into per-trial cell arrays and building/loading cluster metadata.
    if SORTER_HASH
        % SORTER_HASH is a hyphen-delimited hash string identifying the
        % preprocessing / motion-correction / sorter run to pull from
        hashes = split(SORTER_HASH, '-');
        preprocess_hash = hashes{1}; pp_hash = hashes{2}; motion_hash = hashes{3}; %sorter_hash = hashes{4};

        sorting_all = []; spike_times = cell(numel(metadata.probe_type),1);
        for probe = 1:numel(metadata.probe_type)

            %si_path = fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], 'sorting', SORTER_HASH);
            si_path = fullfile(OUT_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], 'sorting', SORTER_HASH);

            if isequal(metadata.probe_type{probe},'plexon')
                ripple_info = loadMetadataJSON(fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], 'ripple_info.json'));

                % Pull out SpikeInterface sorting outputs
                [spikes_perTrial, sorting, ~] = parse_SortingToTbl(tbl, fullfile(si_path, 'sorter_output'), 'Fs', ripple_info.Fs);

            elseif isequal(metadata.probe_type{probe},'neuropixel')
                lfp_meta = readMetaFile(fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], [session_name, '_t0.', metadata.hardware_config{probe}, '.lf.meta']));
                ap_meta = readMetaFile(fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], [session_name, '_t0.', metadata.hardware_config{probe}, '.ap.meta']));

                lfp_fs = lfp_meta.imSampRate;

                % Reconstruct per-trial start/end times (in seconds and LFP
                % samples) from the neuropixel sync pulses for this probe
                if INCLUDE_ITI
                    new_alignTimes = repelem(these_alignTimes, 2); new_alignTimes(end) = [];
                else
                    new_alignTimes = these_alignTimes;
                end
                
                trial_starts_sec = cellfun(@(q,v) q-(v(1)./1000), num2cell(new_alignTimes), tbl.ALIGN_PULSE, 'uni', 1);
                trial_ends_sec = trial_starts_sec + (tbl.time_sec(:,2)-tbl.time_sec(:,1));
    
                trial_starts_lfp_samp = floor(trial_starts_sec * lfp_fs);
                trial_ends_lfp_samp = ceil(trial_ends_sec * lfp_fs);
    
                tbl.imec_sec = [trial_starts_sec trial_ends_sec];
                tbl.imecLFP_samp = [trial_starts_lfp_samp trial_ends_lfp_samp];

                % Pull out SpikeInterface sorting outputs
                [spikes_perTrial,sorting,~] = parse_SortingToTbl(tbl, fullfile(si_path,'sorter_output'), 'NP_ALIGN_PULSES', new_alignTimes, 'Fs', ap_meta.imSampRate);

                % Continuous spikes to match eyedata stream
                sptimes = extract_spikesInWindow(fullfile(si_path,'sorter_output'), np_win_sec, 'Fs', 30000);
                spike_times{probe} = cellfun(@(q) q + (these_alignTimes(1) - alignTimes(1)), sptimes, 'uni', 0);
            end

            tbl.(sprintf('spiketimes_%d',probe)) = spikes_perTrial;  

            % Cluster metadata is only assembled once (on the first file),
            % since it is shared across all tasks/files for this probe
            if nevnum==1
                sorting.probe_index = probe;
                fields = fieldnames(sorting);
                fields(strcmp(fields, 'probe_index')) = [];
                sorting = orderfields(sorting, ['probe_index'; fields]);

                % Prefer precomputed quality metrics if available, else build minimal cluster metadata
                if isfile(fullfile(si_path,'quality_metrics','cluster_metrics.csv'))
                    metrics = parse_clusterMetrics(si_path);
                    sorting.clusters = [metrics removevars(sorting.clusters, 'cluster_id')];
                else
                    sorting.clusters = build_default_probe_clusters(sorting.clusters, metadata, probe);
                end
                sorting.clusters.probe_index = repmat(probe,height(sorting.clusters),1);
                sorting.clusters = movevars(sorting.clusters,{'probe_index'},'After','sess_name');
                sorting.clusters = categoricalize_columns(sorting.clusters);

                motion_path = fullfile(RAW_PATH, session_name, [session_name, '_', metadata.hardware_config{probe}], 'preprocess', preprocess_hash, pp_hash, motion_hash);

                sorting_all = [sorting_all; sorting];
            end
        end

        if nevnum==1
            jsonStr = fileread(fullfile(si_path,'params.json'));
            protocolStruct = jsondecode(jsonStr);
            S1.protocol = protocolStruct;
            S1.sorting = sorting_all;
        end

    end

    % Store this file's results into S1, keyed by task name: task
    % parameters (merged across trials), the ns5 header, raw dat struct,
    % and the final trial table
    if ~isempty(tbl)
        merged_struct = merge_taskParams(tbl);
        S1.(this_task).params = merged_struct;
    end
    S1.(this_task).hdr = out_ns5.hdr;

    S1.(this_task).dat = dat;
    S1.(this_task).tbl = tbl;

    if ismember('neuropixel',metadata.probe_type)
        S1.(this_task).eyedata = eyedata;
        S1.(this_task).eye_times = eye_times;
        S1.(this_task).trial_events = trial_events;
        S1.(this_task).spike_times = spike_times;
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Final steps and saving table
% Merge all per-task tables/metadata in S1 into a single unified struct S,
% then save it with a filename that reflects which sorting (if any) was used:
%   <session>-<short_sorter_id>.mat  when SORTER_HASH is set
%   <session>-nasnet.mat             when NET_PATH (nasnet) sorting was used
%   <session>.mat                    otherwise (unsorted / behavior-only)
S = unify_taskTables(S1,taskTypes);

if ~exist(fullfile(OUT_PATH, session_name, 'tables'), 'dir'), mkdir(fullfile(OUT_PATH, session_name, 'tables')); end

if SORTER_HASH
    short_id = shorten_sorter_path(SORTER_HASH, fullfile(OUT_PATH, session_name, 'tables'));
    save(fullfile(OUT_PATH, session_name, 'tables', sprintf('%s-%s.mat', session_name, short_id)), 'S', '-v7.3');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nevnames, nevpaths, tasks, taskTypes] = parse_and_order_tasks(filePattern)
%PARSE_AND_ORDER_TASKS  Discover .ns5 files, sort chronologically, and tag each with a task name.
%   [nevnames, nevpaths, tasks, taskTypes] = parse_and_order_tasks(filePattern)
%
%   INPUT
%     filePattern : dir()-style glob for the session's .ns5 files
%
%   OUTPUT
%     nevnames  : file names without extension, in chronological order
%     nevpaths  : full paths to each file, same order as nevnames
%     tasks     : per-file task tag (e.g. 'fstm01', 'unit03', or 'unknown')
%     taskTypes : unique task types across all files (e.g. 'fstm', 'unit')

% Create the search pattern to find files that start with 'filename' and end with '.ns5'
raw_files = dir(filePattern);
raw_filenames = {raw_files.name}.';
nevnames = cellfun(@(q) q(1:end-4), raw_filenames, 'uni', 0);
raw_filepaths = arrayfun(@(x) fullfile(x.folder, x.name), raw_files, 'uni', 0);

% Sort files chronologically by their recording start time (not filename)
nsx_hdrs = cellfun(@(q) read_nsx(q,'readdata',false), raw_filepaths, 'uni', 0);
recording_times = cellfun(@(l) l.hdr.timeOrigin, nsx_hdrs, 'uni', 0);
[~,idx] = sort(recording_times);
nevnames = nevnames(idx);
nevpaths = raw_filepaths(idx);

% Define possible task keywords
task_keywords = handle_taskSpecifics();

% Tag each .ns5 file with a task name by matching against known task
% keywords (e.g. 'fstm01', 'unit03'); files that match nothing are 'unknown'
tasks = cell(size(nevnames));
for i = 1:numel(nevnames)
    tasks{i} = match_task_keyword(nevnames{i}, task_keywords);
end

% Reduce tagged task names (e.g. 'fstm01') down to their task type (e.g. 'fstm')
taskTypes = unique(cellfun(@(q) regexp(q, '[a-zA-Z]+', 'match', 'once'), tasks, 'uni', 0));
disp(tasks)

end

function task = match_task_keyword(name, task_keywords)
%MATCH_TASK_KEYWORD  Return the first task keyword (+ trailing chars) found in name, or 'unknown'.
task = 'unknown';
for j = 1:numel(task_keywords)
    pattern = [task_keywords{j}, '\w*'];
    match = regexp(name, pattern, 'match', 'once');
    if ~isempty(match)
        task = match;
        return
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function neural_chans = electrode_group_channels(hw_label)
%ELECTRODE_GROUP_CHANNELS  Map a plexon hardware-config label ('elecA'-'elecD') to its channel block.
if isequal(hw_label,"elecA")
    neural_chans = 1:128;
elseif isequal(hw_label,"elecB")
    neural_chans = 129:256;
elseif isequal(hw_label,"elecC")
    neural_chans = 257:384;
elseif isequal(hw_label,"elecD")
    neural_chans = 385:513;
end
end

function offset = electrode_group_offset(hw_label)
%ELECTRODE_GROUP_OFFSET  Map an FHC hardware-config group label ('elecA'-'elecD') to its channel offset.
if isequal(hw_label,"elecA")
    offset = 0;
elseif isequal(hw_label,"elecB")
    offset = 128;
elseif isequal(hw_label,"elecC")
    offset = 256;
elseif isequal(hw_label,"elecD")
    offset = 384;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tbl = categoricalize_columns(tbl)
%CATEGORICALIZE_COLUMNS  Convert any char/cellstr columns of tbl to categorical, in place.
vars = tbl.Properties.VariableNames;
for v = vars
    col = tbl.(v{1});
    if ischar(col) || (iscell(col) && all(cellfun(@ischar, col)))
        tbl.(v{1}) = categorical(string(col));
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chans_tbl = build_nasnet_clusters_table(chans, metadata, probe)
%BUILD_NASNET_CLUSTERS_TABLE  Build per-cluster metadata table for one nasnet-sorted plexon probe.
chans_tbl = table(repmat(metadata.sess_name, size(chans,1), 1), 'VariableNames', {'sess_name'});
chans_tbl.probe_index = repmat(probe, height(chans_tbl), 1);
chans_tbl.cluster_id = (1:height(chans_tbl))'-1;
chans_tbl.best_channel = chans(:,1);
chans_tbl.sort_code = chans(:,2);
chans_tbl.probe_label = repmat(metadata.probe_label{probe}, height(chans_tbl), 1);
chans_tbl.probe_type = repmat(metadata.probe_type{probe}, height(chans_tbl), 1);
chans_tbl.probe_config = repmat(metadata.probe_config{probe}, height(chans_tbl), 1);
chans_tbl.hardware_config = repmat(metadata.hardware_config{probe}, height(chans_tbl), 1);
chans_tbl.probe_depth_mm = repmat(metadata.probe_depth_mm(probe), height(chans_tbl), 1);
chans_tbl.probe_gridHole = repmat(metadata.probe_gridHole(probe), height(chans_tbl), 1);
end

function chans_tbl = build_fhc_clusters_row(this_task, chans, metadata, probe)
%BUILD_FHC_CLUSTERS_ROW  Build one file's per-unit cluster metadata row(s) for an FHC probe.
% Task name encodes the unit number (e.g. 'unit03') being recorded in this file
out = regexp(this_task, 'unit(\d{2})', 'tokens', 'once');
unitNum = str2double(out{1});

if size(chans,1) == 1
    chans_tbl = table({metadata.sess_name}, 'VariableNames', {'sess_name'});
else
    chans_tbl = table(repmat(metadata.sess_name, size(chans,1), 1), 'VariableNames', {'sess_name'});
end

chans_tbl.sess_name = categorical(string(chans_tbl.sess_name));
chans_tbl.probe_index = repmat(probe, height(chans_tbl), 1);
chans_tbl.unit_id = repmat(unitNum, height(chans_tbl), 1);
chans_tbl.best_channel = chans(:,1);
chans_tbl.sort_code = chans(:,2);
if unitNum <= numel(metadata.probe_depth_mm)
    chans_tbl.probe_depth_mm = repmat(metadata.probe_depth_mm(unitNum), height(chans_tbl), 1);
else
    chans_tbl.probe_depth_mm = NaN(height(chans_tbl), 1);
end
chans_tbl.probe_label = repmat(metadata.probe_label{probe}, height(chans_tbl), 1);
chans_tbl.probe_type = repmat(metadata.probe_type{probe}, height(chans_tbl), 1);
chans_tbl.probe_config = repmat(metadata.probe_config{probe}, height(chans_tbl), 1);
chans_tbl.hardware_config = repmat(metadata.hardware_config{probe}, height(chans_tbl), 1);
chans_tbl.probe_gridHole = repmat(metadata.probe_gridHole{probe}, height(chans_tbl), 1);
end

function chans_tbl_all = finalize_fhc_clusters_table(chans_tbl_all)
%FINALIZE_FHC_CLUSTERS_TABLE  Deduplicate rows, assign cluster ids, and categoricalize the accumulated FHC cluster table.
chans_tbl_all = unique(chans_tbl_all, 'rows');
chans_tbl_all.cluster_id = (1:height(chans_tbl_all))'-1;
chans_tbl_all = movevars(chans_tbl_all,{'cluster_id'},'After','probe_index');
chans_tbl_all = categoricalize_columns(chans_tbl_all);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clusters = build_default_probe_clusters(clusters, metadata, probe)
%BUILD_DEFAULT_PROBE_CLUSTERS  Fill in session/probe metadata columns when no quality-metrics CSV is available.
clusters.sess_name = repmat(metadata.sess_name,height(clusters),1);
clusters.monkey = repmat(metadata.monkey,height(clusters),1);
clusters.experimenter = repmat(metadata.experimenter,height(clusters),1);
clusters.probe_id = repmat(probe-1,height(clusters),1);
clusters.probe_index = repmat(probe,height(clusters),1);
clusters.probe_label = repmat(metadata.probe_label{probe},height(clusters),1);
clusters.probe_type = repmat(metadata.probe_type{probe},height(clusters),1);
clusters.probe_config = repmat(metadata.probe_config{probe},height(clusters),1);
clusters.hardware_config = repmat(metadata.hardware_config{probe},height(clusters),1);
clusters.probe_depth_mm = repmat(metadata.probe_depth_mm(probe),height(clusters),1);
clusters.probe_gridHole = repmat(metadata.probe_gridHole(probe),height(clusters),1);
end
