function [dat, dat_iti, epochEnd, tempdata, channels] = format_datTrials(nev, out_ns5, varargin)
    % format_datTrials - Processes neural and behavioral data for multiple trials, 
    % extracts eye and spike data, and formats the information into a structured array.
    %
    % This function processes and organizes trial-based data from neural event (nev) and neural signal (ns5) files. 
    % It extracts relevant trial-specific information such as eye position, velocity, and pupil data, as well as spike 
    % times for specified neural channels. The function is capable of handling various data epochs and organizes the data 
    % into a structured format that includes trial codes, trial results, parameters, and neural spikes.
    %
    % The function also supports neural spike sorting and handles missing or corrupted trial start/end codes. 
    % The data is downsampled and formatted for further analysis.
    %
    %%%% Required inputs: %%%
    %   nev1               -   nev file containing event data (e.g., neural threshold crossings and digital codes)
    %   out_ns5            -   ns5 file containing the raw 30kHz data (e.g., eye data, pupil, diode, raw neural signals)

    %%%% Optional parameters: %%%
    %   NEURAL_CHANNELS  -  Array of channel IDs to extract neural data from. Leave empty if only behavioral data 
    %                        is needed. Default is an empty array (i.e., []), meaning no neural data is extracted.
    %   EYE_CHAN_LABELS  -  Cell array specifying the eye movement channels. Default is {'10241', '10242'}, 
    %                        representing typical labels for {Eye_HE, Eye_VE}.
    %   DIODE_CHAN_LABEL -  Char of label for photodiode, if recorded. Default is '10243'.
    %
    %   PUPIL_CHAN_LABEL -  Char of label for pupil, if recorded. Default is '10244'.
    %
    %%%% Outputs: %%%
    %   dat      -    Array of structs with the following fields:
    %       - block: The block number of the trial
    %       - time: Start and end times of the trial in seconds
    %       - text: Trial-specific text or annotations
    %       - trialcodes: Trial-specific codes for event markers
    %       - result: The result of the trial, such as "correct" or "incorrect"
    %       - params: The parameters associated with the trial (e.g., reaction times, crossing times)
    %       - eyes: Processed eye position data (in degrees)
    %       - pupil: Pupil data
    %       - diode: Diode signal data
    %       - spiketimes: Spike times for the specified neural channels, if applicable
    %       - net_labels: Spike sorting labels, if spike sorting is enabled
    %
    %%%% Example usage: %%%
    %   [dat,~] = format_datTrials(nev1, out_ns5, 'NEURAL_CHANNELS', [0,1,2,3,4,5,6,7,8,9])
    %
    % This will process the nev1 and out_ns5 files, extract eye and neural data, 
    % and return the data in the structured array `dat_all`, with separate trial information for each epoch.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultEyeChanLabels = {'10241', '10242'};

    % Create an input parser
    p = inputParser;
    addRequired(p, 'nev', @(x) (isnumeric(x)) || isstruct(x));
    addRequired(p, 'out_ns5', @isstruct);
    addParameter(p, 'NEURAL_CHANNELS', [], @isnumeric);
    addParameter(p, 'PROBE_INDEX', 1, @isnumeric);
    addParameter(p, 'EYE_CHAN_LABELS', defaultEyeChanLabels, (@(x) iscell(x))); % eye channel labels
    addParameter(p, 'DIODE_CHAN_LABEL', '10243', @ischar); % diode channel label
    addParameter(p, 'PUPIL_CHAN_LABEL', '10244', @ischar); % pupil channel label
    addParameter(p, 'PREV_TEMPDATA', [], @isstruct);

    % Parse the inputs
    parse(p, nev, out_ns5, varargin{:});

    % Assign parsed values to variables
    nev = p.Results.nev;
    out_ns5 = p.Results.out_ns5;
    neural_channels = p.Results.NEURAL_CHANNELS;
    probe_index = p.Results.PROBE_INDEX;
    eye_channel_labels = p.Results.EYE_CHAN_LABELS;
    pupil_channel_label = p.Results.PUPIL_CHAN_LABEL;
    diode_channel_label = p.Results.DIODE_CHAN_LABEL;
    prev_tempdata = p.Results.PREV_TEMPDATA;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Fs = double(out_ns5.hdr.Fs); % Sampling frequency
    nEpochs = size(out_ns5.hdr.timeStamps, 2); % Number of epochs
    
    starttrial = 1;
    endtrial = 255; % Trial end code
    
    % Identify eye movement channels (Eye_HE, Eye_VE, DIODE, PUPIL)
    eye_channels = find(ismember(out_ns5.hdr.label, eye_channel_labels));
    pupil_channel = find(ismember(out_ns5.hdr.label, pupil_channel_label));
    diode_channel = find(ismember(out_ns5.hdr.label, diode_channel_label));
    
    % Determine if nev is an array of struct; if so, extract nev data
    if isequal(class(nev), 'struct')
        NEV = [nev.nev nev.net_labels']; % nev.waveforms]; % Combine event and neural network labels
        spike_sort = true; % Flag for spike sorting
    else
        NEV = nev;
        spike_sort = false; % No spike sorting if nev is not a struct
    end
    
    dat_all = []; dat_iti_all = []; % Initialize output structure arrays
    past_epochEnd = 0; % Keep track of the end of the previous epoch
    block = 1; % Block number for trial grouping
    
    % Loop through each epoch
    for epoch = 1:nEpochs
        % Extract data for the current epoch
        epochStart = out_ns5.hdr.timeStamps(1, epoch); % Epoch start time (samples)
        epochEnd = out_ns5.hdr.timeStamps(2, epoch); % Epoch end time (samples)
    
        nsStartTime = double(epochStart / Fs); % Convert to seconds
        nsEndTime = double(epochEnd / Fs) + 0.3; % Add a small buffer for the end time
    
        epochDiff = epochEnd - epochStart;
        epochStart_samp = past_epochEnd + 1;
        epochEnd_samp = (epochStart_samp + epochDiff) - 1;
    
        % Extract trial data within the current epoch time range
        this_nev = NEV(NEV(:, 3) >= nsStartTime & NEV(:, 3) <= nsEndTime, :);
        ns5_rng = epochStart_samp:epochEnd_samp; % Range of samples in the epoch

        % Extract digital event codes
        diginnevind = find(this_nev(:, 1) == 0);
        digcodes = this_nev(diginnevind, :);

        channels = unique(NEV(NEV(:,1) ~= 0,1:2),'rows');
        channels = channels(channels(:,1) ~= 0 & ismember(channels(:,1), neural_channels),:);
    
        % Find trial start and end indices
        trialstartindstemp = find(digcodes(:, 2) == starttrial);
        trialstartinds = diginnevind(trialstartindstemp);
        trialstarts = this_nev(trialstartinds, 3); % Trial start times

        trialendindstemp = find(digcodes(:, 2) == endtrial);
        trialendinds = diginnevind(trialendindstemp);
        trialends = this_nev(trialendinds, 3); % Trial end times
    
        % Detect missing start/end codes and handle missing data
        [trialstarts, trialends, trialstartgood, trialendgood] = detectMissingStartEndCode(trialstarts, trialends);
        trialstartinds = trialstartinds(trialstartgood);
        trialendinds = trialendinds(trialendgood);
    
        % Handle mismatched start/end times
        if length(trialstarts) ~= length(trialends) || sum((trialends - trialstarts) < 0)
            if sum(trialstarts(1:end-1) >= trialends) == 0
                trialstarts = trialstarts(1:end-1); % Trim the last trial start
            end
        end
    
        % Convert trial start and end times to samples
        trialstarts_samp = round(trialstarts * Fs) - epochStart;
        trialends_samp = round(trialends * Fs) - epochStart;
        past_epochEnd = epochEnd_samp;
    
        % Get session initial parameters
        predatcodes = digcodes(digcodes(:, 3) < trialstarts(1), :);
        tempdata.text = char(predatcodes(predatcodes(:, 2) >= 256 & predatcodes(:, 2) < 512, 2) - 256)';
        if isempty(tempdata.text)
            tempdata = prev_tempdata;
        else
            tempdata = getDatParams(tempdata); % Get parameters associated with the data
        end
            
        [dat, dat_iti] = init_trial_structs(length(trialstarts), ~isempty(eye_channel_labels));

        dat_iti(end) = [];

        %% Loop through trials and organize data
        for n = 1:length(trialstarts)
            if mod(n, 100) == 0
                fprintf('Processed nev for %i trials of %i...\n', n, length(trialstarts));
            end
            dat(n).block = block;
            dat(n).channels = channels;
            dat(n).time = [trialstarts(n) trialends(n)]; % Store trial start and end times
            this_trial = this_nev(trialstartinds(n):trialendinds(n), :);
            trialdig = this_trial(this_trial(:, 1) == 0, :);
            dat(n).text = char(trialdig(trialdig(:, 2) >= 256 & trialdig(:, 2) < 512, 2) - 256)';
            dat(n).trialcodes = trialdig(trialdig(:, 2) < 256 | (trialdig(:, 2) >= 1000 & trialdig(:, 2) <= 32000), 1:3);

            % Extract trial result
            event = uint32(trialdig);
            dat(n).result = event(event(:, 2) >= 160 & event(:, 2) <= 165, 2);
            if isempty(dat(n).result)
                dat(n).result = event(event(:, 2) >= 150 & event(:, 2) <= 158, 2);
            end
            if isempty(dat(n).result)
                dat(n).result = NaN;
            end

            % Extract block parameters
            blockParams = tempdata.params.trial;
            if isfield(blockParams, 'reactionTime')
                blockParams.crossingTime = blockParams.reactionTime; % Copy reaction time to crossing time
                blockParams = rmfield(blockParams, 'reactionTime'); % Remove old field
            end
                        
            dat(n).params.block = blockParams;
            
            % Check if next trial is in a new block
            if n < length(trialstarts) && trialstartinds(n+1) - trialendinds(n) > 1
                bt = this_nev(trialendinds(n) + 1:trialstartinds(n + 1) - 1, :);
                btdig = bt(bt(:, 1) == 0, :);
                if sum(find(btdig(:, 2) >= 256 & btdig(:, 2) < 512)) > 0
                    tempdata.text = char(btdig(btdig(:, 2) >= 256 & btdig(:, 2) < 512, 2) - 256)';
                    tempdata = getDatParams(tempdata);
                    block = block + 1; % Move to the next block
                end
            end

            %dat(n).ns5_samps = ns5_rng([trialstarts_samp(n),trialends_samp(n)]);

            % Extract and process eye data
            if ~isempty(eye_channel_labels)
                [dat(n).eyedata, dat(n).pupil, dat(n).diode] = extract_eye_pupil_diode(out_ns5, ...
                    ns5_rng(trialstarts_samp(n):trialends_samp(n)), dat(n).params, eye_channels, pupil_channel, diode_channel);
            end
            
            % Process neural spikes, if applicable
            if ~isempty(neural_channels)
                [spks_byUnit, netLabels_byUnit] = extract_spikes_byUnit(this_trial, neural_channels, channels, trialstarts(n), spike_sort);
                dat(n).(sprintf('spiketimes_%d', probe_index)) = spks_byUnit; % Store spike times
                if spike_sort
                    dat(n).(sprintf('netlabels_%d', probe_index)) = netLabels_byUnit; % Store spike sorting labels
                end
            end

            % INTERTRIAL-INTERVALS
            if n < numel(trialstarts)
                dat_iti(n).block = block;
                dat_iti(n).channels = channels;
                dat_iti(n).time = [trialends(n)+(1/Fs) trialstarts(n+1)-(1/Fs)]; % Store trial start and end times
                dat_iti(n).text = dat(n).text;
                this_trial = this_nev(trialendinds(n)+1:trialstartinds(n+1)-1, :);
                dat_iti(n).trialcodes = [0, 1, trialends(n)+(1/Fs); 0, 255, trialstarts(n+1)-(1/Fs)];
                dat_iti(n).result = NaN;
                dat_iti(n).params.block = dat(n).params.block;
    
                % Extract and process eye data
                if ~isempty(eye_channel_labels)
                    [dat_iti(n).eyedata, dat_iti(n).pupil, dat_iti(n).diode] = extract_eye_pupil_diode(out_ns5, ...
                        ns5_rng(trialends_samp(n)+1:trialstarts_samp(n+1)-1), dat_iti(n).params, eye_channels, pupil_channel, diode_channel);
                end
                
                % Process neural spikes, if applicable
                if ~isempty(neural_channels)
                    [spks_byUnit, netLabels_byUnit] = extract_spikes_byUnit(this_trial, neural_channels, channels, trialends(n)+(1/Fs), spike_sort);
                    dat_iti(n).(sprintf('spiketimes_%d', probe_index)) = spks_byUnit; % Store spike times
                    if spike_sort
                        dat_iti(n).(sprintf('netlabels_%d', probe_index)) = netLabels_byUnit; % Store spike sorting labels
                    end
                end
            end

    
        end

        dat = getDatParams(dat); % Final parameter extraction
        dat_iti = getDatParams(dat_iti);
    
        % Concatenate the new structured data to the main array
        dat_all = [dat_all; dat];
        dat_iti_all = [dat_iti_all; dat_iti];
    end

    dat = dat_all;
    dat_iti = dat_iti_all;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dat, dat_iti] = init_trial_structs(n_trials, has_eye_channels)
%INIT_TRIAL_STRUCTS  Preallocate the per-epoch dat/dat_iti struct arrays.
%   Field list matches the original trial/eye/pupil/diode struct exactly;
%   eyedata/pupil/diode fields are only added when eye channels are present.
template = struct('block', [], 'channels', [], 'time', [], 'text', '', ...
    'trialcodes', [], 'result', NaN, 'params', struct());
if has_eye_channels
    template.eyedata = [];
    template.pupil = [];
    template.diode = [];
end
[dat, dat_iti] = deal(repmat(template, n_trials, 1));
end

function [eyedata, pupil, diode] = extract_eye_pupil_diode(out_ns5, samp_range, params, eye_channels, pupil_channel, diode_channel)
%EXTRACT_EYE_PUPIL_DIODE  Extract and downsample eye/pupil/diode signals over one sample range.
eyes = out_ns5.data(eye_channels, samp_range);
eyes_1khz = downsample(eyes', 30)'; % Downsample to 1 kHz
[eyedata, ~] = eye2deg(eyes_1khz(1:2, :), params); % Convert to degrees

pupil = [];
if ~isempty(pupil_channel)
    p = out_ns5.data(pupil_channel, samp_range);
    pupil = downsample(p', 30)'; % Downsample to 1 kHz
end

diode = [];
if ~isempty(diode_channel)
    d = out_ns5.data(diode_channel, samp_range);
    diode = downsample(d', 30)'; % Downsample to 1 kHz
end
end

function [spks_byUnit, netLabels_byUnit] = extract_spikes_byUnit(this_trial, neural_channels, channels, ref_time, spike_sort)
%EXTRACT_SPIKES_BYUNIT  Extract per-channel spike times (ms, relative to ref_time) and sort labels.
spks = this_trial(ismember(this_trial(:, 1), neural_channels), :);
[spks_byUnit, netLabels_byUnit] = deal(cell(1, size(channels, 1)));
for u = 1:size(channels,1)
    spks_byUnit{u} = ((spks(spks(:,1) == channels(u,1) & spks(:,2) == channels(u,2), 3)') - ref_time) .* 1000;
    if spike_sort
        netLabels_byUnit{u} = (spks(spks(:,1) == channels(u,1) & spks(:,2) == channels(u,2), 4)');
    end
end
end
