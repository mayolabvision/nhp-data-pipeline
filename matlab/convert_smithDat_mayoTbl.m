function tbl = convert_smithDat_mayoTbl(dat,varargin)
    % format_dataTable - Processes neural and behavioral data from nev and ns5 files into a structured table.
    %
    % This function processes and formats data from neural and behavioral recording files. It takes in 
    % neural event (nev) and neural data (ns5) files, extracts relevant data based on task parameters, 
    % and returns a structured table with trial information, conditions, and eye movement data.
    %
    % The function supports optional parameters to customize the extraction of neural data, eye movement 
    % signals, and task-specific data. It can also format the output based on specific task types (e.g., 'rfmp', 
    % 'fstm', 'mdir', 'purs') and includes the option to organize data in a table format (because structs suck). 
    %
    %%%% Required inputs: %%%
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
    %%%% Optional parameters: %%%
    %   TASK_NAME   -  Name of the task, if you want to label it something other than what you named is in Ex ('sessionNumber').
    %
    %   LFP         -  Default = empty, if contains something should be cell array with same height as number of trials 
    %
    %%%% Outputs: %%%
    %   tbl  -  A table with many columns, including but not limited to:
    %           - trialName: Name or identifier of the trial
    %           - conditions: Conditions associated with the trial (e.g., jump, pursuitSpeed), with NaN if the condition is missing
    %           - result: Outcome of the trial, in categorical format (e.g., 'correct', 'incorrect')
    %           - params: Parameters of the trial and block (e.g., screen distance, target duration)
    %           - eyePos, eyeVel, eyeAcc: Filtered and smoothed eye movement data (position, velocity, and acceleration)
    %           - eyePos_raw: Raw, unfiltered eye data
    %           - spiketimes: Spike times for the specified neural channels, if applicable
    %           - net_labels: Spike sorting labels, if spike sorting is enabled
    %
    %%%% Example usage: %%%
    %   tbl = format_dataTable(nev, out_ns5, 'CONVERT_TO_TABLE', true)
    %
    % This will process the nev and ns5 files, with data organized in a table and task-specific data converted accordingly.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultHELP_PATH  =  '/Users/kendranoneman/Projects/mayo/helperfunctions';

    % Create an input parser
    p = inputParser;
    addRequired(p, 'dat', @(x) (isnumeric(x)) || isstruct(x));
    addParameter(p, 'TASK_NAME', [], @ischar)
    addParameter(p, 'LFP', []);
    addParameter(p, 'HELPERS_PATH', defaultHELP_PATH, @ischar)

    % Parse the inputs
    parse(p, dat, varargin{:});

    % Assign parsed values to variables
    dat = p.Results.dat;
    TASK_NAME = p.Results.TASK_NAME;
    LFP = p.Results.LFP;
    HELPERS_PATH = p.Results.HELPERS_PATH;

    addpath(fullfile(HELPERS_PATH,'behavior'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(TASK_NAME)
        if isnumeric(dat(1).params.block.sessionNumber)
            TASK_NAME = sprintf('s%d',dat(1).params.block.sessionNumber);
        else
            TASK_NAME = dat(1).params.block.sessionNumber;
        end
    end

    tbl1 = struct2table(dat);
    tbl = table();

    % re-arranging table to be easier to access data 
    tbl.trialName = cellfun(@(q) [TASK_NAME, '.', sprintf('%04d', q)], num2cell(1:height(tbl1))', 'uni', 0);
    tbl.trialName = categorical(string(tbl.trialName));
    if ismember('block', tbl1.Properties.VariableNames), tbl.block = tbl1.block; end
    if ismember('time', tbl1.Properties.VariableNames), tbl.time_sec = tbl1.time; end

    if iscell(tbl1.result)
        tbl.result = convertBetween_eventCodes_eventNames(tbl1.result);
    else
        tbl.result = convertBetween_eventCodes_eventNames(num2cell(tbl1.result));
    end
    tbl.result = categorical(string(tbl.result));

    % Make array of times of start/end time per trial, for aligning with trial codes and indexing eye data
    times_ms = cellfun(@(q) round(q(1)*1000:(q(2)+1)*1000), num2cell(tbl1.time,2), 'uni', 0);
    trialStarts = cellfun(@(q,r) find(q == round(r(r(:,2)==1,3)*1000)), times_ms, tbl1.trialcodes, 'uni', 0);
    eventCodes = tbl1.trialcodes; eventCodes = vertcat(eventCodes{:});
    eventCodes = num2cell(sort(unique(eventCodes(:,2))))';

    eventNames = convertBetween_eventCodes_eventNames(eventCodes);
    trialMarkers = cellfun(@(t) cellfun(@(q,r) find(ismember(q,round(r(r(:,2)==t,3)*1000))), times_ms, tbl1.trialcodes, 'uni', 0), eventCodes, 'uni', 0)';
    trialMarkers = horzcat(trialMarkers{:});
    trialMarkers(cellfun('isempty',trialMarkers)) = {NaN};

    for m=1:length(eventNames)
        if sum(cellfun(@(q) size(q,1)>1, trialMarkers(:,m), 'uni', 1))>0 || sum(cellfun(@(q) size(q,2)>1, trialMarkers(:,m), 'uni', 1))>0
            if numel(unique(cellfun(@numel ,trialMarkers(:,m),'uni',1))) == 1
                codeSplt = cellfun(@(q) num2cell(q), trialMarkers(:,m), 'uni', 0);
                tbl.(eventNames{m}) = vertcat(codeSplt{:});
            else
                tbl.(eventNames{m}) = trialMarkers(:,m);
            end
        else
            tbl.(eventNames{m}) = cell2mat(trialMarkers(:,m));
        end
    end

    tbl.text = tbl1.text;
    tbl.params = tbl1.params;
    tbl.eyedata = tbl1.eyedata; tbl.pupil = tbl1.pupil; tbl.diode = tbl1.diode;

    try
        [~, tbl] = handle_taskSpecifics(tbl, TASK_NAME);
    catch
        disp('------------- task-specific additions failed -------------');
    end

    %tbl.ns5_samps = tbl1.ns5_samps;

    %% Add spiking data if nasnet was used 
    % Find all variable names in tbl1 that start with 'spiketimes_'
    spike_vars = tbl1.Properties.VariableNames(startsWith(tbl1.Properties.VariableNames, 'spiketimes_'));
    
    for i = 1:numel(spike_vars)
        % Extract the numeric suffix (e.g. from 'spiketimes_3' get '3')
        suffix = extractAfter(spike_vars{i}, 'spiketimes_');
        
        % Construct corresponding names for spiketimes and netlabels
        spk_name = ['spiketimes_' suffix];
        net_name = ['netlabels_' suffix];
    
        % Only copy if both variables exist
        if ismember(net_name, tbl1.Properties.VariableNames)
            tbl.(spk_name) = tbl1.(spk_name);
            tbl.(net_name) = tbl1.(net_name);
        end
    end

    if ~isempty(LFP)
        tbl.lfp = LFP;
    end

end
