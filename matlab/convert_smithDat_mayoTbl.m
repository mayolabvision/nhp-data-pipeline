function tbl = convert_smithDat_mayoTbl(dat,dat_iti,varargin)
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
    addRequired(p, 'dat', @isstruct);
    addRequired(p, 'dat_iti', @isstruct);
    addParameter(p, 'LFP', []);
    addParameter(p, 'HELPERS_PATH', defaultHELP_PATH, @ischar);
    addParameter(p, 'INCLUDE_ITI', true, @islogical);

    % Parse the inputs
    parse(p, dat, dat_iti, varargin{:});

    % Assign parsed values to variables
    dat = p.Results.dat;
    dat_iti = p.Results.dat_iti;
    LFP = p.Results.LFP;
    HELPERS_PATH = p.Results.HELPERS_PATH;
    INCLUDE_ITI = p.Results.INCLUDE_ITI;

    addpath(fullfile(HELPERS_PATH,'behavior'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tbl1 = struct2table(dat); tbl2 = struct2table(dat_iti);
    tbl = table(); 

    if INCLUDE_ITI
        % trial names, trl = within-trial and iti = inter-trial
        tbl.trialName = [cellfun(@(q) ['trl.', sprintf('%04d', q)], num2cell(1:height(tbl1))', 'uni', 0); cellfun(@(q) ['iti.', sprintf('%04d', q)], num2cell(1:height(tbl2))', 'uni', 0)];
        tbl.trialName = categorical(string(tbl.trialName));

        tbl1 = [tbl1; tbl2];
    else
        % trial names, trl = within-trial and iti = inter-trial
        tbl.trialName = cellfun(@(q) ['trl.', sprintf('%04d', q)], num2cell(1:height(tbl1))', 'uni', 0);
        tbl.trialName = categorical(string(tbl.trialName));

        tbl1 = tbl1;
    end

    tbl1 = reformatMultiTaskBlock(tbl1);
    tbl.xmlName = categorical(erase(arrayfun(@(p) p.block.xmlFile, tbl1.params, 'uni', 0), '.xml'));

    if ismember('block', tbl1.Properties.VariableNames), tbl.block = tbl1.block; end
    if ismember('time', tbl1.Properties.VariableNames), tbl.time_sec = tbl1.time; end

    if iscell(tbl1.result)
        tbl.result = convertBetween_eventCodes_eventNames(tbl1.result);
    else
        tbl.result = convertBetween_eventCodes_eventNames(num2cell(tbl1.result));
    end
    tbl.result = categorical(string(tbl.result));
    tbl.result(tbl.result=='ALIGN_PULSE') = categorical("NaN");

    % Make array of times of start/end time per trial, for aligning with trial codes and indexing eye data
    times_ms = cellfun(@(q) round(q(1)*1000:(q(2)+1)*1000), num2cell(tbl1.time,2), 'uni', 0);
    trialStarts = cellfun(@(q,r) find(q == round(r(r(:,2)==1,3)*1000)), times_ms, tbl1.trialcodes, 'uni', 0);
    allCodes = tbl1.trialcodes; allCodes = vertcat(allCodes{:});
    allCodes = sort(unique(allCodes(:,2)));

    if isfield(tbl1.params(1,1).block, 'posShiftForCode')
        posShiftForCode = tbl1.params(1,1).block.posShiftForCode;   % UPDATE to 50000 once you switch the ex files over
        maxDotOffsetPix = 5000;    % generous bound - well beyond any realistic screen half-width
        isPosCode = allCodes > (posShiftForCode - maxDotOffsetPix) & ...
                    allCodes < (posShiftForCode + maxDotOffsetPix);
        
        stimPositions = cellfun(@(q) extractStimPositions(q, posShiftForCode), tbl1.trialcodes, 'uni', 0);
        if any(cellfun(@(q) ~isempty(q), stimPositions))
            tbl.stimPos = stimPositions;
        end
        eventCodes = num2cell(allCodes(~isPosCode)');
    else
        eventCodes = num2cell(allCodes');
    end

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
    tbl.text(tbl.result=='NaN') = {''};

    tbl.params = tbl1.params;
    tbl.eyedata = tbl1.eyedata; tbl.pupil = tbl1.pupil; tbl.diode = tbl1.diode;

    if ismember('ALIGN_PULSE', tbl.Properties.VariableNames)
        if INCLUDE_ITI
            names = string(tbl.trialName);
            itiRows = find(contains(names, '.iti.'));
            for i = 1:numel(itiRows)
                r = itiRows(i);
                trlName = strrep(names(r), '.iti.', '.trl.');
                tbl.ALIGN_PULSE(r) = cellfun(@(q) q - (tbl.END_TRIAL(names==trlName)+1), tbl.ALIGN_PULSE(names==trlName), 'uni', 0);
            end
        else
             tbl.ALIGN_PULSE = num2cell(cell2mat(tbl.ALIGN_PULSE), 2);
        end
    end

    [~, idx] = sort(tbl.time_sec(:,1));
    tbl = tbl(idx,:);

    % try
    %     [~, tbl] = handle_taskSpecifics(tbl, TASK_NAME);
    % catch
    %     disp('------------- task-specific additions failed -------------');
    % end

    %tbl.ns5_samps = tbl1.ns5_samps;

    % Add spiking data if nasnet was used 
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

function xyPix = extractStimPositions(trialcodes, posShiftForCode)
%EXTRACTDOTPOSITIONS Recover RF-map distractor (x,y) positions, in order, from a trial's codes.
%   xyPix = extractDotPositions(trialcodes, posShiftForCode) takes a
%   trial's code matrix (e.g. dat(n).trialcodes, or allCodes{n}.codes),
%   with columns [~, code, time], finds every STIM_ON (code 10) event,
%   and reads the two codes immediately following it as that flash's
%   shifted (x,y) position. Returns an N x 2 matrix - one row per
%   STIM_ON, in the order it occurred - with the shift removed so values
%   are actual pixel offsets from screen center.
%
%   posShiftForCode must match whatever shift was used when the codes
%   were generated (10000 for early sessions, 50000 going forward) -
%   pass it explicitly, there's no way to auto-detect which was used.

    STIM_ON = 10;

    codeCol = trialcodes(:,2);
    stimOnIdx = find(codeCol == STIM_ON);

    xyPix = nan(numel(stimOnIdx), 2);

    for k = 1:numel(stimOnIdx)
        i = stimOnIdx(k);

        if i+2 > numel(codeCol)
            warning('extractDotPositions:truncated', ...
                'STIM_ON at row %d has no room for two position codes after it - skipping.', i);
            continue;
        end

        xCode = codeCol(i+1);
        yCode = codeCol(i+2);

        % sanity check - a real position code should sit near the shift,
        % not equal another known event code (e.g. 40 = STIM_OFF), which
        % would mean the trial got cut short right after STIM_ON
        if xCode == 40 || yCode == 40
            warning('extractDotPositions:unexpectedCode', ...
                'STIM_ON at row %d is not followed by two position codes - skipping.', i);
            continue;
        end

        xyPix(k,:) = [xCode - posShiftForCode, yCode - posShiftForCode];
    end

end

function tbl1 = reformatMultiTaskBlock(tbl1)
%REFORMATMULTITASKBLOCK Collapse multi-task tbl1.params.block cell arrays
%down to just the task that actually ran on each trial.
%
%   tbl1 = reformatMultiTaskBlock(tbl1) loops over every row of tbl1. For
%   rows where tbl1.params(row,1).block is a cell array (one struct per
%   interleaved task), it parses 'taskNum=N' out of tbl1.text{row} and
%   replaces tbl1.params(row,1).block with tbl1.params(row,1).block{N} -
%   i.e., just the block struct for whichever task actually ran on that
%   trial. Rows where .block is already a plain struct (single-task
%   sessions) are left untouched.

    for i = 1:height(tbl1)

        blockVal = tbl1.params(i,1).block;

        if ~iscell(blockVal)
            continue
        end

        taskNumTok = regexp(tbl1.text{i}, 'taskNum=(\d+)', 'tokens', 'once');
        if isempty(taskNumTok)
            warning('reformatMultiTaskBlock:noTaskNum', ...
                'Row %d: block is a cell array but no taskNum found in tbl1.text - leaving unchanged.', i);
            continue
        end
        taskNum = str2double(taskNumTok{1});

        if taskNum < 1 || taskNum > numel(blockVal)
            warning('reformatMultiTaskBlock:badTaskNum', ...
                'Row %d: taskNum=%d is out of range for a %d-element block cell array - leaving unchanged.', ...
                i, taskNum, numel(blockVal));
            continue
        end

        tbl1.params(i,1).block = blockVal{taskNum};

    end

end
