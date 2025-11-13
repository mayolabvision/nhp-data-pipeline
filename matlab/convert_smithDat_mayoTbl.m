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
    tbl.block = tbl1.block;
    tbl.time_sec = tbl1.time;

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
    
    % Pull out conditions from trial names, separated by ';' delimeter
    if any(contains(TASK_NAME, {'rfmp', 'rfMapping'}))
        % Define the function to process each string
        process_string = @(input_str) [...
            str2double(cellfun(@(x) x{1}, regexp(input_str, 'xpos=([-0-9]+)', 'tokens'), 'UniformOutput', false))', ...
            str2double(cellfun(@(x) x{1}, regexp(input_str, 'ypos=([-0-9]+)', 'tokens'), 'UniformOutput', false))' ...
        ];

        % Apply the function to each cell in conditions
        conditions = cellfun(@(q) pix2deg(q, tbl1(1,:).params.block.screenDistance, tbl1(1,:).params.block.pixPerCM), cellfun(process_string, tbl1.text, 'uni', 0), 'uni', 0);
        tbl.conditions = cellfun(@(q) num2cell(q,2), conditions, 'uni', 0);

        tbl.STIM_ON(tbl.result~='CORRECT') = cellfun(@(q) q(1:end-1), tbl.STIM_ON(tbl.result~='CORRECT'), 'uni', 0);
        tbl.STIM_OFF(tbl.result~='CORRECT') = cellfun(@(q) q(1:end-1), tbl.STIM_OFF(tbl.result~='CORRECT'), 'uni', 0);
        tbl.conditions = cellfun(@(q,v) q(1:numel(v)), tbl.conditions, tbl.STIM_ON, 'uni', 0);

    else
        pattern = '([^0-9;]+)(?==)';
        matches = cellfun(@(q) regexp(q, pattern, 'match'), tbl1.text, 'uni', 0);

        % Get all unique condition names from all the matches
        cols = sort(unique(horzcat(matches{:})));
        conditions = cellfun(@(x) cellfun(@(q,r) str2double(x(q+1:r-1)), num2cell(strfind(x,'=')), num2cell(strfind(x,';')), 'uni', 0), tbl1.text, 'uni', 0);
        
        ordered_conditions = cell(size(conditions));
        for i = 1:length(conditions)
            temp = repmat({NaN}, 1, length(cols));
            
            current_row_cols = matches{i};
            for j = 1:length(current_row_cols)
                col_idx = find(strcmp(cols, current_row_cols{j}));
                
                temp{col_idx} = conditions{i}{j};
            end
            
            % Store the reordered conditions
            ordered_conditions{i} = temp;
        end
        
        % Convert the ordered conditions to a matrix
        conditions_matrix = cell2mat(vertcat(ordered_conditions{:}));
        
        % Create the table with the final ordered conditions
        for c = 1:length(cols)
            tbl.(cols{c}) = conditions_matrix(:, c);
        end
    end

    if any(contains(TASK_NAME, {'purs','pursuit'}))
        tbl1.result(tbl1.result==0 | tbl1.result==154) = 167;
    end

    tbl.params = tbl1.params;
    tbl.eyedata = tbl1.eyedata; tbl.pupil = tbl1.pupil; tbl.diode = tbl1.diode;

    % Add filtered eye traces and kinematic derivatives
    eyePos = cellfun(@(x) filterEyeTraces_EyeLink(x), tbl.eyedata, 'uni', 0);
    [eyeVel, eyeAcc] = cellfun(@(x) calcDerivative_eyeTraces(x), eyePos, 'uni', 0);
    tbl.eyePos = eyePos; tbl.eyeVel = eyeVel; tbl.eyeAcc = eyeAcc;

    %tbl.ns5_samps = tbl1.ns5_samps;

    if ~iscell(tbl.FIXATE)
        tbl.FIXATE = num2cell(tbl.FIXATE);
    end

    %%%%%%%%%%%%% task-specific re-arranging and calculations %%%%%%%%%%%%%
    if any(contains(TASK_NAME, {'mdir', 'dirmem'}))
        tbl(:, 'MEM_GUIDED_SACC') = [];

        if ~ismember('distance', tbl.Properties.VariableNames)
            tbl.distance = cellfun(@(q) q.distance, {tbl.params.block}.', 'uni', 1);
            tbl = movevars(tbl,{'distance'},'After','angle');
        end
        tbl.distance = cellfun(@(q) round(pix2deg(q,tbl(1,:).params.block.screenDistance,tbl(1,:).params.block.pixPerCM)), num2cell(tbl.distance), 'uni', 1);

        if all(ismember({'targetOnsetDelay', 'delay'}, tbl.Properties.VariableNames))
            tbl.fixDuration = tbl.targetOnsetDelay+tbl(1,:).params.block.targetDuration+tbl.delay;

            tbl = movevars(tbl,{'targetOnsetDelay','delay','fixDuration'},'Before','result');
        end

        tbl.saccades = cell(height(tbl), 1);
        validRows = cellfun(@(f) ~isempty(f) && all(~isnan(f)), tbl.FIXATE);
        tbl.saccades(validRows) = cellfun(@(v,f) ...
            num2cell(detect_saccades(v(:, f(1):end)) + f(1), 2), ...
            tbl.eyeVel(validRows), tbl.FIXATE(validRows), 'uni', 0);

        tbl.saccadeOnset = nan(height(tbl),1);
        for t = 1:height(tbl)
            x = cellfun(@(q) q(1), tbl.saccades{t}, 'uni', 1) - tbl.SACCADE(t);
            x(x>0)=NaN;
            [~,m] = min(abs(x));

            tbl.saccadeOnset(t) = tbl.saccades{t}{m}(1);
            tbl.saccadeOffset(t) = tbl.saccades{t}{m}(2);
        end
        
        tbl.saccadeLatency = tbl.saccadeOnset - cell2mat(tbl.FIX_OFF);

        % Det if first saccade out of fixation window landed in targ win
        inTargets = nan(height(tbl),1); 
        [dThetas,dRhos,dists] = deal(cell(height(tbl),1));
        for t = 1:height(tbl)
            % radial position of eye at s
            [theta_eye, rho_eye] = cart2pol(tbl.eyePos{t}(1,tbl.saccadeOffset(t)+50),tbl.eyePos{t}(2,tbl.saccadeOffset(t)+50));
            rho_targ = tbl.distance(t);
            theta_targ = deg2rad(tbl.angle(t));
            r_window = pix2deg(tbl.params(t).block.targWinRad,tbl.params(t).block.screenDistance,tbl.params(t).block.pixPerCM);

            theta_eye = mod(theta_eye, 2*pi);
            theta_targ = mod(theta_targ, 2*pi);

            % Signed difference: positive = clockwise
            dThetas{t} = rad2deg(- (mod(theta_eye - theta_targ + pi, 2*pi) - pi));
            dRhos{t} = rho_eye-rho_targ;

            % Compute distance using law of cosines
            dist = sqrt(rho_eye.^2 + rho_targ^2 - 2*rho_eye*rho_targ.*cos(theta_eye - theta_targ));
            dists{t} = dist;
            
            % Logical array: true if eye is inside target window
            inTargets(t) = dist <= r_window;  
        end

        tbl.saccadeOffset_dTheta = dThetas;
        tbl.saccadeOffset_dRho = dRhos;
        tbl.saccadeOffset_dist = dists;

    elseif any(contains(TASK_NAME, {'purs','pursuit'}))
        % Define the columns to replace and their new names
        cols_to_replace = {'TARG_ON', 'PURSUIT_TARG', 'TARG_OFF'};
        new_names = {'PURSUIT_TARG_ON', 'PURSUIT_TARG_ON', 'PURSUIT_TARG_OFF'};
        
        % Loop through each column to check and replace
        for i = 1:numel(cols_to_replace)
            if ismember(cols_to_replace{i}, tbl.Properties.VariableNames)
                tbl.Properties.VariableNames{ismember(tbl.Properties.VariableNames, cols_to_replace{i})} = new_names{i};
            end
        end

        if ~ismember('pursuitSpeed',tbl.Properties.VariableNames)
            tbl.pursuitSpeed = repmat(tbl(1,:).params.block.pursuitSpeed,height(tbl),1);
            tbl = movevars(tbl,{'pursuitSpeed'},'Before','fixDuration');
        end
        if ~ismember('jump',tbl.Properties.VariableNames)
            if isfield(tbl(1,:).params.block,'jump')
                tbl.jump = repmat(tbl(1,:).params.block.jump,height(tbl),1);
            else
                tbl.jump = zeros(height(tbl),1);
            end
            tbl = movevars(tbl,{'jump'},'Before','fixDuration');
        end

        tbl.saccades = cell(height(tbl), 1);
        validRows = cellfun(@(f) ~isempty(f) && all(~isnan(f)), tbl.FIXATE);
        tbl.saccades(validRows) = cellfun(@(v,f) ...
            num2cell(detect_saccades(v(:, f(1):end), 'VEL_THRESH', 30, 'ACC_THRESH', 500) + f(1), 2), ...
            tbl.eyeVel(validRows), tbl.FIXATE(validRows), 'uni', 0);

        [pursuitOnset, pursuitLatency] = cellfun(@(u,v,w) detect_pursuitOnset(u, v, w, 'PLOT_TRACES', false), tbl.eyeVel, num2cell(tbl.PURSUIT_TARG_ON), num2cell(tbl.pursuitSpeed), 'uni', 1); 
        tbl.pursuitOnset = pursuitOnset;
        tbl.pursuitLatency = pursuitLatency;

    elseif any(contains(TASK_NAME, {'rfmp', 'rfMapping'}))
        tbl.saccades = cell(height(tbl), 1);
        validRows = cellfun(@(f) ~isempty(f) && all(~isnan(f)), tbl.FIXATE);
        tbl.saccades(validRows) = cellfun(@(v,f) ...
            num2cell(detect_saccades(v(:, f(1):end)) + f(1), 2), ...
            tbl.eyeVel(validRows), tbl.FIXATE(validRows), 'uni', 0);
    end

    for v = ["emptyCnd", "STIM8_ON"]
        if ismember(v, tbl.Properties.VariableNames)
            tbl.(v) = [];
        end
    end

    if ismember('ALIGN_PULSE', tbl.Properties.VariableNames)
        tbl = movevars(tbl,{'ALIGN_PULSE'}, 'After','START_TRIAL');
    end

    if ismember('IGNORED', tbl.Properties.VariableNames) && ismember('CORRECT', tbl.Properties.VariableNames)
        tbl = movevars(tbl,{'IGNORED'},'After','CORRECT');
    end

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
