function [tbl,taskNickname] = handle_taskSpecifics(in_tbl)
% handle_taskSpecifics
%
% USAGE:
%   task_keywords = handle_taskSpecifics();
%   [task_keywords, newTbl, taskNickname] = handle_taskSpecifics(tbl);
%
% Splits tbl by its unique xmlName values (since a single table can mix
% trials from multiple tasks, e.g. interleaved multi-task blocks), runs
% the task-specific processing on each xmlName's rows separately, then
% recombines the pieces. Columns that only apply to some tasks are
% NaN/empty-padded for the rows where they don't apply, and the combined
% table is re-sorted (via trialName) back to tbl's original row order.
%
% taskNickname is the hyphen-joined short nickname for whichever tasks
% are present in tbl.xmlName (e.g. a table mixing dirmem and pursuit
% trials returns 'mdir-purs'); see the task_nicknames dictionary below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dictionary mapping each task's short nickname to the exact xmlName(s)
% that identify it. This drives both the per-task branch dispatch below
% and the taskNickname output.
task_nicknames = { ...
    'rfmp', {'rfMapping_dots'}; ...
    'mdir', {'dirmem_withhelp', 'dirmem_withhelp_varDelays'}; ...
    'purs', {'pursuit_task'}; ...
    'fstm', {'fixAndStim'}; ...
    'sarf', {'dirmemAndRFmap'}; ...
    'purf', {'pursuitAndRFmap'}; ...
    'cfix', {'onlyFixationAndCircularVis'}; ...
    'oflu', {'onlyFixationAndLuminance'}; ...
    'absp', {'ANIblocked_pursuitsacc'}; ...
};

%===========================================================================================

tbl = in_tbl;

% FIX_OFF can come in as a plain numeric column or a cell column
% (depending on whether any trial in the whole session has multiple
% FIX_OFF codes); normalize it once, up front, so every xmlName group
% agrees on its type before we split and recombine below.
if ismember('FIX_OFF', tbl.Properties.VariableNames) && ~iscell(tbl.FIX_OFF)
    tbl.FIX_OFF = num2cell(tbl.FIX_OFF);
end

xmlNames = unique(tbl.xmlName);
tbl_parts = cell(numel(xmlNames), 1);
nicknames = cell(numel(xmlNames), 1);
for k = 1:numel(xmlNames)
    nicknames{k} = lookupNickname(char(xmlNames(k)), task_nicknames);
    tbl_parts{k} = processTaskGroup(tbl(tbl.xmlName == xmlNames(k), :), nicknames{k});
end

tbl = combineTaskGroups(tbl_parts);

% Splitting/recombining by xmlName scrambles row order - restore tbl's
% original order using trialName, which is unique per row.
[~, order] = ismember(in_tbl.trialName, tbl.trialName);
tbl = tbl(order, :);

taskNickname = strjoin(unique(nicknames, 'stable'), '-');

end

function nickname = lookupNickname(taskName, task_nicknames)
%LOOKUPNICKNAME  Look up taskName's short nickname in the task_nicknames dictionary.
nickname = 'unknown';
for r = 1:size(task_nicknames, 1)
    if ismember(taskName, task_nicknames{r,2})
        nickname = task_nicknames{r,1};
        return;
    end
end
end

function tbl = processTaskGroup(tbl, nickname)
%PROCESSTASKGROUP  Apply task-specific processing to one task's rows,
%   dispatched by nickname (see task_nicknames dictionary above).

switch nickname

%---------------------------------------- rfMapping_dots ----------------------------------------%
case 'rfmp'
    process_string = @(input_str) [...
        str2double(cellfun(@(x) x{1}, regexp(input_str, 'xpos=([-0-9]+)', 'tokens'), 'UniformOutput', false))', ...
        str2double(cellfun(@(x) x{1}, regexp(input_str, 'ypos=([-0-9]+)', 'tokens'), 'UniformOutput', false))' ...
    ];

    % Apply the function to each cell in conditions
    conditions = cellfun(@(q) pix2deg(q, tbl(1,:).params.block.screenDistance, tbl(1,:).params.block.pixPerCM), cellfun(process_string, tbl.text, 'uni', 0), 'uni', 0);
    tbl.conditions = cellfun(@(q) num2cell(q,2), conditions, 'uni', 0);

    tbl.STIM_ON(tbl.result~='CORRECT') = cellfun(@(q) q(1:end-1), tbl.STIM_ON(tbl.result~='CORRECT'), 'uni', 0);
    tbl.STIM_OFF(tbl.result~='CORRECT') = cellfun(@(q) q(1:end-1), tbl.STIM_OFF(tbl.result~='CORRECT'), 'uni', 0);
    tbl.conditions = cellfun(@(q,v) q(1:numel(v)), tbl.conditions, tbl.STIM_ON, 'uni', 0);

    tbl.saccades = cell(height(tbl), 1);
    validRows = cellfun(@(f) ~isempty(f) && all(~isnan(f)), tbl.FIXATE);
    tbl.saccades(validRows) = cellfun(@(v,f) ...
        num2cell(detect_saccades(v(:, f(1):end)) + f(1), 2), ...
        tbl.eyeVel(validRows), tbl.FIXATE(validRows), 'uni', 0);

%---------------------------------------- dirmem_withhelp ----------------------------------------%
case 'mdir'
    tbl = expandConditionText(tbl);

    tbl(:, 'MEM_GUIDED_SACC') = [];

    if ~ismember('angle', tbl.Properties.VariableNames)
        tbl.angle = cellfun(@(q) q.distance, {tbl.params.block}.', 'uni', 1);
        tbl.angle(tbl.result=='NaN') = NaN;
    end

    if ~ismember('distance', tbl.Properties.VariableNames)
        tbl.distance = cellfun(@(q) q.distance, {tbl.params.block}.', 'uni', 1);
        tbl.distance(tbl.result=='NaN') = NaN;
    end
    tbl = movevars(tbl,{'distance','angle'},'After','time_sec');
    tbl.distance = cellfun(@(q) round(pix2deg(q,tbl(1,:).params.block.screenDistance,tbl(1,:).params.block.pixPerCM)), num2cell(tbl.distance), 'uni', 1);

    if all(ismember({'targetOnsetDelay', 'delay'}, tbl.Properties.VariableNames))
        tbl.fixDuration = tbl.targetOnsetDelay+tbl(1,:).params.block.targetDuration+tbl.delay;
        tbl.fixDuration(tbl.result=='NaN') = NaN;

        tbl = movevars(tbl,{'targetOnsetDelay','delay','fixDuration'},'Before','result');
    end

    tbl.saccades = cell(height(tbl), 1);
    validRows = cellfun(@(f) ~isempty(f) && all(~isnan(f)), tbl.FIXATE);
    tbl.saccades(validRows) = cellfun(@(v,f) ...
        num2cell(detect_saccades(v(:, f(1):end)) + f(1), 2), ...
        tbl.eyeVel(validRows), tbl.FIXATE(validRows), 'uni', 0);

    [~,rVel] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), tbl.eyeVel, 'uni', 0);

    tbl.saccadeOnset = nan(height(tbl),1); tbl.saccadeOffset = nan(height(tbl),1); tbl.saccadeLatency = nan(height(tbl),1);
    for t = 1:height(tbl)
        if isequal(tbl.result(t),'CORRECT')
            x = cellfun(@(q) q(1), tbl.saccades{t}, 'uni', 1) - tbl.SACCADE(t);
            x(x>0)=NaN;
            [~,m] = min(abs(x));

            tbl.saccadeOnset(t) = tbl.saccades{t}{m}(1);
            tbl.saccadeOffset(t) = tbl.saccades{t}{m}(2);

            eyeSac = rVel{t}(tbl.saccadeOnset(t):tbl.saccadeOnset(t)+200);
            [~,idx] = max(eyeSac);
            eyeSac(idx:end) = NaN;
            saccadeOnset2 = tbl.saccadeOnset(t)+(find(diff(eyeSac)<0,1,'last')+1);
            if ~isempty(saccadeOnset2)
                tbl.saccadeOnset(t) = tbl.saccadeOnset(t)+(find(diff(eyeSac)<0,1,'last')+1);
            end

            tbl.saccadeLatency(t) = tbl.saccadeOnset(t) - tbl.FIX_OFF{t};
        end
    end

    % Det if first saccade out of fixation window landed in targ win
    inTargets = nan(height(tbl),1);
    [dThetas,dRhos,dists] = deal(cell(height(tbl),1));
    for t = 1:height(tbl)
        if isequal(tbl.result(t),'CORRECT')
            % radial position of eye at s
            if (tbl.saccadeOffset(t) + 50) > size(tbl.eyePos{t}, 2)
                continue
            end

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
    end

    tbl.saccadeOffset_dTheta = dThetas;
    tbl.saccadeOffset_dRho = dRhos;
    tbl.saccadeOffset_dist = dists;

%---------------------------------------- pursuit_task ----------------------------------------%
case 'purs'
    tbl = expandConditionText(tbl);

    % Define the columns to replace and their new names
    cols_to_replace = {'TARG_ON', 'PURSUIT_TARG', 'TARG_OFF', 'pursuitSpeed'};
    new_names = {'PURSUIT_TARG_ON', 'PURSUIT_TARG_ON', 'PURSUIT_TARG_OFF', 'speed'};

    % Loop through each column to check and replace
    for i = 1:numel(cols_to_replace)
        if ismember(cols_to_replace{i}, tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{ismember(tbl.Properties.VariableNames, cols_to_replace{i})} = new_names{i};
        end
    end

    if ~ismember('speed',tbl.Properties.VariableNames)
        tbl.speed = repmat(tbl(1,:).params.block.pursuitSpeed,height(tbl),1);
    end
    if ~ismember('angle',tbl.Properties.VariableNames)
        tbl.angle = repmat(tbl(1,:).params.block.angle,height(tbl),1);
    end
    if ~ismember('jump',tbl.Properties.VariableNames)
        if isfield(tbl(1,:).params.block,'jump')
            tbl.jump = repmat(tbl(1,:).params.block.jump,height(tbl),1);
        else
            tbl.jump = zeros(height(tbl),1);
        end
    end
    tbl = movevars(tbl,{'angle','speed','jump'},'After','time_sec');

    tbl.saccades = cell(height(tbl), 1);
    validRows = cellfun(@(f) ~isempty(f) && all(~isnan(f)), tbl.FIXATE);
    tbl.saccades(validRows) = cellfun(@(v,f) ...
        num2cell(detect_saccades(v(:, f(1):end), 'VEL_THRESH', 30, 'ACC_THRESH', 500) + f(1), 2), ...
        tbl.eyeVel(validRows), tbl.FIXATE(validRows), 'uni', 0);

    [pursuitOnset, pursuitLatency] = cellfun(@(u,v,w) detect_pursuitOnset(u, v, w, 'PLOT_TRACES', false), tbl.eyeVel, num2cell(tbl.PURSUIT_TARG_ON), num2cell(tbl.speed), 'uni', 1);
    tbl.pursuitOnset = pursuitOnset;
    tbl.pursuitLatency = pursuitLatency;

%---------------------------------------- Ani Block Saccade-Pursuit ----------------------------------------%
case 'absp'
    tbl = expandConditionText(tbl);

    % Define the columns to replace and their new names
    cols_to_replace = {'thisTrialAngle', 'pursuitSpeed'};
    new_names = {'angle', 'speed'};

    % Loop through each column to check and replace
    for i = 1:numel(cols_to_replace)
        if ismember(cols_to_replace{i}, tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{ismember(tbl.Properties.VariableNames, cols_to_replace{i})} = new_names{i};
        end
    end

    colsToAdd = ["distance","speed","jump","crossingTime","saccadeOnset","saccadeOffset","saccadeLatency","saccadeOffset_dTheta", "saccadeOffset_dRho", "saccadeOffset_dist", "pursuitOnset", "pursuitLatency"];
    for c = colsToAdd
        if ~ismember(c, tbl.Properties.VariableNames)
            tbl.(c) = nan(height(tbl), 1);
        end
    end
    tbl.saccades = cell(height(tbl), 1);

    [~,rVel] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), tbl.eyeVel, 'uni', 0);

    % LOOP THROUGH EACH ROW OF TABLE
    for row = 1:height(tbl)
        if ~isnan(tbl.FIXATE{row})
            tbl.saccades(row) = {num2cell(detect_saccades(tbl.eyeVel{row}(:, tbl.FIXATE{row}(1):end)), 2)};
        end

        if ~isnan(tbl.SACCADE_BLOCK(row)) % SACCADE BLOCK
            tbl.pursuit_fixDuration(row) = NaN;
            if isnan(tbl.distance(row))
                tbl.distance(row) = round(pix2deg(tbl.params(row).block.distance, tbl.params(row).block.screenDistance, tbl.params(row).block.pixPerCM));
            end

            if isequal(tbl.result(row),'CORRECT')
                x = cellfun(@(q) q(1), tbl.saccades{row}, 'uni', 1) - tbl.SACCADE(row);
                x(x>0)=NaN;
                [~,m] = min(abs(x));

                tbl.saccadeOnset(row) = tbl.saccades{row}{m}(1);
                tbl.saccadeOffset(row) = tbl.saccades{row}{m}(2);

                eyeSac = rVel{row}(tbl.saccadeOnset(row):tbl.saccadeOnset(row)+200);
                [~,idx] = max(eyeSac);
                eyeSac(idx:end) = NaN;
                saccadeOnset2 = tbl.saccadeOnset(row)+(find(diff(eyeSac)<0,1,'last')+1);
                if ~isempty(saccadeOnset2)
                    tbl.saccadeOnset(row) = tbl.saccadeOnset(row)+(find(diff(eyeSac)<0,1,'last')+1);
                end

                tbl.saccadeLatency(row) = tbl.saccadeOnset(row) - tbl.FIX_OFF{row};

                [theta_eye, rho_eye] = cart2pol(tbl.eyePos{row}(1,tbl.saccadeOffset(row)+50),tbl.eyePos{row}(2,tbl.saccadeOffset(row)+50));
                rho_targ = tbl.distance(row);
                theta_targ = deg2rad(tbl.angle(row));
                %r_window = pix2deg(tbl.params(row).block.targWinRad,tbl.params(row).block.screenDistance,tbl.params(row).block.pixPerCM);

                theta_eye = mod(theta_eye, 2*pi);
                theta_targ = mod(theta_targ, 2*pi);

                % Signed difference: positive = clockwise
                tbl.saccadeOffset_dTheta(row) = rad2deg(- (mod(theta_eye - theta_targ + pi, 2*pi) - pi));
                tbl.saccadeOffset_dRho(row) = rho_eye-rho_targ;

                % Compute distance using law of cosines
                dist = sqrt(rho_eye.^2 + rho_targ^2 - 2*rho_eye*rho_targ.*cos(theta_eye - theta_targ));
                tbl.saccadeOffset_dist(row) = dist;

            end

        else % PURSUIT BLOCK
            tbl.delay(row) = NaN;
            tbl.sacc_fixDuration(row) = NaN;
            tbl.targetOnsetDelay(row) = NaN;

            if isnan(tbl.jump(row))
                tbl.jump(row) = tbl.params(row).block.jump;
            end
            if isnan(tbl.speed(row))
                tbl.speed(row) = tbl.params(row).block.pursuitSpeed;
            end
            if isnan(tbl.crossingTime(row))
                tbl.crossingTime(row) = tbl.params(row).block.crossingTime;
            end

            [pursuitOnset, pursuitLatency] = detect_pursuitOnset(tbl.eyeVel{row}, tbl.PURSUIT_TARG_ON(row), tbl.speed(row), 'PLOT_TRACES', false);
            tbl.pursuitOnset(row) = pursuitOnset;
            tbl.pursuitLatency(row) = pursuitLatency;
        end
    end
    tbl = movevars(tbl,{'angle','distance','speed','jump','crossingTime'},'After','time_sec');
    tbl = movevars(tbl,{'params','eyedata','pupil','diode','eyePos','eyeVel','eyeAcc'},'After','pursuitLatency');

otherwise % all other tasks, not defined here
    tbl = expandConditionText(tbl);
end

end

function tbl = combineTaskGroups(tbl_parts)
%COMBINETASKGROUPS  Union columns across per-task table pieces and vertcat them.
%   Columns that only exist in some pieces are padded (NaN for numeric,
%   <undefined> for categorical, empty for cell) in the pieces missing
%   them. If a column is cell-typed in at least one piece but not
%   another (e.g. a task-specific field represented as a cell array in
%   one branch and a plain numeric column in another), the numeric
%   version is wrapped into a cell so all pieces agree before vertcat.

allCols = {};
for k = 1:numel(tbl_parts)
    allCols = union(allCols, tbl_parts{k}.Properties.VariableNames, 'stable');
end

for c = allCols
    col = c{1};
    hasCol = cellfun(@(t) ismember(col, t.Properties.VariableNames), tbl_parts);
    isCellSomewhere = any(cellfun(@(t) iscell(t.(col)), tbl_parts(hasCol)));

    if ~isCellSomewhere
        firstPart = tbl_parts{find(hasCol, 1)};
        sample = firstPart.(col);
    end

    for k = 1:numel(tbl_parts)
        h = height(tbl_parts{k});
        if hasCol(k)
            if isCellSomewhere && ~iscell(tbl_parts{k}.(col))
                tbl_parts{k}.(col) = num2cell(tbl_parts{k}.(col), 2);
            end
        elseif isCellSomewhere
            tbl_parts{k}.(col) = cell(h, 1);
        elseif iscategorical(sample)
            tbl_parts{k}.(col) = repmat(categorical(missing), h, 1);
        else
            tbl_parts{k}.(col) = nan(h, 1);
        end
    end
end

for k = 1:numel(tbl_parts)
    tbl_parts{k} = tbl_parts{k}(:, allCols);
end

tbl = vertcat(tbl_parts{:});

end

function tbl = expandConditionText(tbl)
    pattern = '([^0-9;]+)(?==)';
    matches = cellfun(@(q) regexp(q, pattern, 'match'), tbl.text, 'uni', 0);

    % Get all unique condition names
    cols = sort(unique(horzcat(matches{:})));

    % Extract numeric values for each condition
    conditions = cellfun(@(x) ...
        cellfun(@(q,r) str2double(x(q+1:r-1)), ...
        num2cell(strfind(x,'=')), ...
        num2cell(strfind(x,';')), ...
        'uni', 0), ...
        tbl.text, 'uni', 0);

    % Reorder conditions to match column order
    ordered_conditions = cell(size(conditions));
    for i = 1:numel(conditions)
        temp = repmat({NaN}, 1, numel(cols));

        current_row_cols = matches{i};
        for j = 1:numel(current_row_cols)
            col_idx = strcmp(cols, current_row_cols{j});
            temp{col_idx} = conditions{i}{j};
        end

        ordered_conditions{i} = temp;
    end

    % Convert to matrix and add to table
    conditions_matrix = cell2mat(vertcat(ordered_conditions{:}));
    for c = 1:numel(cols)
        tbl.(cols{c}) = conditions_matrix(:, c);
        tbl = movevars(tbl,{cols{c}},'After','result');
    end
end
