% Define the directory containing the .mat files
dataDir = '/Users/kendranoneman/Data/pipeline_parameterization';

% Get list of all .mat files in the directory
files = dir(fullfile(dataDir, '*.mat'));

% Preallocate cell array: one row per file, 3 columns
allSessions = cell(numel(files), 3);

for i = 1:numel(files)
    fprintf('Processing file %d of %d\n', i, numel(files));

    % Full path to file
    datapath = fullfile(files(i).folder, files(i).name);
    
    % Load struct S from file
    load(datapath, 'S');
    fprintf('.......................\n')
    
    % Extract sess_name, protocol, clusters
    sess_name = S.sess_name;
    protocol  = S.protocol;
    
    % Extract clusters
    clusters1 = S.kilosort(1).clusters;
    clusters2 = S.kilosort(2).clusters;

    try
        clusters = [clusters1; clusters2];   
    catch ME
        fprintf('Skipping file %d of %d\n', i, numel(files));
        continue;
    end
    
    % Get union of variable names
    allVars = union(clusters1.Properties.VariableNames, clusters2.Properties.VariableNames);
    
    % Ensure both tables have all variables
    for v = allVars
        varName = v{1};
        if ~ismember(varName, clusters1.Properties.VariableNames)
            clusters1.(varName) = NaN(height(clusters1),1);
        end
        if ~ismember(varName, clusters2.Properties.VariableNames)
            clusters2.(varName) = NaN(height(clusters2),1);
        end
    end
    
    % ----- Add protocol info as columns -----
    
    % motion_correct flag
    if isfield(protocol.motion_correction, 'preprocessing')
        clusters.motion_correct = ones(height(clusters),1);
    else
        clusters.motion_correct = zeros(height(clusters),1);
    end
    
    % whitening_range
    if isfield(protocol, 'sorting') && isfield(protocol.sorting, 'whitening_range')
        wr = protocol.sorting.whitening_range;
        % Replicate the value for all rows of this table
        clusters.whitening_range = repmat(wr, height(clusters), 1);
    else
        clusters.whitening_range = NaN(height(clusters), 1);
    end
    
    % Store in cell array
    allSessions{i,1} = sess_name;
    allSessions{i,2} = protocol;
    allSessions{i,3} = clusters;
end

emptyRows = any(cellfun(@isempty, allSessions), 2);

% Remove those rows
allSessions(emptyRows, :) = [];

%%
% ----- Collapse all tables into one -----

% First, find all unique variable names across the tables
allVarNames = cellfun(@(q) q.Properties.VariableNames, allSessions(:,3), 'uni', 0);
allVarNames = unique(horzcat(allVarNames{:}));

% Standardize tables: ensure each has all columns
for i = 1:size(allSessions,1)
    tbl = allSessions{i,3};
    missingCols = setdiff(allVarNames, tbl.Properties.VariableNames);
    for m = 1:numel(missingCols)
        % Add missing column filled with NaN
        tbl.(missingCols{m}) = NaN(height(tbl),1);
    end
    % Reorder columns consistently
    tbl = tbl(:, allVarNames);
    allSessions{i,3} = tbl;
end

% Now concatenate all into one big table
collapsedTable = vertcat(allSessions{:,3});


%% Plots

sessions = unique(collapsedTable.sess_name);
yvars = {'snr';'firing_rate';'rp_violations'};

f = figure; %('Visible','off');
f.Position = [100 100 1200 900];
tl = tiledlayout(numel(yvars),2);

ax = gobjects(numel(yvars)+1, 2);

ax(1,1) = nexttile;
for s = 1:numel(sessions)
    T = collapsedTable(collapsedTable.motion_correct==1 & collapsedTable.sess_name==sessions(s),:);
    [G, groupVals] = findgroups(T.whitening_range);

    counts  = splitapply(@numel, T.cluster_id, G);

    bar(groupVals, counts, 'FaceAlpha', 0.5)
    hold on;

    ylabel('number of clusters')
    title('with drift correction')
    prettyFig;
end

ax(1,2) = nexttile;
for s = 1:numel(sessions)
    T = collapsedTable(collapsedTable.motion_correct==0 & collapsedTable.sess_name==sessions(s),:);
    [G, groupVals] = findgroups(T.whitening_range);

    counts  = splitapply(@numel, T.cluster_id, G);

    bar(groupVals, counts, 'FaceAlpha', 0.5)
    hold on;

    title('with NO drift correction')
    prettyFig;
end

linkaxes(ax(1,:), 'y');

for y = 1:numel(yvars)
    ax(y+1,1) = nexttile;
    for s = 1:numel(sessions)
        T = collapsedTable(collapsedTable.motion_correct==1 & collapsedTable.sess_name==sessions(s),:);
        [G, groupVals] = findgroups(T.whitening_range);
        
        meanSNR = splitapply(@mean, T.(yvars{y}), G);
        semSNR  = splitapply(@(x) std(x)/sqrt(numel(x)), T.(yvars{y}), G);
        
        errorbar(groupVals, meanSNR, semSNR, 'o-', 'LineWidth', 1.5)
        hold on;

        if y==numel(yvars)-1
            xlabel('whitening_range','interpreter','none');
        end
        ylabel(yvars{y},'Interpreter','none');
        prettyFig;
    end

    ax(y+1,2) = nexttile;
    for s = 1:numel(sessions)
        T = collapsedTable(collapsedTable.motion_correct==0 & collapsedTable.sess_name==sessions(s),:);
        [G, groupVals] = findgroups(T.whitening_range);
        
        meanSNR = splitapply(@mean, T.(yvars{y}), G);
        semSNR  = splitapply(@(x) std(x)/sqrt(numel(x)), T.(yvars{y}), G);
        
        errorbar(groupVals, meanSNR, semSNR, 'o-', 'LineWidth', 1.5)
        hold on;

        if y==numel(yvars)-1
            xlabel('whitening_range','interpreter','none');
            legend(sessions,'interpreter','none')
        end
        prettyFig;
    end

    linkaxes(ax(y+1,:), 'y');

end

linkaxes(ax(:), 'x');



































