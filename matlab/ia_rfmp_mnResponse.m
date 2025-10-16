function ia_rfmp_mnResponse(data,varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    p = inputParser;
    addRequired(p, 'data',  @(x) (ischar(x)) || isstruct(x));
    addParameter(p, 'FIG_PATH', [], @ischar);
    addParameter(p, 'PROBE_INDEX', 1, @isnumeric);

    parse(p, data, varargin{:});
    data = p.Results.data;
    FIG_PATH = p.Results.FIG_PATH;
    PROBE_INDEX = p.Results.PROBE_INDEX;

    fprintf('\n------------------------------\n')
    if ischar(data)
        [~, filename, ~] = fileparts(data);
        load(data,'S');
        fprintf(sprintf('\n----Data loaded for %s----\n',filename))
    else
        filename = data.sess_name;
        S = data;
    end

    probe_label = string(S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.probe_label(1));
    hardware_config = string(S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.hardware_config(1));

    % Find rfmp or rfMapping fields
    fields = fieldnames(S);
    matchingFields = fields(contains(fields, {'rfmp', 'rfMapping'}, 'IgnoreCase', true));

    matchingFields = matchingFields(cellfun(@(f) isfield(S.(f), 'tbl') && ~isempty(S.(f).tbl), matchingFields));

    task_conditions = cell(length(matchingFields), 5);
    for mm = 1:numel(matchingFields)
        stimuli = cell2mat(vertcat(S.(matchingFields{mm}).tbl.conditions{:}));
        task_conditions{mm, 1} = length(unique(stimuli(:,1))) * length(unique(stimuli(:,2))); 
        task_conditions{mm, 2} = S.(matchingFields{mm}).params.frameCount(1); 
        task_conditions{mm, 3} = S.(matchingFields{mm}).params.nStimPerFix(1);
        task_conditions{mm, 4} = S.(matchingFields{mm}).params.bgColor(1,:);
        task_conditions{mm, 5} = [S.(matchingFields{mm}).params.colorR S.(matchingFields{mm}).params.colorG S.(matchingFields{mm}).params.colorB];
    end

    row_strings = cell(size(task_conditions,1),1);
    for i = 1:size(task_conditions,1)
        row_strings{i} = strjoin(cellfun(@mat2str, task_conditions(i,:), 'UniformOutput', false), '|');
    end
    [~, unique_idx, row_conds] = unique(row_strings, 'stable');
    unique_conditions = task_conditions(unique_idx, :);

    if numel(unique(row_conds))==1
        T = []; 
        for mm = 1:numel(matchingFields)
            T = [T; S.(matchingFields{mm}).tbl];
        end

        T.STIM_ON(T.result~='CORRECT') = cellfun(@(q) q(1:end-1), T.STIM_ON(T.result~='CORRECT'), 'uni', 0);
        T.conditions(T.result~="CORRECT") = cellfun(@(q) q(1:end-1), T.conditions(T.result~='CORRECT'), 'uni', 0);
        T = T(~cellfun(@(q) any(isnan(q)), T.STIM_OFF, 'uni', 1),:);
        T = T(T.result=='CORRECT',:);
        T = T(logical(cellfun(@(q) numel(q)==11, T.conditions, 'uni', 1)),:);

        prb_name = ['spiketimes_' num2str(PROBE_INDEX)];

        stim_duration_ms = T.params(1,1).block.frameCount*(1000/120);
        probe_rgb = [T.params(1,1).block.colorR, T.params(1,1).block.colorG, T.params(1,1).block.colorB];
        bg_rgb = T.params(1,1).block.bgColor;

        bg_text = '';
        if isequal(bg_rgb, [0 0 0])
            bg_text = 'black bg';
        elseif isequal(bg_rgb, [255 255 255])
            bg_text = 'white bg';
        end

        if ~isempty(FIG_PATH)
            FIG_PATH2 = fullfile(FIG_PATH, sprintf('%s_%s',hardware_config,probe_label), 'rfmp_mnResponses', matchingFields{1});
            if ~exist("FIG_PATH2", 'dir'), mkdir(FIG_PATH2); end    
        end

        bin_size = 5;
        scale = 1000 / bin_size;
        binned_frs = cellfun(@(w,v) cellfun(@(q) (bin_spktimes(q, 'StartTime', v(1), 'EndTime', v(1)+1000, 'BinSize', bin_size))*scale, w, 'uni', 0), ...
            T.(prb_name), T.STIM_ON, 'uni', 0);

        binned_frs = vertcat(binned_frs{:});
        stacked_frs = arrayfun(@(j) vertcat(binned_frs{:,j}), 1:size(binned_frs,2), 'UniformOutput', false);
        mean_frs = cellfun(@(q) mean(q,1), stacked_frs, 'uni', 0); 

        chunkSize = 25; 
        for startIdx = 1:chunkSize:numel(binned_frs)
            endIdx = min(startIdx+chunkSize-1, numel(binned_frs));
            curCells = mean_frs(startIdx:endIdx);

            f = figure; %('Visible','off');
            f.Position = [100 100 1800 900];
            tl = tiledlayout(5,5);

            for k = 1:numel(curCells)
                nexttile;
                data = curCells{k};   % 212x40 double
                x = linspace(1,1000,numel(data));
                plot(x,data,'linewidth',1,'Color',[0 0 1]);   % example: plot the mean across rows
                title(sprintf('%d', startIdx+k-1));
                prettyFig;
            end

            xlabel(tl, 'time aligned to first stimulus onset of trial (ms)','fontsize',16);
            ylabel(tl, 'mean FR across trials (Hz)','fontsize',16);

            blah = 1;
        end

    else
        for mm = 1:size(unique_conditions,1)
            these_rows = find(row_conds==mm);
            T = []; 
            for r = 1:numel(these_rows)
                T = [T; S.(matchingFields{these_rows(r)}).tbl];
            end
            
            T.STIM_ON(T.result~='CORRECT') = cellfun(@(q) q(1:end-1), T.STIM_ON(T.result~='CORRECT'), 'uni', 0);
            T.conditions(T.result~="CORRECT") = cellfun(@(q) q(1:end-1), T.conditions(T.result~='CORRECT'), 'uni', 0);
            T = T(~cellfun(@(q) any(isnan(q)), T.STIM_OFF, 'uni', 1),:);
            T = T(T.result=='CORRECT',:);
            T = T(logical(cellfun(@(q) numel(q)==11, T.conditions, 'uni', 1)),:);

            prb_name = ['spiketimes_' num2str(PROBE_INDEX)];

            stim_duration_ms = T.params(1,1).block.frameCount*(1000/120);
            probe_rgb = [T.params(1,1).block.colorR, T.params(1,1).block.colorG, T.params(1,1).block.colorB];
            bg_rgb = T.params(1,1).block.bgColor;

            bg_text = '';
            if isequal(bg_rgb, [0 0 0])
                bg_text = 'black bg';
            elseif isequal(bg_rgb, [255 255 255])
                bg_text = 'white bg';
            end

            if ~isempty(FIG_PATH)
                FIG_PATH2 = fullfile(FIG_PATH, sprintf('%s_%s',hardware_config,probe_label), 'rfmp_mnResponses', strjoin(matchingFields(these_rows), '_'));
                if ~exist("FIG_PATH2", 'dir'), mkdir(FIG_PATH2); end    
            end

            bin_size = 1;
            scale = 1000 / bin_size;
            binned_frs = cellfun(@(w,v) cellfun(@(q) (bin_spktimes(q, 'StartTime', v(1), 'EndTime', v(1)+1000, 'BinSize', bin_size)')*scale, w, 'uni', 0), ...
                T.(prb_name), T.STIM_ON, 'uni', 0);

            %binned_frs = cellfun(@(w) cellfun(@(q) bin_spktimes(q,0,round((T.params(1).block.frameCount / (120/1000)) * T.params(1).block.nStimPerFix),T.params(1).block.frameCount/(120/1000)), w, 'uni', 0), T.(prb_name), 'uni', 0);
            binned_frs = vertcat(binned_frs{:});
            binned_frs = arrayfun(@(j) vertcat(binned_frs{:,j}), 1:size(binned_frs,2), 'uni', 0);

            chunkSize = 25; 
            for startIdx = 1:chunkSize:numel(binned_frs)
                endIdx = min(startIdx+chunkSize-1, numel(binned_frs));
                curCells = binned_frs(startIdx:endIdx);

                f = figure; %('Visible','off');
                f.Position = [100 100 1800 900];
                tl = tiledlayout(5,5);

                for k = 1:numel(curCells)
                    nexttile;
                    %x = linspace(1,sum(T.END_TRIAL),numel(data));
                    data = curCells{k};   % 212x40 double
                    plot(data);   % example: plot the mean across rows
                    title(sprintf('%d', startIdx+k-1));
                    prettyFig;
                end

                blah = 1;
            end

            
        end

    end

end