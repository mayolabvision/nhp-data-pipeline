function metrics = parse_clusterMetrics(kilosort4_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'kilosort4_path', @ischar);
%addParameter(p, 'NP_ALIGN_PULSES', [], @isnumeric);
%addParameter(p, 'Fs', 30000, @isnumeric);

parse(p, kilosort4_path); %, varargin{:});
ks_path = p.Results.kilosort4_path;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

metrics = readtable(fullfile(kilosort4_path,'quality_metrics','cluster_metrics.csv'), 'Delimiter', ',');

metrics.sess_name = categorical(metrics.sess_name);
metrics.monkey    = categorical(metrics.monkey);
metrics.experimenter = categorical(metrics.experimenter);
metrics.probe_label = categorical(metrics.probe_label);
metrics.probe_type = categorical(metrics.probe_type);
metrics.probe_config = categorical(metrics.probe_config);
metrics.hardware_config = categorical(metrics.hardware_config);
% Only run if the column exists
if ismember('probe_gridHole', metrics.Properties.VariableNames)
    % Convert each row
    metrics.probe_gridHole = cellfun(@(x) ...
        strsplit(regexprep(x, {'[\[\]'' ]'}, ''), ','), ...
        metrics.probe_gridHole, 'UniformOutput', false);

    % Convert to column cell arrays instead of row
    metrics.probe_gridHole = cellfun(@(x) x(:), metrics.probe_gridHole, ...
        'UniformOutput', false);
end

% Convert any table column whose cells are char strings of numpy-like arrays
vars = metrics.Properties.VariableNames;

for k = 1:numel(vars)
    v = vars{k};
    col = metrics.(v);

    if iscell(col) && ~isempty(col)
        isCharCell = cellfun(@(x) ischar(x) || isstring(x), col);
        if all(isCharCell | cellfun('isempty', col))
            metrics.(v) = cellfun(@numpyStr2Mat, col, 'UniformOutput', false);
        end
    end
end

function A = numpyStr2Mat(s)
    % Turn numpy-like string arrays into MATLAB numeric arrays.
    % Handles 1-D, 2-D, and 3-D+ array string representations.
    % If the string contains "..." (numpy summarization), returns [].

    if isempty(s); A = []; return; end
    if isstring(s); s = char(s); end
    s = strtrim(s);

    % Strip outer quotes if present
    if ~isempty(s) && (s(1) == '''' || s(1) == '"') && s(end) == s(1)
        s = s(2:end-1);
    end

    % If truncated numpy repr ("..."), just bail
    if contains(s, '...')
        A = [];
        return
    end

    % Normalize whitespace/newlines and commas
    s = regexprep(s, '\r?\n', ' ');
    s = regexprep(s, '\s+', ' ');
    s = strrep(s, ',', ' ');

    % Core trick:
    % Replace inner closing+opening brackets with semicolons at *all* depths.
    % e.g., "][", "] [", "]] [[" → use semicolons to split rows/planes.
    while contains(s, '][')
        s = regexprep(s, '\]\s*\[', '; ');
    end

    % Remove any remaining outer brackets
    s = regexprep(s, '^\[|\]$', '');

    % Wrap back into MATLAB matrix syntax
    s = ['[', s, ']'];

    A = str2num(s); %#ok<ST2NM>
    
end

end
