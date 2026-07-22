function plot_sortingQC(data, figs_path)
%PLOT_SORTINGQC Generate spike-sorting QC summary figures for each probe.
%   PLOT_SORTINGQC(DATA, FIGS_PATH) produces, for every probe in
%   S.sorting, a set of population-level spike-sorting QC figures: drift
%   map, amplitude stability, ISI violations, spike counts, template
%   similarity, unit locations, autocorrelograms, mean waveforms, and
%   firing rate vs depth. Every plot summarizes all units in one figure
%   -- nothing is broken out per individual cluster. Figures are saved
%   as PNGs under a per-probe sub-directory of FIGS_PATH.
%
%   INPUTS
%     data      : either
%                   (1) char/string path to a processed session .mat
%                       file containing a struct (named S, or the sole
%                       variable in the file) with a S.sorting array
%                       (one element per probe), or
%                   (2) the already-loaded struct S itself
%     figs_path : root directory under which per-probe QC figure
%                 sub-directories are created
%
%   Examples
%     plot_sortingQC('/path/to/session-hash.mat', '/path/to/figs')
%     plot_sortingQC(S, '/path/to/figs')

if nargin < 2
    error('plot_sortingQC:missingInputs', 'Both data and figs_path are required.');
end
if ~exist(figs_path, 'dir')
    mkdir(figs_path);
end

if ischar(data) || isstring(data)
    S = local_loadStruct(char(data));
elseif isstruct(data)
    S = data;
else
    error('plot_sortingQC:invalidInput', 'data must be a char/string path or a struct.');
end

if ~isfield(S, 'sorting')
    error('plot_sortingQC:noSorting', 'Loaded struct has no ''sorting'' field.');
end

if isfield(S, 'sess_name') && ~isempty(S.sess_name)
    sess_name = char(S.sess_name);
elseif ischar(data) || isstring(data)
    [~, sess_name] = fileparts(data);
else
    sess_name = 'session';
end

nProbes = numel(S.sorting);
for pIdx = 1:nProbes
    sorting = S.sorting(pIdx);

    % ----- locate this probe's existing figs sub-directory (e.g. imec0_lFEF) -----
    % probe_index is 1-based (imec0 <-> probe_index==1), matching SpikeGLX imec numbering.
    if isfield(sorting, 'probe_index') && ~isempty(sorting.probe_index)
        imec_tag = sprintf('imec%d', sorting.probe_index - 1);
    else
        imec_tag = sprintf('imec%d', pIdx - 1);
    end

    probe_label = '';
    if isfield(sorting, 'clusters') && istable(sorting.clusters) && height(sorting.clusters) > 0 && ...
            ismember('probe_label', sorting.clusters.Properties.VariableNames)
        lbl = sorting.clusters.probe_label(1);
        if iscell(lbl); lbl = lbl{1}; end
        if iscategorical(lbl) || isstring(lbl); lbl = char(lbl); end
        probe_label = lbl;
    end
    probe_tag = imec_tag;
    if ~isempty(probe_label)
        probe_tag = sprintf('%s_%s', imec_tag, probe_label);
    end

    existing = dir(fullfile(figs_path, [imec_tag '_*']));
    existing = existing([existing.isdir]);
    if numel(existing) >= 1
        if numel(existing) > 1
            warning('plot_sortingQC:ambiguousProbeDir', ...
                'Multiple directories match ''%s_*'' in %s; using ''%s''.', ...
                imec_tag, figs_path, existing(1).name);
        end
        probe_dir = fullfile(existing(1).folder, existing(1).name);
    else
        probe_dir = fullfile(figs_path, probe_tag);
        if ~exist(probe_dir, 'dir')
            mkdir(probe_dir);
        end
    end

    save_dir = fullfile(probe_dir, 'sortingQC');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    % ----- pull + condition the raw fields for this probe -----
    % NOTE: spike_times is stored in samples (see parse_SortingToTbl.m),
    % so it is converted to seconds here using the AP sample rate.
    Fs = local_sampleRate(sorting, 30000);
    spike_times       = double(sorting.spike_times(:))' ./ Fs;
    spike_clusters    = double(sorting.spike_clusters(:))';
    templates         = double(sorting.templates);
    amplitudes        = double(sorting.amplitudes(:))';
    spike_positions   = double(sorting.spike_positions);
    channel_positions = double(sorting.channel_positions);

    fprintf('[plot_sortingQC] %s / %s: %d spikes, %d units, Fs=%.1f Hz\n', ...
        sess_name, probe_tag, numel(spike_times), numel(unique(spike_clusters)), Fs);

    % ----- per-unit sorting label (good/mua/noise), ascending unit id -----
    unit_ids = unique(spike_clusters);
    cluster_labels = repmat({'unlabeled'}, numel(unit_ids), 1);
    if isfield(sorting, 'clusters') && istable(sorting.clusters) && ...
            all(ismember({'cluster_id', 'KSLabel_cc'}, sorting.clusters.Properties.VariableNames))
        ct = sorting.clusters;
        for u = 1:numel(unit_ids)
            row = find(ct.cluster_id == unit_ids(u), 1);
            if ~isempty(row)
                lbl = ct.KSLabel_cc(row);
                if iscell(lbl); lbl = lbl{1}; end
                cluster_labels{u} = char(lbl);
            end
        end
    end

    % ----- build + save each QC figure -----
    save_qcFig(plot_driftMap(spike_times, spike_positions, amplitudes), save_dir, '01_drift_map');
    save_qcFig(plot_ampOverTime(spike_times, spike_clusters, amplitudes), save_dir, '02_amplitude_over_time');
    save_qcFig(plot_isiViolationHist(spike_times, spike_clusters, cluster_labels), save_dir, '03_isi_violations');
    save_qcFig(plot_spikeCountHist(spike_clusters), save_dir, '04_spike_count_hist');
    save_qcFig(plot_templateCorrHeatmap(templates, spike_clusters, spike_positions), save_dir, '05_template_correlation');
    save_qcFig(plot_unitPositionScatter(spike_clusters, spike_positions, channel_positions), save_dir, '06_unit_positions');
    save_qcFig(plot_acgGallery(spike_times, spike_clusters), save_dir, '07_acg_gallery');
    save_qcFig(plot_waveformGallery(templates, spike_clusters, spike_positions), save_dir, '08_waveform_gallery');
    save_qcFig(plot_frVsDepth(spike_times, spike_clusters, spike_positions), save_dir, '09_firing_rate_vs_depth');

    fprintf('[plot_sortingQC] saved figures to %s\n', save_dir);
end

end


%% ======================= plot sub-functions ============================

function fig = plot_driftMap(spike_times, spike_positions, amplitudes)
%PLOT_DRIFTMAP Spike depth vs time, colored by amplitude, all units pooled.
fig = new_darkFigure([1000 600]);
ax = axes(fig);
logAmp = log10(max(amplitudes, eps));
scatter(ax, spike_times/60, spike_positions(:,2), 4, logAmp, 'filled', ...
    'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0);
colormap(ax, 'parula');
cb = colorbar(ax); cb.Color = 'w'; cb.Label.String = 'log_{10} amplitude'; cb.Label.Color = 'w';
xlabel(ax, 'time (min)'); ylabel(ax, 'depth ({\mu}m)');
title(ax, 'Drift map');
style_darkAxes(ax);
end


function fig = plot_ampOverTime(spike_times, spike_clusters, amplitudes)
%PLOT_AMPOVERTIME Per-unit amplitude trend over time (median + 10-90th pct band), all units overlaid.
nBins = 30;
edges = linspace(min(spike_times), max(spike_times), nBins+1);
binCenters = (edges(1:end-1) + edges(2:end)) / 2 / 60; % minutes

unit_ids = unique(spike_clusters);
nUnits = numel(unit_ids);

medAmpByUnit = nan(nUnits,1);
for i = 1:nUnits
    medAmpByUnit(i) = median(amplitudes(spike_clusters == unit_ids(i)), 'omitnan');
end
[~, rankOrder] = sort(medAmpByUnit);
unitRank = nan(nUnits,1);
unitRank(rankOrder) = 1:nUnits;
cmap = turbo(max(nUnits,2));

fig = new_darkFigure([1100 650]);
ax = axes(fig); hold(ax, 'on');
for i = 1:nUnits
    mask = spike_clusters == unit_ids(i);
    ts = spike_times(mask);
    amps = amplitudes(mask);
    if numel(ts) < 5
        continue
    end
    binIdx = discretize(ts, edges);
    med = nan(1, nBins); p10 = nan(1, nBins); p90 = nan(1, nBins);
    for b = 1:nBins
        v = amps(binIdx == b);
        if isempty(v)
            continue
        end
        med(b) = median(v);
        p10(b) = local_prctile(v, 10);
        p90(b) = local_prctile(v, 90);
    end
    valid = ~isnan(med);
    if nnz(valid) < 2
        continue
    end
    c = cmap(unitRank(i), :);
    fill(ax, [binCenters(valid), fliplr(binCenters(valid))], [p10(valid), fliplr(p90(valid))], c, ...
        'FaceAlpha', 0.05, 'EdgeColor', 'none');
    plot(ax, binCenters(valid), med(valid), 'Color', [c 0.5], 'LineWidth', 0.75);
end
xlabel(ax, 'time (min)'); ylabel(ax, 'amplitude');
title(ax, sprintf('Amplitude stability (%d units, color = median amplitude rank)', nUnits));
style_darkAxes(ax);
end


function fig = plot_isiViolationHist(spike_times, spike_clusters, cluster_labels)
%PLOT_ISIVIOLATIONHIST Histogram of per-unit ISI violation rate, split by sorting label.
%   cluster_labels : cellstr, one label per unique(spike_clusters), ascending unit-id order
refractory_s = 0.0015; % 1.5 ms

unit_ids = unique(spike_clusters);
nUnits = numel(unit_ids);
viol_rate = zeros(nUnits, 1);
for i = 1:nUnits
    ts = sort(spike_times(spike_clusters == unit_ids(i)));
    if numel(ts) < 2
        continue
    end
    isis = diff(ts);
    viol_rate(i) = mean(isis < refractory_s);
end

if nargin < 3 || isempty(cluster_labels)
    cluster_labels = repmat({'unlabeled'}, nUnits, 1);
end
labels = categorical(cluster_labels(:));
catList = categories(labels);

fig = new_darkFigure([900 550]);
ax = axes(fig); hold(ax, 'on');
edges = linspace(0, max(0.02, max(viol_rate)), 40);
cmap = lines(numel(catList));
for c = 1:numel(catList)
    mask = labels == catList{c};
    histogram(ax, viol_rate(mask), edges, 'FaceColor', cmap(c,:), 'FaceAlpha', 0.65, ...
        'EdgeColor', 'none', 'DisplayName', sprintf('%s (n=%d)', catList{c}, nnz(mask)));
end
legend(ax, 'TextColor', 'w', 'Color', 'none', 'EdgeColor', [0.5 0.5 0.5], 'Location', 'northeast');
xlabel(ax, 'ISI violation rate'); ylabel(ax, '# units');
title(ax, 'ISI violations by sorting label');
style_darkAxes(ax);
end


function fig = plot_spikeCountHist(spike_clusters)
%PLOT_SPIKECOUNTHIST Log-scale histogram of spike count per unit.
[G, ~] = findgroups(spike_clusters(:));
counts = splitapply(@numel, spike_clusters(:), G);

fig = new_darkFigure([800 550]);
ax = axes(fig);
edges = logspace(log10(max(1, min(counts))), log10(max(counts)), 30);
histogram(ax, counts, edges, 'FaceColor', [0.3 0.7 0.9], 'EdgeColor', 'none');
set(ax, 'XScale', 'log');
xlabel(ax, 'spike count (log scale)'); ylabel(ax, '# units');
title(ax, sprintf('Spike count distribution (%d units)', numel(counts)));
style_darkAxes(ax);
end


function fig = plot_templateCorrHeatmap(templates, spike_clusters, spike_positions)
%PLOT_TEMPLATECORRHEATMAP Pairwise template correlation, units ordered by depth.
[unit_ids, ~] = local_depthSortedUnits(templates, spike_clusters, spike_positions);
nUnits = numel(unit_ids);

flatWave = zeros(nUnits, size(templates,2) * size(templates,3));
for i = 1:nUnits
    tmpl = squeeze(templates(unit_ids(i) + 1, :, :));
    flatWave(i, :) = tmpl(:)';
end
C = corrcoef(flatWave');

fig = new_darkFigure([850 800]);
ax = axes(fig);
imagesc(ax, C, [-1 1]);
axis(ax, 'square');
colormap(ax, 'turbo');
cb = colorbar(ax); cb.Color = 'w'; cb.Label.String = 'correlation'; cb.Label.Color = 'w';
xlabel(ax, 'unit (sorted by depth)'); ylabel(ax, 'unit (sorted by depth)');
title(ax, sprintf('Pairwise template correlation (%d units)', nUnits));
style_darkAxes(ax);
end


function fig = plot_unitPositionScatter(spike_clusters, spike_positions, channel_positions)
%PLOT_UNITPOSITIONSCATTER Unit x vs y location, marker size/color = spike count.
[G, gid] = findgroups(spike_clusters(:));
count = splitapply(@numel, spike_clusters(:), G);
x = splitapply(@(v) median(v, 'omitnan'), spike_positions(:,1), G);
y = splitapply(@(v) median(v, 'omitnan'), spike_positions(:,2), G);

countRange = max(count) - min(count);
if countRange == 0
    sz = 40 * ones(size(count));
else
    sz = 10 + 90 * (count - min(count)) / countRange;
end

fig = new_darkFigure([650 900]);
ax = axes(fig); hold(ax, 'on');
scatter(ax, x, y, sz, log10(count), 'filled', 'MarkerFaceAlpha', 0.8);
colormap(ax, 'parula');
cb = colorbar(ax); cb.Color = 'w'; cb.Label.String = 'log_{10} spike count'; cb.Label.Color = 'w';
if nargin >= 3 && ~isempty(channel_positions)
    xlim(ax, [min(channel_positions(:,1))-20, max(channel_positions(:,1))+20]);
    ylim(ax, [min(channel_positions(:,2))-20, max(channel_positions(:,2))+20]);
end
xlabel(ax, 'x ({\mu}m)'); ylabel(ax, 'y ({\mu}m, depth)');
title(ax, sprintf('Unit positions (%d units)', numel(gid)));
style_darkAxes(ax);
end


function fig = plot_acgGallery(spike_times, spike_clusters)
%PLOT_ACGGALLERY Autocorrelogram gallery for the top 30 units by spike count.
[G, gid] = findgroups(spike_clusters(:));
counts = splitapply(@numel, spike_clusters(:), G);
[~, order] = sort(counts, 'descend');
topN = min(30, numel(gid));
topUnits = gid(order(1:topN));
topCounts = counts(order(1:topN));

binSize = 0.001; % 1 ms
maxLag = 0.05;   % 50 ms
edges = -maxLag:binSize:maxLag;
centers_ms = (edges(1:end-1) + edges(2:end)) / 2 * 1000;

nCols = 6;
nRows = ceil(topN / nCols);
fig = new_darkFigure([220*nCols, 160*nRows]);
tl = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:topN
    ts = sort(spike_times(spike_clusters == topUnits(i)));
    acg = local_acg(ts, edges);
    ax = nexttile(tl);
    bar(ax, centers_ms, acg, 1, 'FaceColor', [0.4 0.75 0.95], 'EdgeColor', 'none');
    hold(ax, 'on'); xline(ax, 0, 'Color', [0.6 0.6 0.6]);
    title(ax, sprintf('u%d (n=%d)', topUnits(i), topCounts(i)), 'FontSize', 8);
    set(ax, 'XTick', [], 'YTick', []);
    style_darkAxes(ax);
end
title(tl, 'Autocorrelograms (top 30 units by spike count)', 'Color', 'w');
end


function fig = plot_waveformGallery(templates, spike_clusters, spike_positions)
%PLOT_WAVEFORMGALLERY Best-channel mean waveform gallery, units ordered by depth.
[unit_ids, depth] = local_depthSortedUnits(templates, spike_clusters, spike_positions);
nUnits = numel(unit_ids);

nCols = 8;
nRows = ceil(nUnits / nCols);
fig = new_darkFigure([160*nCols, 110*nRows]);
tl = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nUnits
    tmpl = squeeze(templates(unit_ids(i) + 1, :, :));
    [~, bestCh] = max(range(tmpl, 1));
    wf = tmpl(:, bestCh);
    ax = nexttile(tl);
    plot(ax, wf, 'Color', [0.4 0.85 0.6], 'LineWidth', 1);
    title(ax, sprintf('u%d (%d{\\mu}m)', unit_ids(i), round(depth(i))), 'FontSize', 7);
    axis(ax, 'off');
end
title(tl, 'Mean waveforms, best channel (sorted shallow \rightarrow deep)', 'Color', 'w');
end


function fig = plot_frVsDepth(spike_times, spike_clusters, spike_positions)
%PLOT_FRVSDEPTH Mean firing rate vs depth, one point per unit.
[G, gid] = findgroups(spike_clusters(:));
count = splitapply(@numel, spike_clusters(:), G);
depth = splitapply(@(v) median(v, 'omitnan'), spike_positions(:,2), G);
duration = max(spike_times) - min(spike_times);
fr = count / duration;

fig = new_darkFigure([650 850]);
ax = axes(fig);
scatter(ax, fr, depth, 20, [0.95 0.6 0.3], 'filled', 'MarkerFaceAlpha', 0.7);
set(ax, 'XScale', 'log');
xlabel(ax, 'mean firing rate (Hz, log scale)'); ylabel(ax, 'depth ({\mu}m)');
title(ax, sprintf('Firing rate vs depth (%d units)', numel(gid)));
style_darkAxes(ax);
end


%% ============================ helpers ===================================

function S = local_loadStruct(data_path)
loaded = load(data_path);
if isfield(loaded, 'S')
    S = loaded.S;
    return
end
fn = fieldnames(loaded);
if numel(fn) == 1
    S = loaded.(fn{1});
    return
end
error('plot_sortingQC:noStructFound', ...
    'Could not find a single struct variable (or one named ''S'') in %s', data_path);
end


function Fs = local_sampleRate(sorting, defaultFs)
Fs = defaultFs;
if ~isfield(sorting, 'ap_meta') || isempty(sorting.ap_meta)
    return
end
meta = sorting.ap_meta;
candidateFields = {'imSampRate', 'sample_rate', 'sRateHz', 'fs', 'Fs'};
if isstruct(meta)
    for k = 1:numel(candidateFields)
        if isfield(meta, candidateFields{k}) && isnumeric(meta.(candidateFields{k})) && isscalar(meta.(candidateFields{k}))
            Fs = double(meta.(candidateFields{k}));
            return
        end
    end
elseif istable(meta)
    for k = 1:numel(candidateFields)
        if ismember(candidateFields{k}, meta.Properties.VariableNames)
            val = meta.(candidateFields{k});
            if isnumeric(val) && isscalar(val)
                Fs = double(val);
                return
            end
        end
    end
end
end


function [unit_ids_sorted, depth_sorted] = local_depthSortedUnits(templates, spike_clusters, spike_positions)
%LOCAL_DEPTHSORTEDUNITS Unit ids present in both spike_clusters and templates, ordered by depth.
unit_ids = unique(spike_clusters);
nTemplateRows = size(templates, 1);
valid = (unit_ids >= 0) & (unit_ids + 1 <= nTemplateRows);
if any(~valid)
    warning('plot_sortingQC:idMismatch', ...
        '%d unit(s) have cluster IDs outside the templates array and were skipped.', nnz(~valid));
end
unit_ids = unit_ids(valid);

[G, gid] = findgroups(spike_clusters(:));
depthAll = splitapply(@(v) median(v, 'omitnan'), spike_positions(:,2), G);
depth = depthAll(ismember(gid, unit_ids));

[depth_sorted, order] = sort(depth);
unit_ids_sorted = unit_ids(order);
end


function acgCounts = local_acg(ts, edges)
%LOCAL_ACG Symmetric autocorrelogram counts for a sorted spike-time vector.
maxLag = edges(end);
n = numel(ts);
acgCounts = zeros(1, numel(edges) - 1);
if n < 2
    return
end
dtAll = cell(n, 1);
for i = 1:n
    j = i + 1;
    dtLocal = [];
    while j <= n && (ts(j) - ts(i)) <= maxLag
        dtLocal(end+1) = ts(j) - ts(i); %#ok<AGROW>
        j = j + 1;
    end
    dtAll{i} = dtLocal;
end
dt = [dtAll{:}];
acgCounts = histcounts([dt, -dt], edges);
end


function y = local_prctile(v, p)
%LOCAL_PRCTILE Linear-interpolation percentile, avoids Statistics Toolbox dependency.
v = sort(v(~isnan(v)));
n = numel(v);
if n == 0
    y = NaN;
    return
end
if n == 1
    y = v(1);
    return
end
idx = (p/100) * (n-1) + 1;
lo = floor(idx); hi = ceil(idx);
if lo == hi
    y = v(lo);
else
    y = v(lo) + (v(hi) - v(lo)) * (idx - lo);
end
end


function fig = new_darkFigure(sz)
fig = figure('Color', [0.08 0.08 0.08], 'Position', [100 100 sz(1) sz(2)], 'Visible', 'off');
end


function style_darkAxes(ax)
set(ax, 'Color', [0.12 0.12 0.12], 'XColor', [0.85 0.85 0.85], 'YColor', [0.85 0.85 0.85], ...
    'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.3, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
if ~isempty(ax.Title); ax.Title.Color = 'w'; end
end


function save_qcFig(fig, save_dir, name)
outFile = fullfile(save_dir, [name '.png']);
exportgraphics(fig, outFile, 'Resolution', 200, 'BackgroundColor', fig.Color);
close(fig);
end
