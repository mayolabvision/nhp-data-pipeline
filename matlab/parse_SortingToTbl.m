function [spikes_perTrial, sorting, trlAvg_frs] = parse_SortingToTbl(tbl,si_path,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'si_path', @ischar);
addParameter(p, 'NP_ALIGN_PULSES', [], @isnumeric);
addParameter(p, 'Fs', 30000, @isnumeric);

parse(p, tbl, si_path, varargin{:});
tbl = p.Results.tbl;
si_path = p.Results.si_path;
NP_ALIGN_PULSES = p.Results.NP_ALIGN_PULSES;
Fs = p.Results.Fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in MATLAB table
tic
sorting = load_sortingOutputs(si_path);

spike_times_sec = double(sorting.spike_times)./Fs;
nClusters = max(sorting.spike_clusters) + 1; % cluster_id is 0-indexed (kilosort/phy); spikes_perTrial{t}{k} = cluster_id (k-1)

spikes_perTrial = cell(height(tbl),1);
for t = 1:height(tbl)
    if mod(t, 100) == 0
        disp(['Trial: ', num2str(t), '/', num2str(height(tbl))]);
    end

    if isempty(NP_ALIGN_PULSES)
        rp = tbl.time_sec(t,1);
        et = tbl.time_sec(t,2);

        spike_units = sorting.spike_clusters(spike_times_sec>=rp & spike_times_sec<=et);
        spike_times = ((spike_times_sec(spike_times_sec>=rp & spike_times_sec<=et)) - rp)*1000;
    else
        np = NP_ALIGN_PULSES(t);
        rp = tbl.ALIGN_PULSE{t,1}(1);
        et =  tbl.END_TRIAL(t);

        spike_units = sorting.spike_clusters(spike_times_sec>=(np-(rp./1000)) & spike_times_sec<=(np+((et-rp)./1000)));
        spike_times = ((spike_times_sec((spike_times_sec>=(np-(rp./1000)) & spike_times_sec<=(np+((et-rp)./1000)))) - np)*1000) + rp;   
    end

    spikes_perTrial{t} = cellfun(@(u) spike_times(spike_units==u), num2cell((0:nClusters-1)'), 'uni', 0);
end

spike_counts = cellfun(@(w) cellfun(@(q) numel(q), w, 'uni', 1), spikes_perTrial, 'uni', 0);
unit_frs = cellfun(@(w,v) w./v, spike_counts, num2cell((tbl.END_TRIAL - tbl.START_TRIAL)./1000), 'uni', 0);
unit_frs = unit_frs(cellfun(@sum, unit_frs, 'uni', 1)>0);
trlAvg_frs = mean(vertcat(unit_frs{:}),1,'omitnan')';

end
