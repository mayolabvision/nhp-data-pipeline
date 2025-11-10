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
kilo_files1 = dir(fullfile(si_path, '*.mat'));
for i = 1:length(kilo_files1)
    kpath = fullfile(si_path, kilo_files1(i).name);
    [~, kname, ~] = fileparts(kpath);
    
    loaded = load(kpath);  % Load into struct
    varNames = fieldnames(loaded);  % Get the name(s) of loaded variables

    loaded_var = loaded.(varNames{1});
    if isa(loaded_var, 'single') || isa(loaded_var, 'int32') || isa(loaded_var, 'int64')
        loaded_var = double(loaded_var);
    end
    
    % Assign the first variable in the .mat file to the struct
    sorting.(kname) = loaded_var;
end

% .tsv files
kilo_files2 = dir(fullfile(si_path, '*.tsv'));
[~, sort_idx] = sort({kilo_files2.name}.');
kilo_files2 = kilo_files2(sort_idx);

for ii = 1:length(kilo_files2)
    cfpath = fullfile(kilo_files2(ii).folder, kilo_files2(ii).name);
    %[~, cfname, ~] = fileparts(kilo_files2(ii).name);
    if ii==1
        clusts = readtable(cfpath, 'FileType', 'text', 'Delimiter', '\t');
    else
        cc = readtable(cfpath, 'FileType', 'text', 'Delimiter', '\t');
        cc2 = join(clusts, cc, 'Keys', 'cluster_id');
        clusts = cc2;
    end
end
sorting.clusters = clusts; 

spike_times_sec = double(sorting.spike_times)./Fs;

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
        rp = tbl.ALIGN_PULSE{t,1};
        et =  tbl.END_TRIAL(t);
    
        spike_units = sorting.spike_clusters(spike_times_sec>=(np-(rp./1000)) & spike_times_sec<=(np+((et-rp)./1000)));
        spike_times = ((spike_times_sec((spike_times_sec>=(np-(rp./1000)) & spike_times_sec<=(np+((et-rp)./1000)))) - np)*1000) + rp;   
    end

    spikes_perTrial{t} = cellfun(@(u) spike_times(spike_units==u), num2cell(double(unique(sorting.spike_clusters)+1)), 'uni', 0);
end

spike_counts = cellfun(@(w) cellfun(@(q) numel(q), w, 'uni', 1), spikes_perTrial, 'uni', 0);
unit_frs = cellfun(@(w,v) w./v, spike_counts, num2cell((tbl.END_TRIAL - tbl.START_TRIAL)./1000), 'uni', 0);
unit_frs = unit_frs(cellfun(@sum, unit_frs, 'uni', 1)>0);
trlAvg_frs = mean(vertcat(unit_frs{:}),1,'omitnan')';

end
