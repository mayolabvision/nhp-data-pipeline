function sorting = load_sortingOutputs(si_path)
%LOAD_SORTINGOUTPUTS  Load SpikeInterface/kilosort output files from si_path into a struct.
%   sorting = load_sortingOutputs(si_path)
%
%   Every .mat file in si_path is loaded as a top-level field of the
%   returned struct (e.g. sorting.spike_times, sorting.spike_clusters).
%   Every .tsv file (e.g. cluster_KSLabel.tsv, cluster_group.tsv) is
%   joined on 'cluster_id' into sorting.clusters.
%
% INPUT
%   si_path : path to the SpikeInterface/kilosort sorter_output folder
%
% OUTPUT
%   sorting : struct with one field per .mat file, plus sorting.clusters
%             (a table joining all .tsv files on cluster_id)

% .mat files
kilo_files1 = dir(fullfile(si_path, '*.mat'));
for i = 1:length(kilo_files1)
    kpath = fullfile(si_path, kilo_files1(i).name);
    [~, kname, ~] = fileparts(kpath);

    loaded = load(kpath); % Load into struct
    varNames = fieldnames(loaded); % Get the name(s) of loaded variables

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
    if ii == 1
        clusts = readtable(cfpath, 'FileType', 'text', 'Delimiter', '\t');
    else
        cc = readtable(cfpath, 'FileType', 'text', 'Delimiter', '\t');
        clusts = join(clusts, cc, 'Keys', 'cluster_id');
    end
end
sorting.clusters = clusts;

end
