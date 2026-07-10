function spike_times = extract_spikesInWindow(si_path, np_win_sec, varargin)
%EXTRACT_SPIKESINWINDOW  Pull spike times within one continuous window from SpikeInterface/kilosort output.
%
%   spike_times = extract_spikesInWindow(si_path, np_win_sec, 'Fs', 30000)
%
%   Unlike parse_SortingToTbl (which parses spikes into a per-trial cell
%   array using task alignment pulses/trial codes), this just pulls every
%   spike that occurred within one window [np_win_sec(1), np_win_sec(2)]
%   (seconds, on the neuropixel/imec clock) - no sync pulses or task codes
%   involved.
%
% INPUT
%   si_path    : path to the SpikeInterface/kilosort sorter_output folder
%   np_win_sec : [1 x 2] window [start, end], in seconds, on the same clock
%                sorting.spike_times is relative to (e.g. the imec clock)
%
% NAME-VALUE PARAMETERS
%   Fs  -  Sampling rate (Hz) that sorting.spike_times (in samples) is
%          relative to. Default 30000.
%
% OUTPUT
%   spike_times : [nClusters x 1] cell array. spike_times{k} is a vector of
%                 spike times (sec, relative to np_win_sec(1)) for
%                 cluster_id = k-1 (kilosort/phy cluster ids are 0-indexed;
%                 this is just that same id shifted by 1 for MATLAB).
%                 nClusters = max(cluster_id)+1 across the WHOLE recording
%                 (not just those with spikes in np_win_sec), so any
%                 cluster_id with no spikes in the window - or no spikes at
%                 all, e.g. a gap left by kilosort/phy curation - still
%                 gets its own entry at spike_times{cluster_id+1}, just an
%                 empty one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'si_path', @ischar);
addRequired(p, 'np_win_sec', @isnumeric);
addParameter(p, 'Fs', 30000, @isnumeric);

parse(p, si_path, np_win_sec, varargin{:});
si_path = p.Results.si_path;
np_win_sec = p.Results.np_win_sec;
Fs = p.Results.Fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sorting = load_sortingOutputs(si_path);

spike_times_sec = double(sorting.spike_times) ./ Fs;
spike_clusters = double(sorting.spike_clusters);

% Cluster count comes from the whole recording, not just the window, so
% clusters silent during np_win_sec still get an (empty) entry. Indexed
% directly by cluster_id+1 (not by rank), so gaps in cluster_id (e.g. from
% kilosort/phy curation) leave empty slots rather than shifting everything
% after them.
nClusters = max(spike_clusters) + 1;

in_window = spike_times_sec >= np_win_sec(1) & spike_times_sec <= np_win_sec(2);
win_times = spike_times_sec(in_window) - np_win_sec(1); % relative to window start
win_clusters = spike_clusters(in_window);

spike_times = cell(nClusters, 1);
for cluster_id = 0:(nClusters - 1)
    spike_times{cluster_id + 1} = win_times(win_clusters == cluster_id);
end

end
