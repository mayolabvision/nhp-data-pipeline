function [eyedata, eye_times, trial_events] = extract_eyeStream(nev, out_ns5, varargin)
%EXTRACT_EYESTREAM  Build one continuous 1 kHz eye-position stream for a recording.
%
%   [eyedata, eye_times, trial_events] = extract_eyeStream(nev, out_ns5, 'EYE_CHAN_LABELS', eye_chan_labels)
%
%   Unlike format_datTrials (which segments eye data into a struct array,
%   one entry per trial/ITI), this pulls one continuous block of eye data
%   spanning the whole recording - starting 1 second before the first sync
%   pulse (first digital code == 0) and ending at the last recorded sample
%   in out_ns5 - and downsamples it once from 30 kHz to 1 kHz.
%
%%%% Required inputs: %%%
%   nev      -  nev event data (same format accepted by format_datTrials:
%               either a numeric [channel, code, time] array, or a struct
%               with .nev/.net_labels fields when spike-sorted)
%   out_ns5  -  ns5 struct containing the raw 30 kHz data (out_ns5.data)
%               and header (out_ns5.hdr)
%
%%%% Optional parameters: %%%
%   EYE_CHAN_LABELS -  Cell array specifying the eye movement channels.
%                       Default {'10241', '10242'} ({Eye_HE, Eye_VE}).
%   PREV_TEMPDATA   -  Carry-over calibration params from a previous file,
%                       used only if this recording has no pre-trial setup
%                       text to parse (mirrors format_datTrials). Default [].
%
%%%% Outputs: %%%
%   eyedata      -  [N x 2] continuous eye position (degrees), at 1 kHz
%   eye_times    -  [N x 1] elapsed time (sec) of each row of eyedata;
%                   always starts at 0 = 1 second before the first sync
%                   pulse (or the start of that sync pulse's recording
%                   epoch, if the epoch started less than 1 second earlier),
%                   i.e. [0; 0.001; 0.002; ...]
%   trial_events -  [N x 2] one row per eye_times sample:
%                     column 1 = trial number (NaN during inter-trial intervals)
%                     column 2 = digital event code that occurred in that
%                                1 ms bin (NaN if no code occurred). Only
%                                codes < 256 or in [1000, 32000] are kept
%                                (same event-code filter format_datTrials
%                                uses to build trialcodes); codes 256-511
%                                are text/annotation payload, not events.
%
% CALIBRATION NOTE: eye2deg needs a calibration params struct. Because this
% function does one bulk extraction instead of looping trial-by-trial, it
% only uses the FIRST block's calibration (parsed from the pre-trial setup
% text) for the entire stream - it does not re-parse calibration if a block
% transition happens partway through the recording, unlike format_datTrials.
%
% TRIAL NUMBERING: sync pulses (code == 0) are sent twice per trial, ~100 ms
% apart. So trial 1 starts at the 1st occurrence of code 0, trial 2 starts
% at the 3rd occurrence, trial k starts at the (2k-1)-th occurrence. Each
% trial's labeled region runs from that start time through its own
% end-trial (255) code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultEyeChanLabels = {'10241', '10242'};

p = inputParser;
addRequired(p, 'nev', @(x) (isnumeric(x)) || isstruct(x));
addRequired(p, 'out_ns5', @isstruct);
addParameter(p, 'EYE_CHAN_LABELS', defaultEyeChanLabels, @(x) iscell(x));
addParameter(p, 'PREV_TEMPDATA', [], @isstruct);

parse(p, nev, out_ns5, varargin{:});

nev = p.Results.nev;
out_ns5 = p.Results.out_ns5;
eye_channel_labels = p.Results.EYE_CHAN_LABELS;
prev_tempdata = p.Results.PREV_TEMPDATA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs = double(out_ns5.hdr.Fs); % Sampling frequency
endtrial = 255; % Trial end code
syncCode = 0; % Sync pulse code

eye_channels = find(ismember(out_ns5.hdr.label, eye_channel_labels));

% Determine if nev is an array of struct; if so, extract nev data
if isequal(class(nev), 'struct')
    NEV = [nev.nev nev.net_labels'];
else
    NEV = nev;
end

digcodes = NEV(NEV(:, 1) == 0, :); % All digital codes across the whole recording

% Continuous window: starts 1 second before the first sync pulse (but never
% before the start of the recording epoch that sync pulse falls in - can't
% extract data that wasn't recorded), and ends at the very last recorded
% sample in out_ns5
syncTimes = digcodes(digcodes(:, 2) == syncCode, 3);
if isempty(syncTimes)
    error('extract_eyeStream:noSyncPulse', 'No sync pulse (code %d) found in this recording.', syncCode);
end
sync_time = syncTimes(1);

sync_native = round(sync_time * Fs);
containing_epoch = find(sync_native >= out_ns5.hdr.timeStamps(1,:) & sync_native <= out_ns5.hdr.timeStamps(2,:), 1, 'first');
if isempty(containing_epoch)
    error('extract_eyeStream:syncOutsideEpochs', 'First sync pulse does not fall within any recorded epoch.');
end
epoch_start_time = out_ns5.hdr.timeStamps(1, containing_epoch) / Fs;
t_start = max(sync_time - 1, epoch_start_time);
t_end = out_ns5.hdr.timeStamps(2, end) / Fs; % last recorded time bin in out_ns5

endTimes = digcodes(digcodes(:, 2) == endtrial, 3);
if isempty(endTimes)
    error('extract_eyeStream:noEndTrial', 'No end-trial code (%d) found in this recording.', endtrial);
end

% Map [t_start, t_end] to sample indices into out_ns5.data (accounting for
% recording epochs, since out_ns5.data only stores samples for recorded
% epochs and skips any gap between them)
samp_start = time_to_ns5_sample(t_start, out_ns5, Fs);
samp_end = size(out_ns5.data, 2); % last recorded sample column

% Session calibration params, parsed once from the pre-trial setup text
% (see CALIBRATION NOTE above). Note this uses sync_time, not t_start, as
% the cutoff - the setup text is sent before the sync pulse, not necessarily
% before the (now earlier) stream start boundary.
predatcodes = digcodes(digcodes(:, 3) < sync_time, :);
tempdata.text = char(predatcodes(predatcodes(:, 2) >= 256 & predatcodes(:, 2) < 512, 2) - 256)';
if isempty(tempdata.text)
    tempdata = prev_tempdata;
else
    tempdata = getDatParams(tempdata);
end
params.block = tempdata.params.trial;
if isfield(params.block, 'reactionTime')
    params.block.crossingTime = params.block.reactionTime;
    params.block = rmfield(params.block, 'reactionTime');
end

% One bulk extraction + downsample over the entire continuous range
eyedeg = extract_eye_pupil_diode(out_ns5, samp_start:samp_end, params, eye_channels, [], []);
eyedata = eyedeg';
eye_times = (0:size(eyedata, 1)-1)' / 1000;
N = size(eyedata, 1);

% Digital event codes, one per 1 ms bin (excludes the 256-511 text-payload
% range, same filter format_datTrials uses to build trialcodes)
is_event_code = digcodes(:, 2) < 256 | (digcodes(:, 2) >= 1000 & digcodes(:, 2) <= 32000);
event_codes = digcodes(is_event_code & digcodes(:, 3) >= t_start & digcodes(:, 3) <= t_end, :);
code_bins = round((event_codes(:, 3) - t_start) * 1000) + 1;
valid = code_bins >= 1 & code_bins <= N;

codes_col = NaN(N, 1);
codes_col(code_bins(valid)) = event_codes(valid, 2); % last code wins if two land in the same bin

% Trial numbers: trial k starts at the (2k-1)-th sync pulse and ends at
% its own end-trial (255) code
nTrials = min(floor(numel(syncTimes) / 2), numel(endTimes));
trialnum_col = NaN(N, 1);
for k = 1:nTrials
    trial_start_time = syncTimes(2*k - 1);
    trial_end_time = endTimes(k);
    bin_start = max(1, round((trial_start_time - t_start) * 1000) + 1);
    bin_end = min(N, round((trial_end_time - t_start) * 1000) + 1);
    trialnum_col(bin_start:bin_end) = k;
end

trial_events = [trialnum_col, codes_col];

end

function samp = time_to_ns5_sample(t, out_ns5, Fs)
%TIME_TO_NS5_SAMPLE  Map an absolute time (sec) to its column index in out_ns5.data.
%   out_ns5.data only stores samples for recorded epochs, so gaps between
%   epochs (if any) must be skipped when converting an absolute sample
%   number into an index into out_ns5.data.
nEpochs = size(out_ns5.hdr.timeStamps, 2);
past_epochEnd = 0;
native_samp = round(t * Fs);
for epoch = 1:nEpochs
    epochStart = out_ns5.hdr.timeStamps(1, epoch);
    epochEnd = out_ns5.hdr.timeStamps(2, epoch);
    epochStart_samp = past_epochEnd + 1;
    epochEnd_samp = epochStart_samp + (epochEnd - epochStart) - 1;

    if native_samp >= epochStart && native_samp <= epochEnd
        samp = epochStart_samp + (native_samp - epochStart);
        return
    end
    past_epochEnd = epochEnd_samp;
end
error('extract_eyeStream:timeOutsideEpochs', 'Time %.4f sec does not fall within any recorded epoch.', t);
end
