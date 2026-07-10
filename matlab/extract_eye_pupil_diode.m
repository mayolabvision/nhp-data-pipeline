function [eyedata, pupil, diode] = extract_eye_pupil_diode(out_ns5, samp_range, params, eye_channels, pupil_channel, diode_channel)
%EXTRACT_EYE_PUPIL_DIODE  Extract and downsample eye/pupil/diode signals over one sample range.
%   [eyedata, pupil, diode] = extract_eye_pupil_diode(out_ns5, samp_range, params, eye_channels, pupil_channel, diode_channel)
%
%   INPUT
%     out_ns5        : ns5 struct with raw 30 kHz data (out_ns5.data)
%     samp_range      : sample indices (into out_ns5.data) to extract
%     params           : calibration params struct passed through to eye2deg
%     eye_channels    : row indices of the eye (HE/VE) channels
%     pupil_channel   : row index of the pupil channel, or [] to skip
%     diode_channel   : row index of the diode channel, or [] to skip
%
%   OUTPUT
%     eyedata : [2 x N] eye position in degrees, downsampled to 1 kHz
%     pupil   : [1 x N] pupil signal downsampled to 1 kHz, or [] if pupil_channel is empty
%     diode   : [1 x N] diode signal downsampled to 1 kHz, or [] if diode_channel is empty

eyes = out_ns5.data(eye_channels, samp_range);
eyes_1khz = downsample(eyes', 30)'; % Downsample to 1 kHz
[eyedata, ~] = eye2deg(eyes_1khz(1:2, :), params); % Convert to degrees

pupil = [];
if ~isempty(pupil_channel)
    p = out_ns5.data(pupil_channel, samp_range);
    pupil = downsample(p', 30)'; % Downsample to 1 kHz
end

diode = [];
if ~isempty(diode_channel)
    d = out_ns5.data(diode_channel, samp_range);
    diode = downsample(d', 30)'; % Downsample to 1 kHz
end
end
