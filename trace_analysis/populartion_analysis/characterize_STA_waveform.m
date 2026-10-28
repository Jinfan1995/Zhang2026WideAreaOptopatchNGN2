function STA_features = characterize_STA_waveform(STA_trace, fps, peak_idx, onset_cfg)
% CHARACTERIZE_STA_WAVEFORM
% Characterize one peak-aligned STA waveform using the same definitions as
% Layer-0 AP shape extraction.
%
% INPUT
%   STA_trace   1 x n waveform, peak aligned.
%   fps         acquisition frame rate.
%   peak_idx    known AP peak index in STA_trace.
%   onset_cfg   onset settings passed to detect_AP_onset.
%
% OUTPUT fields
%   available
%   onset_idx
%   width20_ms
%   width50_ms
%   width80_ms
%   width20_rise_ms
%   width20_fall_ms
%   width50_rise_ms
%   width50_fall_ms
%   width80_rise_ms
%   width80_fall_ms
%   width20_rise_frames
%   width20_fall_frames
%   width50_rise_frames
%   width50_fall_frames
%   width80_rise_frames
%   width80_fall_frames
%   height
%   AHP
%   AHP_delay_ms
%   upstroke
%   downstroke
%
% Conventions
%   AHP = post-peak minimum - onset value (usually negative)
%   upstroke/downstroke = normalized-signal units per millisecond

STA_features = local_empty_features();

STA_trace = double(STA_trace(:)');

if isempty(STA_trace) || ...
        ~isscalar(fps) || ~isfinite(fps) || fps <= 0 || ...
        ~isscalar(peak_idx) || ~isfinite(peak_idx)
    return
end

peak_idx = round(double(peak_idx));
if peak_idx < 2 || peak_idx > numel(STA_trace)
    return
end

if any(~isfinite(STA_trace))
    return
end

[onset_idx, onset_info] = detect_AP_onset( ...
    STA_trace, ...
    peak_idx, ...
    onset_cfg);

if onset_info.fail || ~isfinite(onset_idx)
    return
end

onset_idx = round(double(onset_idx));
if onset_idx < 1 || onset_idx >= peak_idx
    return
end

onset_value = STA_trace(onset_idx);
peak_value = STA_trace(peak_idx);
pulse_height = peak_value - onset_value;

if ~isfinite(pulse_height) || pulse_height <= 0
    return
end

[width20_frames,width20_rise_frames,width20_fall_frames] = ...
    local_width_at_level( ...
        STA_trace, onset_idx, peak_idx, ...
        onset_value + 0.20*pulse_height);

[width50_frames,width50_rise_frames,width50_fall_frames] = ...
    local_width_at_level( ...
        STA_trace, onset_idx, peak_idx, ...
        onset_value + 0.50*pulse_height);

[width80_frames,width80_rise_frames,width80_fall_frames] = ...
    local_width_at_level( ...
        STA_trace, onset_idx, peak_idx, ...
        onset_value + 0.80*pulse_height);

post_peak = STA_trace(peak_idx:end);
[AHP_min_value,AHP_local_idx] = min(post_peak,[],'omitnan');

if isempty(AHP_local_idx) || ~isfinite(AHP_min_value)
    AHP_value = NaN;
    AHP_delay_ms = NaN;
else
    AHP_idx = peak_idx + AHP_local_idx - 1;
    AHP_value = AHP_min_value - onset_value;
    AHP_delay_ms = (AHP_idx - peak_idx) ./ fps .* 1000;
end

% Neighboring-frame forward derivative, converted to /ms.
d1_per_ms = nan(size(STA_trace));
d1_per_ms(1:end-1) = diff(STA_trace) .* fps ./ 1000;

upstroke_start = onset_idx;
upstroke_end = min(peak_idx-1,numel(d1_per_ms)-1);
if upstroke_end >= upstroke_start
    upstroke_value = max( ...
        d1_per_ms(upstroke_start:upstroke_end), ...
        [],'omitnan');
else
    upstroke_value = NaN;
end

downstroke_start = peak_idx;
downstroke_end = min(peak_idx+20,numel(d1_per_ms)-1);
if downstroke_end >= downstroke_start
    downstroke_value = min( ...
        d1_per_ms(downstroke_start:downstroke_end), ...
        [],'omitnan');
else
    downstroke_value = NaN;
end

STA_features.available = true;
STA_features.onset_idx = onset_idx;
STA_features.width20_ms = width20_frames ./ fps .* 1000;
STA_features.width50_ms = width50_frames ./ fps .* 1000;
STA_features.width80_ms = width80_frames ./ fps .* 1000;
STA_features.width20_rise_ms = width20_rise_frames ./ fps .* 1000;
STA_features.width20_fall_ms = width20_fall_frames ./ fps .* 1000;
STA_features.width50_rise_ms = width50_rise_frames ./ fps .* 1000;
STA_features.width50_fall_ms = width50_fall_frames ./ fps .* 1000;
STA_features.width80_rise_ms = width80_rise_frames ./ fps .* 1000;
STA_features.width80_fall_ms = width80_fall_frames ./ fps .* 1000;

% Raw interpolated timing in frames. These are durations relative to the
% known AP peak, so values can be fractional because threshold crossings
% are linearly interpolated between neighboring samples.
STA_features.width20_rise_frames = width20_rise_frames;
STA_features.width20_fall_frames = width20_fall_frames;
STA_features.width50_rise_frames = width50_rise_frames;
STA_features.width50_fall_frames = width50_fall_frames;
STA_features.width80_rise_frames = width80_rise_frames;
STA_features.width80_fall_frames = width80_fall_frames;

STA_features.height = pulse_height;
STA_features.AHP = AHP_value;
STA_features.AHP_delay_ms = AHP_delay_ms;
STA_features.upstroke = upstroke_value;
STA_features.downstroke = downstroke_value;

end


function S = local_empty_features()
S = struct( ...
    'available',false, ...
    'onset_idx',NaN, ...
    'width20_ms',NaN, ...
    'width50_ms',NaN, ...
    'width80_ms',NaN, ...
    'width20_rise_ms',NaN, ...
    'width20_fall_ms',NaN, ...
    'width50_rise_ms',NaN, ...
    'width50_fall_ms',NaN, ...
    'width80_rise_ms',NaN, ...
    'width80_fall_ms',NaN, ...
    'width20_rise_frames',NaN, ...
    'width20_fall_frames',NaN, ...
    'width50_rise_frames',NaN, ...
    'width50_fall_frames',NaN, ...
    'width80_rise_frames',NaN, ...
    'width80_fall_frames',NaN, ...
    'height',NaN, ...
    'AHP',NaN, ...
    'AHP_delay_ms',NaN, ...
    'upstroke',NaN, ...
    'downstroke',NaN);
end


function [width_frames,rise_to_peak_frames,peak_to_fall_frames] = ...
    local_width_at_level(trace,onset_idx,peak_idx,threshold)
% Same linearly interpolated threshold-crossing definition used in
% traces_extract_features.

width_frames = NaN;
rise_to_peak_frames = NaN;
peak_to_fall_frames = NaN;

if onset_idx >= peak_idx
    return
end

% Rising edge: last crossing before peak.
left_trace = trace(onset_idx:peak_idx);
left_idx = onset_idx:peak_idx;
left_crossing = ...
    (left_trace(1:end-1)-threshold) .* ...
    (left_trace(2:end)-threshold) <= 0;
left_cross = find(left_crossing,1,'last');

if isempty(left_cross)
    return
end

x1 = left_idx(left_cross);
x2 = left_idx(left_cross+1);
y1 = trace(x1);
y2 = trace(x2);

if ~isfinite(y1) || ~isfinite(y2) || y2 == y1
    return
end

left_position = x1 + (threshold-y1) ./ (y2-y1);

% Falling edge: first crossing after peak.
right_trace = trace(peak_idx:end);
right_idx = peak_idx:numel(trace);
right_crossing = ...
    (right_trace(1:end-1)-threshold) .* ...
    (right_trace(2:end)-threshold) <= 0;
right_cross = find(right_crossing,1,'first');

if isempty(right_cross)
    return
end

x1 = right_idx(right_cross);
x2 = right_idx(right_cross+1);
y1 = trace(x1);
y2 = trace(x2);

if ~isfinite(y1) || ~isfinite(y2) || y2 == y1
    return
end

right_position = x1 + (threshold-y1) ./ (y2-y1);

rise_to_peak_frames = peak_idx - left_position;
peak_to_fall_frames = right_position - peak_idx;
width_frames = right_position - left_position;

end
