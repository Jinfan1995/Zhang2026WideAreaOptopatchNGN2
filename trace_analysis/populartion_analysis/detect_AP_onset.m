function [onset_idx, onset_info] = detect_AP_onset(trace, peak_idx, cfg)
% DETECT_AP_ONSET
% Robust AP onset detection within a restricted pre-peak window.
%
% INPUT
% -------------------------------------------------------------------------
% trace
%       1-D AP waveform.
%
% peak_idx
%       Index of the AP peak within trace.
%
% cfg
%       analysis_config.onset structure containing:
%
%       cfg.smooth_window_frames
%       cfg.search_prepeak_frames
%       cfg.dvdt_fraction
%       cfg.min_sustain_frames
%       cfg.d2_refine_radius_frames
%
% OUTPUT
% -------------------------------------------------------------------------
% onset_idx
%       Detected onset index within trace.
%
% onset_info.fail
%       true when onset could not be identified.
%
% onset_info.used_fallback
%       true when the sustained-upstroke criterion failed and the
%       derivative-based fallback was used.
%
% onset_info.dvdt_threshold
%       First-derivative threshold used for onset detection.
%
% onset_info.max_dvdt
%       Maximum positive first derivative in the search window.
%
% LOGIC
% -------------------------------------------------------------------------
% 1. Lightly smooth the trace for onset detection only.
%
% 2. First derivative is the FORWARD neighboring-frame difference:
%
%       d1(i) = trace(i+1) - trace(i)
%
%    Therefore:
%       d1(i) > 0  -> signal rises from frame i to i+1
%       d1(i) < 0  -> signal falls from frame i to i+1
%
% 3. Second derivative is the CENTERED second difference:
%
%       d2(i) = trace(i+1) - 2*trace(i) + trace(i-1)
%
%    A large positive d2 indicates strong upward curvature centered
%    at frame i.
%
% 4. Search only within:
%
%       peak_idx - cfg.search_prepeak_frames : peak_idx - 1
%
%    With search_prepeak_frames = 3, onset can only be peak-3,
%    peak-2, or peak-1.
%
% 5. Find the earliest frame in this window whose forward derivative
%    is sufficiently positive and sustained toward the peak.
%
% 6. Refine locally using the maximum positive second derivative,
%    while requiring d1 > 0 so a downward-going frame cannot be chosen.
%
% IMPORTANT
% -------------------------------------------------------------------------
% Smoothing is used ONLY for onset localization. The original trace should
% still be used for AP height, width, AHP, and waveform measurements.
% -------------------------------------------------------------------------

%% ========================================================================
% 1. Prepare input
% =========================================================================
trace = double(trace(:)');
n = numel(trace);

onset_idx = NaN;

onset_info = struct;
onset_info.fail = true;
onset_info.used_fallback = false;
onset_info.dvdt_threshold = NaN;
onset_info.max_dvdt = NaN;

% Basic validation
if n < 3 || ...
        ~isscalar(peak_idx) || ...
        ~isfinite(peak_idx) || ...
        peak_idx < 2 || ...
        peak_idx > n
    return
end

peak_idx = round(peak_idx);

%% ========================================================================
% 2. Light smoothing
% =========================================================================
trace_smooth = smoothdata( ...
    trace, ...
    'gaussian', ...
    cfg.smooth_window_frames);

%% ========================================================================
% 3. First derivative: forward neighboring-frame difference
% =========================================================================
%
% d1(i) describes what happens immediately AFTER frame i.
%
% Example:
%       trace(12) < trace(11)  -> d1(11) < 0
%       trace(13) > trace(12)  -> d1(12) > 0
%
d1 = nan(1,n);

d1(1:n-1) = ...
    trace_smooth(2:n) - ...
    trace_smooth(1:n-1);

%% ========================================================================
% 4. Second derivative: centered second difference
% =========================================================================
%
% d2(i) is centered at frame i:
%
%       trace(i+1) - 2*trace(i) + trace(i-1)
%
d2 = nan(1,n);

d2(2:n-1) = ...
    trace_smooth(3:n) ...
    - 2 .* trace_smooth(2:n-1) ...
    + trace_smooth(1:n-2);

%% ========================================================================
% 5. Define allowed onset-search region
% =========================================================================
search_start = max(2, ...
    peak_idx - cfg.search_prepeak_frames);

search_end = min( ...
    peak_idx - 1, ...
    n - 1);

if search_start > search_end
    return
end

search_idx = search_start:search_end;

%% ========================================================================
% 6. Determine first-derivative threshold
% =========================================================================
max_d1 = max(d1(search_idx), [], 'omitnan');

onset_info.max_dvdt = max_d1;

if ~isfinite(max_d1) || max_d1 <= 0
    return
end

dvdt_threshold = ...
    cfg.dvdt_fraction * max_d1;

onset_info.dvdt_threshold = dvdt_threshold;

%% ========================================================================
% 7. Find first sustained positive upstroke
% =========================================================================
onset_candidate = NaN;

for k = search_idx

    % Normally require cfg.min_sustain_frames consecutive positive
    % forward differences.
    %
    % Near the AP peak, fewer pre-peak transitions may remain.
    % In that case use all transitions still available up to the peak.
    k_end = min( ...
        k + cfg.min_sustain_frames - 1, ...
        search_end);

    if all( ...
            isfinite(d1(k:k_end)) & ...
            d1(k:k_end) >= dvdt_threshold)

        onset_candidate = k;
        break
    end
end

%% ========================================================================
% 8. Fallback if sustained criterion was not met
% =========================================================================
if isnan(onset_candidate)

    % Find all frames that are followed by an actual increase.
    valid_idx = search_idx( ...
        isfinite(d1(search_idx)) & ...
        d1(search_idx) > 0);

    % ---------------------------------------------------------------
    % Case 1: exactly one rising frame
    % ---------------------------------------------------------------
    % If there is only one positive-going transition before the peak,
    % define that frame directly as the onset.
    if numel(valid_idx) == 1

        onset_idx = valid_idx(1);

        onset_info.used_fallback = true;
        onset_info.fail = false;

        return
    end

    % ---------------------------------------------------------------
    % Case 2: multiple rising frames
    % ---------------------------------------------------------------
    % Choose the one with the strongest upward curvature.
    if numel(valid_idx) > 1

        [~, temp_idx] = ...
            max(d2(valid_idx), [], 'omitnan');

        onset_idx = valid_idx(temp_idx);

        onset_info.used_fallback = true;
        onset_info.fail = false;

        return
    end

    % ---------------------------------------------------------------
    % Case 3: no positive-going frame
    % ---------------------------------------------------------------
    % Do not fail onset detection. If no frame has d1 > 0, choose the
    % frame with the largest forward difference (the least-negative /
    % most-positive transition) within the allowed pre-peak window.
    [~, temp_idx] = ...
        max(d1(search_idx), [], 'omitnan');

    onset_idx = search_idx(temp_idx);

    onset_info.used_fallback = true;
    onset_info.fail = false;

    return
end

%% ========================================================================
% 9. Refine onset using centered second derivative
% =========================================================================
radius = cfg.d2_refine_radius_frames;

refine_idx = ...
    max(search_start, onset_candidate-radius) : ...
    min(search_end, onset_candidate+radius);

% Critical constraint:
% only frames followed by a positive-going transition can be selected.
refine_idx = refine_idx( ...
    isfinite(d1(refine_idx)) & ...
    d1(refine_idx) > 0);

if isempty(refine_idx)

    onset_idx = onset_candidate;

else

    % Select strongest upward curvature within the local neighborhood.
    [~, temp_idx] = ...
        max(d2(refine_idx), [], 'omitnan');

    onset_idx = refine_idx(temp_idx);

end

%% ========================================================================
% 10. Final validation
% =========================================================================
if ~isfinite(onset_idx) || ...
        onset_idx < search_start || ...
        onset_idx > search_end || ...
        onset_idx >= peak_idx

    onset_idx = NaN;
    return
end

onset_info.fail = false;

end
