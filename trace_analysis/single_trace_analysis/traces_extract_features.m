function [features, spike_data, feature_qc] = traces_extract_features(traces_nm, kernel, protocol, analysis_config)
% TRACES_EXTRACT_FEATURES
% Protocol-driven Optopatch feature extraction.
% Layer-0 spike/event extraction plus Layer-1 single-cell features.
%
% Performance changes in this version:
%   - STA waveforms are not stored/calculated here
%   - Layer-1 AP/firing/subthreshold/recovery summaries are calculated here
%   - fast closed-form ordinary least-squares instead of repeated POLYFIT
%   - HISTCOUNTS for whole-trace AP-count bins
%   - vectorized step-height and ramp-area calculations
%   - post-blue-off decay windows/times precomputed once from protocol
%   - protocol epoch boundaries precomputed outside the cell loop
%   - no full-trace conversion to double inside the cell loop
%   - kernel converted once to double before convolution (required by imfilter)
%   - FINDPEAKS uses MinPeakHeight directly
%   - batch vectorized offset-exponential tau fitting (no FIT objects)
%   - 8 AP-adaptation regressions solved together per stimulus epoch
%   - vectorized local peak refinement around detected spikes
%   - MOVMEDIAN and serial cell processing are intentionally unchanged
% INPUT
% -------------------------------------------------------------------------
% traces_nm
%       nCells x nFrames normalized traces.
% kernel
%       nCells x kernelLength STA kernels.
% protocol
%       Stimulus protocol structure.
%       Required:
%       protocol.fps
%       protocol.baseline.frames
%       protocol.analysis.median_filter_size
%       protocol.epoch(i).name
%       protocol.epoch(i).type
%       protocol.epoch(i).start_frame
%       protocol.epoch(i).end_frame
%       epoch.type should include:
%           "step"
%           "ramp"
%           "blue_off"
% analysis_config
%       Feature extraction parameters.
% OUTPUT
% -------------------------------------------------------------------------
% features
%       nCells x nFeatures MATLAB table.
%       One row per cell containing Layer-1 AP adaptation, firing,
%       subthreshold, rheobase, and post-blue-off recovery features.
% spike_data
%       nCells x 1 Layer-0 struct containing event-aligned spike data.
%       spike_frames and onset_frames are int16. Unknown onset/AHP-delay
%       frame values use -1. Width measurements are floating point because
%       crossing positions are linearly interpolated.
%
%       Upstroke/downstroke convention:
%           d1(i) = [trace(i+1) - trace(i)] / dt_ms
%
%       where:
%           dt_ms = 1000 / protocol.fps
%
%       Therefore spike_data.upstroke and spike_data.downstroke are stored
%       in normalized-signal units per millisecond.
%
% feature_qc
%       QC flags. No. cells are deleted.
% -------------------------------------------------------------------------

%% ========================================================================
% 1. Basic information
% =========================================================================
[nCells, nFrames] = size(traces_nm);
fps = protocol.fps;
baseline_frames = unique(protocol.baseline.frames(:));
baseline_frames = baseline_frames(baseline_frames >= 1 & baseline_frames <= nFrames);
if isempty(baseline_frames)
    error('No valid baseline frames were supplied.');
end
if size(kernel,1) ~= nCells
    error('kernel and traces_nm must contain the same number of cells.');
end

%% ========================================================================
% 2. Analysis configuration
% =========================================================================
median_filter_size = protocol.analysis.median_filter_size;
cfg = analysis_config.features;

spike_SNR_th = cfg.spike_SNR_threshold;
min_peak_distance = cfg.min_peak_distance_frames;
peak_refine_radius = cfg.peak_refine_radius_frames;
edge_pre = cfg.spike_edge_pre_frames;
edge_post = cfg.spike_edge_post_frames;
extract_pre = cfg.spike_extract_pre_frames;
extract_post = cfg.spike_extract_post_frames;
onset_cfg = analysis_config.onset;
close_spike_threshold = cfg.close_spike_fail_frames;

% Layer-1 defaults. These can be overridden in analysis_config.features.
if isfield(cfg,'AP_count_bin_frames')
    AP_count_bin_frames = cfg.AP_count_bin_frames;
else
    AP_count_bin_frames = 100;
end

if isfield(cfg,'min_APs_for_adaptation_fit')
    min_APs_for_adaptation_fit = cfg.min_APs_for_adaptation_fit;
else
    min_APs_for_adaptation_fit = 3;
end

if isfield(cfg,'min_spikes_for_frequency_fit')
    min_spikes_for_frequency_fit = cfg.min_spikes_for_frequency_fit;
else
    min_spikes_for_frequency_fit = 4;
end

if isfield(cfg,'decay_start_offset_frames')
    decay_start_offset_frames = cfg.decay_start_offset_frames;
else
    decay_start_offset_frames = 10;
end

if isfield(cfg,'decay_end_offset_frames')
    decay_end_offset_frames = cfg.decay_end_offset_frames;
else
    decay_end_offset_frames = 200;
end

% Fast decay-tau solver. The model is still:
%       y = a*exp(-t/tau) + c
% but tau is searched on a log-spaced grid while a and c are solved
% analytically for every tau candidate. A 3-point parabolic refinement in
% log(tau) improves precision without a nonlinear FIT call.
if isfield(cfg,'decay_tau_min_ms')
    decay_tau_min_ms = cfg.decay_tau_min_ms;
else
    decay_tau_min_ms = max(1000 ./ fps, 1);
end

% The decay window is short (~10 to 200 frames after blue off), so very
% long tau values are poorly identifiable. Cap the search at 1000 ms.
% Values that still land near either search boundary are rejected below.
decay_tau_hard_max_ms = 1000;

if isfield(cfg,'decay_tau_max_ms')
    decay_tau_max_ms = min(double(cfg.decay_tau_max_ms), decay_tau_hard_max_ms);
else
    decay_tau_max_ms = decay_tau_hard_max_ms;
end

if isfield(cfg,'decay_tau_grid_points')
    decay_tau_grid_points = cfg.decay_tau_grid_points;
else
    decay_tau_grid_points = 192;
end

if decay_tau_min_ms <= 0 || ...
        decay_tau_max_ms <= decay_tau_min_ms || ...
        decay_tau_grid_points < 16
    error('Invalid fast decay-tau grid configuration.');
end

decay_tau_grid_ms = logspace( ...
    log10(decay_tau_min_ms), ...
    log10(decay_tau_max_ms), ...
    round(decay_tau_grid_points));

if analysis_config.progress.enable
    progress_interval = max(1, round(nCells * analysis_config.progress.percent_interval / 100));
else
    progress_interval = inf;
end

%% ========================================================================
% 3. Read protocol epochs
% =========================================================================
if ~isfield(protocol,'epoch') || isempty(protocol.epoch)
    error('protocol.epoch is missing.');
end

% Sort epochs according to their actual time.
epoch_start = [protocol.epoch.start_frame];
[~, epoch_order] = sort(epoch_start);
epochs = protocol.epoch(epoch_order);
epoch_types = string({epochs.type});

step_idx = find(epoch_types == "step");
ramp_idx = find(epoch_types == "ramp");
blue_off_idx = find(epoch_types == "blue_off");
stim_idx = find(epoch_types == "step" | epoch_types == "ramp");

nSteps = numel(step_idx);
nRamps = numel(ramp_idx);
nStim = numel(stim_idx);

step_start = [epochs(step_idx).start_frame];
step_end = [epochs(step_idx).end_frame];
ramp_start = [epochs(ramp_idx).start_frame];
ramp_end = [epochs(ramp_idx).end_frame];
stim_start = [epochs(stim_idx).start_frame];
stim_end = [epochs(stim_idx).end_frame];

% For each blue-on stimulus, identify the first following blue-off epoch.
decay_blue_epoch = nan(1,nStim);
for iStim = 1:nStim
    candidates = blue_off_idx([epochs(blue_off_idx).start_frame] > stim_end(iStim));
    if ~isempty(candidates)
        decay_blue_epoch(iStim) = candidates(1);
    end
end

% Precompute post-blue-off decay fit windows and time axes once.
% These depend only on the protocol, not on cell identity.
decay_fit_frames = cell(1,nStim);
decay_time_sec = cell(1,nStim);

for iStim = 1:nStim
    if ~isfinite(decay_blue_epoch(iStim))
        continue
    end

    blue_e = decay_blue_epoch(iStim);
    blue_off_start = epochs(blue_e).start_frame;
    blue_off_end = epochs(blue_e).end_frame;

    fit_start = blue_off_start + decay_start_offset_frames;
    fit_end = min([ ...
        blue_off_start + decay_end_offset_frames, ...
        blue_off_end, ...
        nFrames]);

    if fit_end <= fit_start
        continue
    end

    decay_fit_frames{iStim} = fit_start:fit_end;
    decay_time_sec{iStim} = ...
        ((double(decay_fit_frames{iStim}) - double(blue_off_start)) ./ fps)';
end

%% ========================================================================
% 4. Generate whole-trace AP-count bins
% =========================================================================
% Layer-1 firing trajectory is represented as AP COUNTS, not frequency,
% in consecutive fixed-width bins spanning the entire trace.
count_bin_start = 1:AP_count_bin_frames:nFrames;
count_bin_end = min(count_bin_start + AP_count_bin_frames - 1, nFrames);
nCountBins = numel(count_bin_start);

% HISTCOUNTS edges give bins [1,101), [101,201), ... so integer
% spike-frame indices are counted exactly once, including the final short bin.
count_bin_edges = double([count_bin_start, nFrames + 1]);

%% ========================================================================
% 5. Preprocess traces
% =========================================================================
lowpass_traces = movmedian(traces_nm, median_filter_size, 2, 'omitnan');
highpass_traces = traces_nm - lowpass_traces;

%% ========================================================================
% 6. Convolve high-pass traces with STA kernels
% =========================================================================
traces_conv = nan(size(traces_nm), 'like', traces_nm);

% MATLAB imfilter requires the filter kernel h to be DOUBLE.
% Convert the kernel once outside the cell loop to avoid repeated casts.
kernel_work = double(kernel);

for iCell = 1:nCells
    if all(~isfinite(traces_nm(iCell,:)))
        continue
    end
    if any(~isfinite(kernel_work(iCell,:)))
        continue
    end

    traces_conv(iCell,:) = ...
        imfilter( ...
            highpass_traces(iCell,:), ...
            kernel_work(iCell,:));
end

%% ========================================================================
% 7. Estimate noise
% =========================================================================
noise_std = std(traces_conv(:,baseline_frames), 0, 2, 'omitnan');

% Input-validity mask used both in the cell loop and for vectorized
% subthreshold calculations. This preserves the previous behavior: cells
% with invalid inputs keep Layer-1 values as NaN.
valid_cell_input = ...
    ~all(~isfinite(traces_nm),2) & ...
    isfinite(noise_std) & ...
    noise_std > 0;

%% ========================================================================
% 8. Vectorized subthreshold features
% =========================================================================
% These operations depend only on low-pass traces and fixed protocol
% windows, so calculate all cells together instead of inside the cell loop.
step_height = nan(nCells,nSteps);
for iStep = 1:nSteps
    s = step_start(iStep);
    eFrame = step_end(iStep);
    step_height(:,iStep) = ...
        median(lowpass_traces(:,s:eFrame), 2, 'omitnan');
end
step_height(~valid_cell_input,:) = NaN;

ramp_area = nan(nCells,nRamps);
for iRamp = 1:nRamps
    s = ramp_start(iRamp);
    eFrame = ramp_end(iRamp);
    ramp_area(:,iRamp) = ...
        trapz(lowpass_traces(:,s:eFrame), 2);
end
ramp_area(~valid_cell_input,:) = NaN;

%% ========================================================================
% 9. Preallocate remaining Layer-1 feature matrices
% =========================================================================
% Whole-trace AP counts
AP_count_bin = zeros(nCells,nCountBins,'uint16');

% Overall AP phenotype across all valid events in each cell
AP_n_spikes = zeros(nCells,1,'uint16');
AP_width20_median_ms = nan(nCells,1);
AP_width50_median_ms = nan(nCells,1);
AP_width80_median_ms = nan(nCells,1);
AP_width50_std_ms = nan(nCells,1);
AP_width50_CV = nan(nCells,1);
AP_width50_rise_median_ms = nan(nCells,1);
AP_width50_fall_median_ms = nan(nCells,1);
AP_width20_rise_median_ms = nan(nCells,1);
AP_width20_fall_median_ms = nan(nCells,1);
AP_width80_rise_median_ms = nan(nCells,1);
AP_width80_fall_median_ms = nan(nCells,1);
AP_height_median = nan(nCells,1);
AP_height_CV = nan(nCells,1);
AP_AHP_median = nan(nCells,1);
AP_AHP_delay_median_ms = nan(nCells,1);
AP_upstroke_median = nan(nCells,1);
AP_downstroke_median = nan(nCells,1);

% STA-derived scalar AP phenotype. The STA waveform itself is NOT stored.
AP_STA_width_ms = nan(nCells,1);
AP_STA_height = nan(nCells,1);
AP_STA_AHP = nan(nCells,1);
AP_STA_upstroke = nan(nCells,1);
AP_STA_downstroke = nan(nCells,1);

% Subthreshold/excitability feature not vectorized above
rheobase_frame = nan(nCells,nRamps);

% Firing timing and adaptation for each blue-on step/ramp
first_spike_delay_ms = nan(nCells,nStim);
first_second_ISI_ms = nan(nCells,nStim);
penultimate_last_ISI_ms = nan(nCells,nStim);
last_spike_delay_ms = nan(nCells,nStim);
frequency_time_slope_Hz_per_s = nan(nCells,nStim);
frequency_time_intercept_Hz = nan(nCells,nStim);
max_spike_rate_Hz = nan(nCells,nStim);

% AP waveform adaptation: y = slope * time_ms + intercept
AP_height_slope = nan(nCells,nStim);
AP_height_intercept = nan(nCells,nStim);
AP_width20_slope = nan(nCells,nStim);
AP_width20_intercept = nan(nCells,nStim);
AP_width50_slope = nan(nCells,nStim);
AP_width50_intercept = nan(nCells,nStim);
AP_width80_slope = nan(nCells,nStim);
AP_width80_intercept = nan(nCells,nStim);
AP_width20_rise_slope = nan(nCells,nStim);
AP_width20_rise_intercept = nan(nCells,nStim);
AP_width20_fall_slope = nan(nCells,nStim);
AP_width20_fall_intercept = nan(nCells,nStim);
AP_width50_rise_slope = nan(nCells,nStim);
AP_width50_rise_intercept = nan(nCells,nStim);
AP_width50_fall_slope = nan(nCells,nStim);
AP_width50_fall_intercept = nan(nCells,nStim);
AP_width80_rise_slope = nan(nCells,nStim);
AP_width80_rise_intercept = nan(nCells,nStim);
AP_width80_fall_slope = nan(nCells,nStim);
AP_width80_fall_intercept = nan(nCells,nStim);
AP_AHP_delay_slope = nan(nCells,nStim);
AP_AHP_delay_intercept = nan(nCells,nStim);
AP_AHP_depth_slope = nan(nCells,nStim);
AP_AHP_depth_intercept = nan(nCells,nStim);

% Post-blue-off exponential recovery:
%       y = a*exp(-t/tau) + c
% t is seconds after blue-off onset; fitting window starts at +10 frames
% and ends at +200 frames (or earlier if the blue-off epoch ends).
% Only tau (ms) is retained as a Layer-1 feature.
decay_tau_ms = nan(nCells,nStim);

%% ========================================================================
% 10. Spike data structure (Layer 0)
% =========================================================================
% One index = one refined spike event. All event-level vectors have the
% same length and refer to spike_data(iCell).spike_frames.
%
% Integer frame fields use int16. A value of -1 means unavailable for
% event properties such as onset_frames or AHP_delay_frames. Width
% measurements remain floating point because threshold crossings are
% linearly interpolated and can occur at fractional frame positions.
if nFrames > double(intmax('int16'))
    error(['nFrames exceeds int16 capacity. Use int32/uint16 for frame ' ...
        'indices before analyzing this protocol.']);
end

empty_spike = struct( ...
    'spike_frames', int16([]), ...
    'peak_aligned_traces', single([]), ...
    'onset_frames', int16([]), ...
    'width20_frames', single([]), ...
    'width20_rise_frames', single([]), ...
    'width20_fall_frames', single([]), ...
    'width50_frames', single([]), ...
    'width50_rise_frames', single([]), ...
    'width50_fall_frames', single([]), ...
    'width80_frames', single([]), ...
    'width80_rise_frames', single([]), ...
    'width80_fall_frames', single([]), ...
    'height', single([]), ...
    'AHP', single([]), ...
    'AHP_delay_frames', int16([]), ...
    'upstroke', single([]), ...
    'downstroke', single([]));

spike_data = ...
    repmat(empty_spike, nCells, 1);


%% ========================================================================
% 11. Feature QC
% =========================================================================
feature_qc.invalid_input = false(nCells,1);
feature_qc.no_spike = false(nCells,1);
feature_qc.close_spikes = false(nCells,1);
feature_qc.spike_shape_unavailable = false(nCells,1);
feature_qc.STA_onset_fail = false(nCells,1);
feature_qc.individual_onset_fail_count = zeros(nCells,1);
feature_qc.decay_fit_fail = false(nCells,nStim);

%% ------------------------------------------------------------
% Progress reporting
% ------------------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf('Extracting electrophysiology features\n');
fprintf('============================================================\n');
fprintf('Number of cells: %d\n', nCells);
fprintf('Number of epochs: %d\n', numel(protocol.epoch));

nStep = sum([protocol.epoch.type] == "step");
nRamp = sum([protocol.epoch.type] == "ramp");
nOff  = sum([protocol.epoch.type] == "blue_off");

fprintf('Step epochs:     %d\n', nStep);
fprintf('Ramp epochs:     %d\n', nRamp);
fprintf('Blue-off epochs: %d\n', nOff);
fprintf('------------------------------------------------------------\n');


t_start = tic;

%% ========================================================================
% 12. Fast batch post-blue-off decay-tau fitting
% =========================================================================
% The expensive general-purpose FIT call is removed from the cell loop.
% For each stimulus, all valid cells are fitted together using the same
% offset-exponential model:
%
%       y(t) = a*exp(-t/tau) + c
%
% For every candidate tau, a and c are ordinary least-squares parameters.
% This lets MATLAB evaluate all cells with matrix multiplication.
for iStim = 1:nStim
    fit_frames = decay_fit_frames{iStim};
    decay_t = decay_time_sec{iStim};

    if isempty(fit_frames) || isempty(decay_t)
        feature_qc.decay_fit_fail(valid_cell_input,iStim) = true;
        continue
    end

    decay_Y = double(lowpass_traces(:,fit_frames));

    [tau_this, ok_this] = ...
        local_fast_decay_tau_batch( ...
            decay_Y, ...
            decay_t, ...
            valid_cell_input, ...
            decay_tau_grid_ms);

    % Reject poorly identifiable fits that land at/near either edge of the
    % allowed tau search range. A boundary hit means the data do not contain
    % enough information to determine tau within the requested window.
    % Reject the first/last two grid positions rather than reporting an
    % arbitrary extreme tau value.
    if numel(decay_tau_grid_ms) >= 5
        tau_lower_reject = decay_tau_grid_ms(2);
        tau_upper_reject = decay_tau_grid_ms(end-1);

        boundary_hit = ...
            isfinite(tau_this) & ...
            (tau_this <= tau_lower_reject | ...
             tau_this >= tau_upper_reject);

        tau_this(boundary_hit) = NaN;
        ok_this(boundary_hit) = false;
    end

    decay_tau_ms(:,iStim) = tau_this;

    feature_qc.decay_fit_fail(valid_cell_input,iStim) = ...
        ~ok_this(valid_cell_input);
end

%% ========================================================================
% 13. Analyze each cell
% =========================================================================
for iCell = 1:nCells
    %% --------------------------------------------------------------------
    % 12A. Input QC
    % ---------------------------------------------------------------------
    if ~valid_cell_input(iCell)
        feature_qc.invalid_input(iCell) = true;
        continue
    end
    %% --------------------------------------------------------------------
    % 12B. Second-round spike detection
    % ---------------------------------------------------------------------
    % Keep the full trace in its existing datatype (typically single).
    detection_trace = traces_conv(iCell,:);
    detection_trace(~isfinite(detection_trace)) = 0;

    threshold = spike_SNR_th * noise_std(iCell);

    % MATLAB FINDPEAKS warns when MinPeakHeight exceeds every sample.
    % Treat this expected case directly as a valid non-spiking cell.
    if max(detection_trace,[],'omitnan') <= threshold
        local_spikes = [];
    else
        [~, local_spikes] = ...
            findpeaks( ...
                detection_trace, ...
                'MinPeakHeight', threshold, ...
                'MinPeakDistance', min_peak_distance);
    end
    local_spikes = local_spikes(:)';
    % If no spikes detected
    if isempty(local_spikes)
        feature_qc.no_spike(iCell) = true;
        % DO NOT continue. Keep processing the cell so:
        % firing frequency = 0
        % step height       = valid measurement
        % ramp area         = valid measurement
        % while spike-dependent features remain NaN.
    end
    %% --------------------------------------------------------------------
    % 12C. Refine spike position using original high-pass trace
    % ---------------------------------------------------------------------
    local_highpass = highpass_traces(iCell,:);

    if ~isempty(local_spikes)
        refine_offsets = -peak_refine_radius:peak_refine_radius;

        refine_idx = ...
            double(local_spikes(:)) + double(refine_offsets);

        refine_idx = ...
            max(1, min(nFrames, refine_idx));

        refine_values = ...
            reshape( ...
                local_highpass(refine_idx), ...
                size(refine_idx));

        [~, refine_choice] = ...
            max(refine_values, [], 2);

        linear_choice = ...
            sub2ind( ...
                size(refine_idx), ...
                (1:size(refine_idx,1))', ...
                refine_choice);

        refined_spikes = ...
            refine_idx(linear_choice)';

        % Remove exact duplicate locations created by refinement.
        local_spikes = unique(refined_spikes, 'stable');
    end

    %% --------------------------------------------------------------------
    % 12D. Close-spike QC
    % closely spaced detections may indicate mixed traces.
    % ---------------------------------------------------------------------
    if numel(local_spikes) > 1
        if min(diff(sort(local_spikes))) <= close_spike_threshold
            feature_qc.close_spikes(iCell) = true;
        end
    end

    %% ====================================================================
    % 13. Individual AP waveform measurements (Layer 0)
    % =====================================================================
    % Every refined spike receives exactly one event index. If a particular
    % waveform measurement cannot be obtained, its entry remains NaN (or -1
    % for integer frame fields) rather than shortening that feature vector.
    nLocalSpikes = numel(local_spikes);
    waveform_length = extract_pre + extract_post + 1;

    local_peak_aligned = ...
        nan(nLocalSpikes, waveform_length, 'single');

    local_onset = ...
        repmat(int16(-1), 1, nLocalSpikes);

    local_width20 = nan(1,nLocalSpikes,'single');
    local_width20_rise = nan(1,nLocalSpikes,'single');
    local_width20_fall = nan(1,nLocalSpikes,'single');

    local_width50 = nan(1,nLocalSpikes,'single');
    local_width50_rise = nan(1,nLocalSpikes,'single');
    local_width50_fall = nan(1,nLocalSpikes,'single');

    local_width80 = nan(1,nLocalSpikes,'single');
    local_width80_rise = nan(1,nLocalSpikes,'single');
    local_width80_fall = nan(1,nLocalSpikes,'single');

    local_height = nan(1,nLocalSpikes,'single');
    local_AHP = nan(1,nLocalSpikes,'single');
    local_AHP_delay = repmat(int16(-1),1,nLocalSpikes);
    local_upstroke = nan(1,nLocalSpikes,'single');
    local_downstroke = nan(1,nLocalSpikes,'single');

    peak_idx = extract_pre + 1;

    for iSpike = 1:nLocalSpikes
        spike_frame = local_spikes(iSpike);

        % Edge exclusion
        if spike_frame <= edge_pre
            continue
        end
        if spike_frame >= nFrames-edge_post
            continue
        end

        % Extract AP trace. spike_frame has already been refined to the
        % original high-pass peak, so peak_idx is the known AP-peak index.
        trace_start = spike_frame - extract_pre;
        trace_end = spike_frame + extract_post;
        if trace_start < 1 || trace_end > nFrames
            continue
        end

        single_peak_trace = ...
            double(traces_nm(iCell,trace_start:trace_end));

        if any(~isfinite(single_peak_trace))
            continue
        end

        % Store one PEAK-aligned waveform per detected event. Layer 1 can
        % later generate STA/median waveforms from these rows.
        local_peak_aligned(iSpike,:) = ...
            single(single_peak_trace);

        %% ---------------------------------------------------------------
        % Robust individual AP onset detection
        % ----------------------------------------------------------------
        [I, onset_info] = detect_AP_onset( ...
            single_peak_trace, ...
            peak_idx, ...
            onset_cfg);

        if onset_info.fail || ~isfinite(I)
            feature_qc.individual_onset_fail_count(iCell) = ...
                feature_qc.individual_onset_fail_count(iCell) + 1;
            continue
        end

        onset_frame = trace_start + I - 1;
        local_onset(iSpike) = int16(onset_frame);

        %% ---------------------------------------------------------------
        % Spike height
        % ----------------------------------------------------------------
        onset_value = single_peak_trace(I);
        peak_value = single_peak_trace(peak_idx);
        pulse_height = peak_value - onset_value;

        %% ---------------------------------------------------------------
        % AP widths at 20%, 50%, and 80% of AP height
        % ----------------------------------------------------------------
        [width20, rise20, fall20, ...
         width50, rise50, fall50, ...
         width80, rise80, fall80] = ...
            local_AP_widths_20_50_80( ...
                single_peak_trace, ...
                I, ...
                peak_idx, ...
                onset_value, ...
                pulse_height);

        %% ---------------------------------------------------------------
        % AHP magnitude and timing
        % ----------------------------------------------------------------
        % Search only after the AP peak so a pre-peak fluctuation cannot be
        % mislabeled as the AHP minimum.
        post_peak_trace = single_peak_trace(peak_idx:end);
        [AHP_min_value, AHP_local_idx] = ...
            min(post_peak_trace, [], 'omitnan');

        AHP_idx = peak_idx + AHP_local_idx - 1;
        AHP_value = AHP_min_value - onset_value;
        AHP_delay_frames = AHP_idx - peak_idx;

        %% ---------------------------------------------------------------
        % Upstroke / downstroke
        % ----------------------------------------------------------------
        % d1(i) = [trace(i+1)-trace(i)] / dt_ms
        % Units: normalized-signal units per millisecond.
        d1_per_ms = nan(size(single_peak_trace));
        d1_per_ms(1:end-1) = ...
            diff(single_peak_trace) .* fps ./ 1000;

        upstroke_start = I;
        upstroke_end = min(peak_idx - 1, numel(d1_per_ms) - 1);

        if upstroke_end >= upstroke_start
            upstroke_value = ...
                max(d1_per_ms(upstroke_start:upstroke_end), [], 'omitnan');
        else
            upstroke_value = NaN;
        end

        downstroke_start = peak_idx;
        downstroke_end = min(peak_idx + 20, numel(d1_per_ms) - 1);

        if downstroke_end >= downstroke_start
            downstroke_value = ...
                min(d1_per_ms(downstroke_start:downstroke_end), [], 'omitnan');
        else
            downstroke_value = NaN;
        end

        %% ---------------------------------------------------------------
        % Store event measurements at the SAME event index
        % ----------------------------------------------------------------
        local_width20(iSpike) = single(width20);
        local_width20_rise(iSpike) = single(rise20);
        local_width20_fall(iSpike) = single(fall20);

        local_width50(iSpike) = single(width50);
        local_width50_rise(iSpike) = single(rise50);
        local_width50_fall(iSpike) = single(fall50);

        local_width80(iSpike) = single(width80);
        local_width80_rise(iSpike) = single(rise80);
        local_width80_fall(iSpike) = single(fall80);

        local_height(iSpike) = single(pulse_height);
        local_AHP(iSpike) = single(AHP_value);
        local_AHP_delay(iSpike) = int16(AHP_delay_frames);
        local_upstroke(iSpike) = single(upstroke_value);
        local_downstroke(iSpike) = single(downstroke_value);
    end

    %% Store Layer-0 individual AP data
    spike_data(iCell).spike_frames = int16(local_spikes);
    spike_data(iCell).peak_aligned_traces = local_peak_aligned;
    spike_data(iCell).onset_frames = local_onset;

    spike_data(iCell).width20_frames = local_width20;
    spike_data(iCell).width20_rise_frames = local_width20_rise;
    spike_data(iCell).width20_fall_frames = local_width20_fall;

    spike_data(iCell).width50_frames = local_width50;
    spike_data(iCell).width50_rise_frames = local_width50_rise;
    spike_data(iCell).width50_fall_frames = local_width50_fall;

    spike_data(iCell).width80_frames = local_width80;
    spike_data(iCell).width80_rise_frames = local_width80_rise;
    spike_data(iCell).width80_fall_frames = local_width80_fall;

    spike_data(iCell).height = local_height;
    spike_data(iCell).AHP = local_AHP;
    spike_data(iCell).AHP_delay_frames = local_AHP_delay;
    spike_data(iCell).upstroke = local_upstroke;
    spike_data(iCell).downstroke = local_downstroke;

    %% ====================================================================
    % 14. Spike-shape availability QC
    % =====================================================================
    valid_peak_rows = all(isfinite(local_peak_aligned),2);
    if ~any(valid_peak_rows)
        feature_qc.spike_shape_unavailable(iCell) = true;
    end

    %% ====================================================================
    % 15. Overall Layer-1 AP phenotype + STA scalar features
    % =====================================================================
    AP_n_spikes(iCell) = uint16(nLocalSpikes);

    width20_ms_all = double(local_width20) ./ fps .* 1000;
    width50_ms_all = double(local_width50) ./ fps .* 1000;
    width80_ms_all = double(local_width80) ./ fps .* 1000;
    width20_rise_ms_all = double(local_width20_rise) ./ fps .* 1000;
    width20_fall_ms_all = double(local_width20_fall) ./ fps .* 1000;
    width50_rise_ms_all = double(local_width50_rise) ./ fps .* 1000;
    width50_fall_ms_all = double(local_width50_fall) ./ fps .* 1000;
    width80_rise_ms_all = double(local_width80_rise) ./ fps .* 1000;
    width80_fall_ms_all = double(local_width80_fall) ./ fps .* 1000;

    AHP_delay_ms_all = double(local_AHP_delay);
    AHP_delay_ms_all(AHP_delay_ms_all < 0) = NaN;
    AHP_delay_ms_all = AHP_delay_ms_all ./ fps .* 1000;

    AP_width20_median_ms(iCell) = median(width20_ms_all,'omitnan');
    AP_width50_median_ms(iCell) = median(width50_ms_all,'omitnan');
    AP_width80_median_ms(iCell) = median(width80_ms_all,'omitnan');
    AP_width50_std_ms(iCell) = local_nan_std(width50_ms_all);
    AP_width50_CV(iCell) = local_nan_cv(width50_ms_all);
    AP_width20_rise_median_ms(iCell) = median(width20_rise_ms_all,'omitnan');
    AP_width20_fall_median_ms(iCell) = median(width20_fall_ms_all,'omitnan');
    AP_width50_rise_median_ms(iCell) = median(width50_rise_ms_all,'omitnan');
    AP_width50_fall_median_ms(iCell) = median(width50_fall_ms_all,'omitnan');
    AP_width80_rise_median_ms(iCell) = median(width80_rise_ms_all,'omitnan');
    AP_width80_fall_median_ms(iCell) = median(width80_fall_ms_all,'omitnan');

    AP_height_median(iCell) = median(double(local_height),'omitnan');
    AP_height_CV(iCell) = local_nan_cv(double(local_height));
    AP_AHP_median(iCell) = median(double(local_AHP),'omitnan');
    AP_AHP_delay_median_ms(iCell) = median(AHP_delay_ms_all,'omitnan');
    AP_upstroke_median(iCell) = median(double(local_upstroke),'omitnan');
    AP_downstroke_median(iCell) = median(double(local_downstroke),'omitnan');

    % Build a robust median STA from valid PEAK-aligned event waveforms,
    % detect onset once on the STA, and store only scalar STA features.
    if any(valid_peak_rows)
        STA_peak_aligned = median(double(local_peak_aligned(valid_peak_rows,:)),1,'omitnan');
        [STA_onset_idx, STA_onset_info] = detect_AP_onset( ...
            STA_peak_aligned, peak_idx, onset_cfg);

        if STA_onset_info.fail || ~isfinite(STA_onset_idx) || STA_onset_idx >= peak_idx
            feature_qc.STA_onset_fail(iCell) = true;
        else
            STA_onset_value = STA_peak_aligned(STA_onset_idx);
            STA_peak_value = STA_peak_aligned(peak_idx);
            STA_height_value = STA_peak_value - STA_onset_value;
            STA_half_height = STA_onset_value + 0.5 * STA_height_value;

            [STA_width_frames,~,~] = local_width_at_level( ...
                STA_peak_aligned, STA_onset_idx, peak_idx, STA_half_height);

            STA_post_peak = STA_peak_aligned(peak_idx:end);
            STA_AHP_value = min(STA_post_peak,[],'omitnan') - STA_onset_value;

            STA_d1_per_ms = nan(size(STA_peak_aligned));
            STA_d1_per_ms(1:end-1) = diff(STA_peak_aligned) .* fps ./ 1000;

            STA_upstroke_start = STA_onset_idx;
            STA_upstroke_end = min(peak_idx-1,numel(STA_d1_per_ms)-1);
            if STA_upstroke_end >= STA_upstroke_start
                STA_upstroke_value = max( ...
                    STA_d1_per_ms(STA_upstroke_start:STA_upstroke_end),[], 'omitnan');
            else
                STA_upstroke_value = NaN;
            end

            STA_downstroke_start = peak_idx;
            STA_downstroke_end = min(peak_idx+20,numel(STA_d1_per_ms)-1);
            if STA_downstroke_end >= STA_downstroke_start
                STA_downstroke_value = min( ...
                    STA_d1_per_ms(STA_downstroke_start:STA_downstroke_end),[], 'omitnan');
            else
                STA_downstroke_value = NaN;
            end

            AP_STA_width_ms(iCell) = STA_width_frames ./ fps .* 1000;
            AP_STA_height(iCell) = STA_height_value;
            AP_STA_AHP(iCell) = STA_AHP_value;
            AP_STA_upstroke(iCell) = STA_upstroke_value;
            AP_STA_downstroke(iCell) = STA_downstroke_value;
        end
    end

    %% ====================================================================
    % 16. Layer-1 AP adaptation within each blue-on epoch
    % =====================================================================
    % x = AP peak time in milliseconds after the current stimulus starts.
    % Each AP feature uses its own finite-value mask. At least 3 valid APs
    % are required for a linear fit.
    for iStim = 1:nStim
        s = stim_start(iStim);
        eFrame = stim_end(iStim);
        event_idx = local_spikes >= s & local_spikes <= eFrame;

        if ~any(event_idx)
            continue
        end

        AP_time_ms = ...
            (double(local_spikes(event_idx)) - double(s)) ./ fps .* 1000;

        width20_ms = double(local_width20(event_idx)) ./ fps .* 1000;
        width50_ms = double(local_width50(event_idx)) ./ fps .* 1000;
        width80_ms = double(local_width80(event_idx)) ./ fps .* 1000;
        width20_rise_ms = double(local_width20_rise(event_idx)) ./ fps .* 1000;
        width20_fall_ms = double(local_width20_fall(event_idx)) ./ fps .* 1000;
        width50_rise_ms = double(local_width50_rise(event_idx)) ./ fps .* 1000;
        width50_fall_ms = double(local_width50_fall(event_idx)) ./ fps .* 1000;
        width80_rise_ms = double(local_width80_rise(event_idx)) ./ fps .* 1000;
        width80_fall_ms = double(local_width80_fall(event_idx)) ./ fps .* 1000;

        AHP_delay_ms = double(local_AHP_delay(event_idx));
        AHP_delay_ms(AHP_delay_ms < 0) = NaN;
        AHP_delay_ms = AHP_delay_ms ./ fps .* 1000;

        % Fit all 12 AP-adaptation features together. Each feature keeps its
        % own finite-value mask and still requires >= min_APs_for_adaptation_fit.
        AP_adaptation_Y = [ ...
            double(local_height(event_idx))', ...
            width20_ms(:), ...
            width50_ms(:), ...
            width80_ms(:), ...
            width20_rise_ms(:), ...
            width20_fall_ms(:), ...
            width50_rise_ms(:), ...
            width50_fall_ms(:), ...
            width80_rise_ms(:), ...
            width80_fall_ms(:), ...
            AHP_delay_ms(:), ...
            double(local_AHP(event_idx))'];

        [adapt_slope, adapt_intercept] = ...
            local_linear_fit_columns( ...
                AP_time_ms(:), ...
                AP_adaptation_Y, ...
                min_APs_for_adaptation_fit);

        AP_height_slope(iCell,iStim) = adapt_slope(1);
        AP_height_intercept(iCell,iStim) = adapt_intercept(1);

        AP_width20_slope(iCell,iStim) = adapt_slope(2);
        AP_width20_intercept(iCell,iStim) = adapt_intercept(2);

        AP_width50_slope(iCell,iStim) = adapt_slope(3);
        AP_width50_intercept(iCell,iStim) = adapt_intercept(3);

        AP_width80_slope(iCell,iStim) = adapt_slope(4);
        AP_width80_intercept(iCell,iStim) = adapt_intercept(4);

        AP_width20_rise_slope(iCell,iStim) = adapt_slope(5);
        AP_width20_rise_intercept(iCell,iStim) = adapt_intercept(5);

        AP_width20_fall_slope(iCell,iStim) = adapt_slope(6);
        AP_width20_fall_intercept(iCell,iStim) = adapt_intercept(6);

        AP_width50_rise_slope(iCell,iStim) = adapt_slope(7);
        AP_width50_rise_intercept(iCell,iStim) = adapt_intercept(7);

        AP_width50_fall_slope(iCell,iStim) = adapt_slope(8);
        AP_width50_fall_intercept(iCell,iStim) = adapt_intercept(8);

        AP_width80_rise_slope(iCell,iStim) = adapt_slope(9);
        AP_width80_rise_intercept(iCell,iStim) = adapt_intercept(9);

        AP_width80_fall_slope(iCell,iStim) = adapt_slope(10);
        AP_width80_fall_intercept(iCell,iStim) = adapt_intercept(10);

        AP_AHP_delay_slope(iCell,iStim) = adapt_slope(11);
        AP_AHP_delay_intercept(iCell,iStim) = adapt_intercept(11);

        AP_AHP_depth_slope(iCell,iStim) = adapt_slope(12);
        AP_AHP_depth_intercept(iCell,iStim) = adapt_intercept(12);
    end

    %% ====================================================================
    % 17. Whole-trace AP counts in consecutive fixed-width bins
    % =====================================================================
    AP_count_bin(iCell,:) = ...
        uint16(histcounts(double(local_spikes), count_bin_edges));

    %% ====================================================================
    % 18. Step subthreshold features
    % =====================================================================
    % step_height was calculated for all cells before the main loop.

    %% ====================================================================
    % 19. Ramp subthreshold/excitability features
    % =====================================================================
    for iRamp = 1:nRamps
        s = ramp_start(iRamp);
        eFrame = ramp_end(iRamp);
        spikes_epoch = local_spikes(local_spikes >= s & local_spikes <= eFrame);

        % ramp_area was calculated for all cells before the main loop.
        if ~isempty(spikes_epoch)
            rheobase_frame(iCell,iRamp) = spikes_epoch(1);
        end
    end

    %% ====================================================================
    % 20. Firing timing, frequency-vs-time adaptation, and maximum rate
    % =====================================================================
    for iStim = 1:nStim
        s = stim_start(iStim);
        eFrame = stim_end(iStim);
        spikes_epoch = local_spikes(local_spikes >= s & local_spikes <= eFrame);
        nEpochSpikes = numel(spikes_epoch);

        if nEpochSpikes >= 1
            first_spike_delay_ms(iCell,iStim) = ...
                (double(spikes_epoch(1)) - double(s)) ./ fps .* 1000;
            last_spike_delay_ms(iCell,iStim) = ...
                (double(spikes_epoch(end)) - double(s)) ./ fps .* 1000;
        end

        if nEpochSpikes >= 2
            isi_frames = diff(double(spikes_epoch));
            first_second_ISI_ms(iCell,iStim) = isi_frames(1) ./ fps .* 1000;
            penultimate_last_ISI_ms(iCell,iStim) = isi_frames(end) ./ fps .* 1000;
            max_spike_rate_Hz(iCell,iStim) = fps ./ min(isi_frames);
        end

        % Instantaneous frequency is assigned to the midpoint of each ISI.
        % Require >=4 spikes -> >=3 frequency/time points for the linear fit.
        if nEpochSpikes >= min_spikes_for_frequency_fit
            isi_frames = diff(double(spikes_epoch));
            inst_frequency_Hz = fps ./ isi_frames;
            isi_midpoint_frames = ...
                (double(spikes_epoch(1:end-1)) + double(spikes_epoch(2:end))) ./ 2;
            frequency_time_sec = ...
                (isi_midpoint_frames - double(s)) ./ fps;

            [frequency_time_slope_Hz_per_s(iCell,iStim), ...
             frequency_time_intercept_Hz(iCell,iStim)] = ...
                local_linear_fit(frequency_time_sec, inst_frequency_Hz, 3);
        end

        % Post-blue-off decay tau is fitted in one batch for all cells
        % before this cell loop (Section 12).
    end

% ========================================================================
% Progress report
% =========================================================================
    if analysis_config.progress.enable && (mod(iCell, progress_interval) == 0 || iCell == nCells)
    elapsed_time = toc(t_start);
    percent_done = 100 * iCell / nCells;
    estimated_total_time = elapsed_time / iCell * nCells;
    remaining_time = estimated_total_time - elapsed_time;
    fprintf( ...
        'Feature extraction: %5.1f%%  (%d/%d cells)  ', ...
        percent_done, iCell, nCells);
    fprintf( ...
        'Elapsed: %.1f min  Remaining: %.1f min\n', ...
        elapsed_time / 60, ...
        remaining_time / 60);
    end
end
%% ========================================================================
% 19. Overall feature QC
% =========================================================================
feature_qc.fail = feature_qc.invalid_input';

%% ========================================================================
% 20. Build Layer-1 feature TABLE
% =========================================================================
features = table((1:nCells)','VariableNames',{'cell_index'});

%% Whole-trace AP-count columns
for iBin = 1:nCountBins
    variable_name = sprintf('AP_count_bin_%03d',iBin);
    features.(variable_name) = AP_count_bin(:,iBin);
end

%% Overall AP phenotype
features.AP_n_spikes = AP_n_spikes;
features.AP_width20_median_ms = AP_width20_median_ms;
features.AP_width50_median_ms = AP_width50_median_ms;
features.AP_width80_median_ms = AP_width80_median_ms;
features.AP_width50_std_ms = AP_width50_std_ms;
features.AP_width50_CV = AP_width50_CV;
features.AP_width20_rise_median_ms = AP_width20_rise_median_ms;
features.AP_width20_fall_median_ms = AP_width20_fall_median_ms;
features.AP_width50_rise_median_ms = AP_width50_rise_median_ms;
features.AP_width50_fall_median_ms = AP_width50_fall_median_ms;
features.AP_width80_rise_median_ms = AP_width80_rise_median_ms;
features.AP_width80_fall_median_ms = AP_width80_fall_median_ms;
features.AP_height_median = AP_height_median;
features.AP_height_CV = AP_height_CV;
features.AP_AHP_median = AP_AHP_median;
features.AP_AHP_delay_median_ms = AP_AHP_delay_median_ms;
features.AP_upstroke_median = AP_upstroke_median;
features.AP_downstroke_median = AP_downstroke_median;

%% STA scalar AP phenotype (STA waveform is not stored)
features.AP_STA_width_ms = AP_STA_width_ms;
features.AP_STA_height = AP_STA_height;
features.AP_STA_AHP = AP_STA_AHP;
features.AP_STA_upstroke = AP_STA_upstroke;
features.AP_STA_downstroke = AP_STA_downstroke;

%% AP adaptation features: one slope/intercept pair per blue-on epoch
for iStim = 1:nStim
    e = stim_idx(iStim);
    epoch_name = matlab.lang.makeValidName(char(epochs(e).name));

    features.(sprintf('AP_height_slope_%s_per_ms',epoch_name)) = AP_height_slope(:,iStim);
    features.(sprintf('AP_height_intercept_%s',epoch_name)) = AP_height_intercept(:,iStim);

    features.(sprintf('AP_width20_slope_%s_ms_per_ms',epoch_name)) = AP_width20_slope(:,iStim);
    features.(sprintf('AP_width20_intercept_%s_ms',epoch_name)) = AP_width20_intercept(:,iStim);

    features.(sprintf('AP_width50_slope_%s_ms_per_ms',epoch_name)) = AP_width50_slope(:,iStim);
    features.(sprintf('AP_width50_intercept_%s_ms',epoch_name)) = AP_width50_intercept(:,iStim);

    features.(sprintf('AP_width80_slope_%s_ms_per_ms',epoch_name)) = AP_width80_slope(:,iStim);
    features.(sprintf('AP_width80_intercept_%s_ms',epoch_name)) = AP_width80_intercept(:,iStim);

    features.(sprintf('AP_width20_rise_slope_%s_ms_per_ms',epoch_name)) = AP_width20_rise_slope(:,iStim);
    features.(sprintf('AP_width20_rise_intercept_%s_ms',epoch_name)) = AP_width20_rise_intercept(:,iStim);

    features.(sprintf('AP_width20_fall_slope_%s_ms_per_ms',epoch_name)) = AP_width20_fall_slope(:,iStim);
    features.(sprintf('AP_width20_fall_intercept_%s_ms',epoch_name)) = AP_width20_fall_intercept(:,iStim);

    features.(sprintf('AP_width50_rise_slope_%s_ms_per_ms',epoch_name)) = AP_width50_rise_slope(:,iStim);
    features.(sprintf('AP_width50_rise_intercept_%s_ms',epoch_name)) = AP_width50_rise_intercept(:,iStim);

    features.(sprintf('AP_width50_fall_slope_%s_ms_per_ms',epoch_name)) = AP_width50_fall_slope(:,iStim);
    features.(sprintf('AP_width50_fall_intercept_%s_ms',epoch_name)) = AP_width50_fall_intercept(:,iStim);

    features.(sprintf('AP_width80_rise_slope_%s_ms_per_ms',epoch_name)) = AP_width80_rise_slope(:,iStim);
    features.(sprintf('AP_width80_rise_intercept_%s_ms',epoch_name)) = AP_width80_rise_intercept(:,iStim);

    features.(sprintf('AP_width80_fall_slope_%s_ms_per_ms',epoch_name)) = AP_width80_fall_slope(:,iStim);
    features.(sprintf('AP_width80_fall_intercept_%s_ms',epoch_name)) = AP_width80_fall_intercept(:,iStim);

    features.(sprintf('AP_AHP_delay_slope_%s_ms_per_ms',epoch_name)) = AP_AHP_delay_slope(:,iStim);
    features.(sprintf('AP_AHP_delay_intercept_%s_ms',epoch_name)) = AP_AHP_delay_intercept(:,iStim);

    features.(sprintf('AP_AHP_depth_slope_%s_per_ms',epoch_name)) = AP_AHP_depth_slope(:,iStim);
    features.(sprintf('AP_AHP_depth_intercept_%s',epoch_name)) = AP_AHP_depth_intercept(:,iStim);
end

%% Firing timing/adaptation features for each blue-on epoch
for iStim = 1:nStim
    e = stim_idx(iStim);
    epoch_name = matlab.lang.makeValidName(char(epochs(e).name));

    features.(sprintf('first_spike_delay_%s_ms',epoch_name)) = first_spike_delay_ms(:,iStim);
    features.(sprintf('first_second_ISI_%s_ms',epoch_name)) = first_second_ISI_ms(:,iStim);
    features.(sprintf('penultimate_last_ISI_%s_ms',epoch_name)) = penultimate_last_ISI_ms(:,iStim);
    features.(sprintf('last_spike_delay_%s_ms',epoch_name)) = last_spike_delay_ms(:,iStim);
    features.(sprintf('frequency_time_slope_%s_Hz_per_s',epoch_name)) = frequency_time_slope_Hz_per_s(:,iStim);
    features.(sprintf('frequency_time_intercept_%s_Hz',epoch_name)) = frequency_time_intercept_Hz(:,iStim);
    features.(sprintf('max_spike_rate_%s_Hz',epoch_name)) = max_spike_rate_Hz(:,iStim);
end

%% Step subthreshold features
for iStep = 1:nSteps
    e = step_idx(iStep);
    epoch_name = matlab.lang.makeValidName(char(epochs(e).name));
    features.(sprintf('step_height_%s',epoch_name)) = step_height(:,iStep);
end

%% Ramp subthreshold/excitability features
for iRamp = 1:nRamps
    e = ramp_idx(iRamp);
    epoch_name = matlab.lang.makeValidName(char(epochs(e).name));
    features.(sprintf('ramp_area_%s',epoch_name)) = ramp_area(:,iRamp);
    features.(sprintf('rheobase_frame_%s',epoch_name)) = rheobase_frame(:,iRamp);
end

%% Post-blue-off exponential recovery features
for iStim = 1:nStim
    e = stim_idx(iStim);
    epoch_name = matlab.lang.makeValidName(char(epochs(e).name));

    features.(sprintf('decay_tau_%s_ms',epoch_name)) = decay_tau_ms(:,iStim);
end

%% ============================================================
% Final feature extraction report
% =============================================================
total_time = toc(t_start);

fprintf('\n');
fprintf('============================================================\n');
fprintf('Feature extraction complete\n');
fprintf('============================================================\n');

fprintf('Cells analyzed: %d\n', nCells);
fprintf('Total time:     %.2f min\n', total_time / 60);
fprintf('Time per cell:  %.3f sec\n', total_time / nCells);
end



%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function [width20, rise20, fall20, ...
          width50, rise50, fall50, ...
          width80, rise80, fall80] = ...
    local_AP_widths_20_50_80( ...
        trace, onset_idx, peak_idx, onset_value, pulse_height)
% LOCAL_AP_WIDTHS_20_50_80
% Compute 20%, 50%, and 80% AP widths from one short waveform.
% The crossing helper is shared so Layer-0 keeps identical width semantics
% while avoiding repeated setup in the calling code.

level20 = onset_value + 0.20 * pulse_height;
level50 = onset_value + 0.50 * pulse_height;
level80 = onset_value + 0.80 * pulse_height;

[width20, rise20, fall20] = ...
    local_width_at_level(trace, onset_idx, peak_idx, level20);

[width50, rise50, fall50] = ...
    local_width_at_level(trace, onset_idx, peak_idx, level50);

[width80, rise80, fall80] = ...
    local_width_at_level(trace, onset_idx, peak_idx, level80);
end


function [width_frames, rise_to_peak_frames, peak_to_fall_frames] = ...
    local_width_at_level(trace, onset_idx, peak_idx, threshold)
% LOCAL_WIDTH_AT_LEVEL
% Estimate AP width at a specified signal threshold using linearly
% interpolated rising- and falling-edge crossings.
%
% Outputs are in frames and may be fractional:
%   width_frames         = falling crossing - rising crossing
%   rise_to_peak_frames  = peak index - rising crossing
%   peak_to_fall_frames  = falling crossing - peak index
%
% Returns NaN when either crossing cannot be identified.

width_frames = NaN;
rise_to_peak_frames = NaN;
peak_to_fall_frames = NaN;

if onset_idx >= peak_idx
    return
end

%% Rising-edge crossing
left_trace = trace(onset_idx:peak_idx);
left_idx = onset_idx:peak_idx;

left_crossing = ...
    (left_trace(1:end-1)-threshold) .* ...
    (left_trace(2:end)-threshold) <= 0;

% Use the last crossing before the peak to reduce sensitivity to small
% pre-spike fluctuations.
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

left_position = ...
    x1 + (threshold-y1) / (y2-y1);

%% Falling-edge crossing
right_trace = trace(peak_idx:end);
right_idx = peak_idx:numel(trace);

right_crossing = ...
    (right_trace(1:end-1)-threshold) .* ...
    (right_trace(2:end)-threshold) <= 0;

% First crossing after the AP peak.
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

right_position = ...
    x1 + (threshold-y1) / (y2-y1);

%% Width components
rise_to_peak_frames = peak_idx - left_position;
peak_to_fall_frames = right_position - peak_idx;
width_frames = right_position - left_position;
end



function [slope, intercept] = local_linear_fit_columns(x,Y,min_points)
% LOCAL_LINEAR_FIT_COLUMNS
% Fast OLS for multiple y features sharing the same x:
%
%       Y(:,j) = slope(j)*x + intercept(j)
%
% Each Y column has its own finite-value mask.

x = double(x(:));
Y = double(Y);

if isempty(Y)
    slope = nan(1,0);
    intercept = nan(1,0);
    return
end

if size(Y,1) ~= numel(x)
    error('local_linear_fit_columns: x and Y have incompatible sizes.');
end

nFeat = size(Y,2);

valid_x = isfinite(x);
valid = isfinite(Y) & valid_x;

n = sum(valid,1);

x_safe = x;
x_safe(~valid_x) = 0;

x_masked = x_safe .* valid;

Y_safe = Y;
Y_safe(~valid) = 0;

sum_x = sum(x_masked,1);
sum_y = sum(Y_safe,1);
sum_xx = sum(x_masked.^2,1);
sum_xy = sum(x_masked .* Y_safe,1);

n_safe = max(n,1);

denom = ...
    sum_xx - (sum_x.^2) ./ n_safe;

numer = ...
    sum_xy - (sum_x .* sum_y) ./ n_safe;

slope = numer ./ denom;
intercept = (sum_y - slope .* sum_x) ./ n_safe;

mean_x = sum_x ./ n_safe;

bad = ...
    n < min_points | ...
    ~isfinite(denom) | ...
    denom <= eps(max(abs(mean_x),1));

slope(bad) = NaN;
intercept(bad) = NaN;

slope = reshape(slope,1,nFeat);
intercept = reshape(intercept,1,nFeat);
end


function [tau_ms, fit_ok] = ...
    local_fast_decay_tau_batch(Y,t_sec,valid_rows,tau_grid_ms)
% LOCAL_FAST_DECAY_TAU_BATCH
% Fast batch fitting of:
%
%       y(t) = a*exp(-t/tau) + c
%
% Only tau is returned. For each candidate tau, a and c are solved by
% ordinary least squares. Cells with fully finite decay windows are fitted
% together using matrix multiplication. Rare rows containing NaNs fall back
% to the same grid solver one row at a time.
%
% The best log-grid point is refined by 3-point parabolic interpolation in
% log(tau), so output precision is finer than the discrete grid spacing.

nRows = size(Y,1);
tau_ms = nan(nRows,1);
fit_ok = false(nRows,1);

t_sec = double(t_sec(:));
tau_grid_ms = double(tau_grid_ms(:)');

if size(Y,2) ~= numel(t_sec)
    error('Decay trace width and decay time vector do not match.');
end

if numel(t_sec) < 5 || numel(tau_grid_ms) < 3
    return
end

rows = find(valid_rows(:));

if isempty(rows)
    return
end

tau_grid_sec = tau_grid_ms ./ 1000;

% Basis matrix: time x tau candidate.
E = exp(-t_sec ./ tau_grid_sec);

% Center each basis because c is a free intercept.
E_centered = E - mean(E,1);
E_ss = sum(E_centered.^2,1);

basis_ok = isfinite(E_ss) & E_ss > eps;
if nnz(basis_ok) < 3
    return
end

E_centered = E_centered(:,basis_ok);
E_ss = E_ss(basis_ok);
tau_use_ms = tau_grid_ms(basis_ok);

finite_all = all(isfinite(Y(rows,:)),2);
rows_full = rows(finite_all);
rows_partial = rows(~finite_all);

%% Fast matrix path for fully finite rows
if ~isempty(rows_full)
    Y_full = double(Y(rows_full,:));
    Y_centered = Y_full - mean(Y_full,2);

    Y_ss = sum(Y_centered.^2,2);

    % Explained SS for each tau candidate:
    %       covariance(y,e)^2 / SS(e)
    % Maximizing this is exactly equivalent to minimizing least-squares SSE
    % over amplitude a and offset c for that fixed tau.
    cross_term = Y_centered * E_centered;
    score = (cross_term.^2) ./ E_ss;

    [best_score,best_idx] = max(score,[],2);

    good = ...
        isfinite(Y_ss) & ...
        Y_ss > eps & ...
        isfinite(best_score) & ...
        best_score > 0;

    if any(good)
        good_rows_local = find(good);
        idx_good = best_idx(good);
        best_tau = tau_use_ms(idx_good);
        best_tau = best_tau(:);

        % Parabolic refinement on the uniformly spaced log(tau) grid.
        log_tau = log(tau_use_ms);
        interior = idx_good > 1 & idx_good < numel(tau_use_ms);

        if any(interior)
            rg = good_rows_local(interior);
            j = idx_good(interior);

            lin1 = sub2ind(size(score), rg, j-1);
            lin2 = sub2ind(size(score), rg, j);
            lin3 = sub2ind(size(score), rg, j+1);

            s1 = score(lin1);
            s2 = score(lin2);
            s3 = score(lin3);

            denom_q = s1 - 2.*s2 + s3;
            delta = zeros(size(s2));

            valid_q = isfinite(denom_q) & abs(denom_q) > eps;
            delta(valid_q) = ...
                0.5 .* (s1(valid_q) - s3(valid_q)) ./ denom_q(valid_q);

            delta = max(-1,min(1,delta));

            dlog = log_tau(2) - log_tau(1);
            log_tau_j = reshape(log_tau(j),[],1);
            refined = exp(log_tau_j + delta .* dlog);

            best_tau(interior) = refined;
        end

        tau_ms(rows_full(good)) = best_tau;
        fit_ok(rows_full(good)) = true;
    end
end

%% Rare fallback path for rows containing NaNs
for k = 1:numel(rows_partial)
    r = rows_partial(k);

    y = double(Y(r,:))';
    valid = isfinite(y) & isfinite(t_sec);

    if nnz(valid) < 5
        continue
    end

    [tau_one,ok_one] = ...
        local_fast_decay_tau_single( ...
            y(valid), ...
            t_sec(valid), ...
            tau_use_ms);

    if ok_one
        tau_ms(r) = tau_one;
        fit_ok(r) = true;
    end
end
end


function [tau_ms,ok] = ...
    local_fast_decay_tau_single(y,t_sec,tau_grid_ms)
% Same offset-exponential tau grid fit as the batch helper, for one row
% with a reduced finite time mask.

tau_ms = NaN;
ok = false;

y = double(y(:));
t_sec = double(t_sec(:));
tau_grid_ms = double(tau_grid_ms(:)');

if numel(y) < 5 || numel(tau_grid_ms) < 3
    return
end

y_centered = y - mean(y);
y_ss = sum(y_centered.^2);

if ~isfinite(y_ss) || y_ss <= eps
    return
end

E = exp(-t_sec ./ (tau_grid_ms./1000));
Ec = E - mean(E,1);
Ess = sum(Ec.^2,1);

good_basis = isfinite(Ess) & Ess > eps;
if nnz(good_basis) < 3
    return
end

Ec = Ec(:,good_basis);
Ess = Ess(good_basis);
tau_use = tau_grid_ms(good_basis);

cross_term = y_centered' * Ec;
score = (cross_term.^2) ./ Ess;

[best_score,j] = max(score);

if ~isfinite(best_score) || best_score <= 0
    return
end

tau_ms = tau_use(j);

if j > 1 && j < numel(tau_use)
    s1 = score(j-1);
    s2 = score(j);
    s3 = score(j+1);
    denom_q = s1 - 2*s2 + s3;

    if isfinite(denom_q) && abs(denom_q) > eps
        delta = 0.5 * (s1-s3) / denom_q;
        delta = max(-1,min(1,delta));

        log_tau = log(tau_use);
        dlog = log_tau(2)-log_tau(1);
        tau_ms = exp(log_tau(j) + delta*dlog);
    end
end

ok = isfinite(tau_ms) && tau_ms > 0;
end


function [slope, intercept] = local_linear_fit(x,y,min_points)
% LOCAL_LINEAR_FIT Fast ordinary least-squares fit y = slope*x + intercept.
% Uses the closed-form OLS solution to avoid repeated POLYFIT overhead.
% Returns NaN/NaN when fewer than min_points valid pairs are available or
% when all valid x values are identical.

x = double(x(:));
y = double(y(:));
valid = isfinite(x) & isfinite(y);

if nnz(valid) < min_points
    slope = NaN;
    intercept = NaN;
    return
end

x = x(valid);
y = y(valid);

x_mean = mean(x);
y_mean = mean(y);
dx = x - x_mean;
denom = sum(dx.^2);

if ~isfinite(denom) || denom <= eps(max(abs(x_mean),1))
    slope = NaN;
    intercept = NaN;
    return
end

slope = sum(dx .* (y - y_mean)) ./ denom;
intercept = y_mean - slope .* x_mean;
end


function value = local_nan_std(x)
% LOCAL_NAN_STD Sample standard deviation using finite values only.
x = double(x(:));
x = x(isfinite(x));
if numel(x) < 2
    value = NaN;
else
    value = std(x);
end
end


function value = local_nan_cv(x)
% LOCAL_NAN_CV Positive coefficient of variation: std / abs(mean).
x = double(x(:));
x = x(isfinite(x));
if numel(x) < 2
    value = NaN;
    return
end
mu = mean(x);
if ~isfinite(mu) || abs(mu) < eps
    value = NaN;
else
    value = std(x) ./ abs(mu);
end
end
