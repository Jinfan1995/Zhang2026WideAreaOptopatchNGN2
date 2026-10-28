function analysis_config = make_trace_analysis_config_v2()
% MAKE_TRACE_ANALYSIS_CONFIG_V2
% Define all analysis parameters used by the Optopatch
% Trace Analysis V2 pipeline.
% This structure describes HOW the data should be analyzed.
% Experimental information such as:
%   - fps
%   - stimulus timing
%   - baseline frames
%   - ramp/step epochs
%
% Usage:
%   analysis_config = make_trace_analysis_config_v2();
% Output:
%   analysis_config
%       MATLAB structure containing analysis/QC parameters.
%
% ==========================================================
% Design principle
% ==========================================================
% protocol
%       = what experiment was performed
% analysis_config
%       = how the data should be analyzed
% ==========================================================


%% =========================================================
% 1. Configuration identity
% ==========================================================
analysis_config = struct;
analysis_config.name = "Optopatch_trace_analysis";
analysis_config.version = "TraceAnalysis_v2.4_residual_drift_QC";
analysis_config.created_date = datetime('today');

%% =========================================================
% 2. Exponential drift correction
% ==========================================================
% Primary bleaching/drift model.
analysis_config.exp_fit.model = "exp2";

% If exp2 becomes numerically/physically unstable, try a simpler exp1 fit
% before declaring the cell an exponential-fit failure.
analysis_config.exp_fit.use_exp1_fallback = true;

% Reject a fitted baseline if its absolute value anywhere in the movie falls
% below this fraction of the robust observed baseline fluorescence.
% This prevents division by a small denominator from generating huge spikes.
analysis_config.exp_fit.min_fitted_baseline_fraction = 0.05;

% After correction, reject the trace if the largest absolute corrected value
% is more than this fold above the typical corrected baseline magnitude.
% Corrected fluorescence should normally stay close to order 1.
analysis_config.exp_fit.max_corrected_fold = 50;

%% =========================================================
% 3. Residual baseline-drift QC and correction
% ==========================================================
%
% This module runs AFTER exponential correction and normalization.
%
% For each baseline period:
%   1. calculate the median normalized trace value
%   2. calculate baseline_drift_range =
%          max(period medians) - min(period medians)
%
% Traces with range >= max_baseline_range fail QC and are not corrected.
% Traces below the threshold are corrected by subtracting a linear trend
% fitted only to the baseline-period medians.

analysis_config.drift.enable = true;

% User-defined maximum acceptable normalized baseline drift.
% "Meaningful" traces must have baseline drift strictly < 0.15.
analysis_config.drift.max_baseline_range = 0.15;

% Minimum number of valid baseline-period anchors needed for linear fitting.
analysis_config.drift.min_anchor_periods = 3;

% Current correction model.
analysis_config.drift.method = "linear";

% Failed drift traces are returned as NaN so they cannot accidentally enter
% STA/kernel/feature extraction.
analysis_config.drift.nan_failed_traces = true;

%% =========================================================
% 4. SNR analysis / QC
% ==========================================================
% Minimum step-response SNR
analysis_config.snr.step_threshold = 4;

% Minimum spike SNR
analysis_config.snr.spike_threshold = 8;

% ---------------------------------------------------------
% Extreme / artifact thresholds
% ---------------------------------------------------------
% Extremely large step SNR can indicate an abnormal/artifact trace.
analysis_config.snr.extreme_step_threshold = 100;

% Extremely large spike SNR can indicate an abnormal/artifact trace.
analysis_config.snr.extreme_spike_threshold = 100;

% Extremely large high-frequency baseline noise.
% noise_std is calculated from the high-pass trace over baseline frames.
analysis_config.snr.extreme_noise_std_threshold = 1;

% Exclude frames close to beginning/end when calculating
% whole-trace SNR metrics
analysis_config.snr.edge_exclusion_frames = 400;

%% =========================================================
% 5. STA kernel generation
% ==========================================================
% Initial spike threshold used to generate STA kernel
analysis_config.sta.spike_SNR_threshold = 8;
% Minimum separation used by findpeaks
analysis_config.sta.min_peak_distance_frames = 2;
% Minimum allowed interval between accepted spikes
analysis_config.sta.min_spike_interval_frames = 10;
% ---------------------------------------------------------
% Trace-edge exclusion
% ---------------------------------------------------------
analysis_config.sta.edge_start_exclusion_frames = 200;
analysis_config.sta.edge_end_exclusion_frames = 31;
% ---------------------------------------------------------
% Initial AP extraction around detected spike
% ---------------------------------------------------------
analysis_config.sta.extract_pre_frames = 15;
analysis_config.sta.extract_post_frames = 30;
% ---------------------------------------------------------
% Spike-onset localization
% ---------------------------------------------------------
% Light smoothing used ONLY for derivative calculation
analysis_config.onset.smooth_window_frames = 3;
% Restrict onset search to ONLY the 3 frames immediately before the AP peak
% Allowed onset indices relative to the peak: peak-3, peak-2, peak-1
analysis_config.onset.search_prepeak_frames = 3;
% Define beginning of upstroke as fraction of max dV/dt
analysis_config.onset.dvdt_fraction = 0.15;
% dV/dt must remain above threshold for this many frames
analysis_config.onset.min_sustain_frames = 2;
% After detecting the rising phase, refine onset using d2V/dt2
analysis_config.onset.d2_refine_radius_frames = 1;
% ---------------------------------------------------------
% STA kernel
% ---------------------------------------------------------
% STA kernel construction is now PEAK-aligned. Onset is detected once on
% the averaged STA for diagnostics, not used to align individual spikes.
analysis_config.sta.kernel_half_width_frames = 15;

%% =========================================================
% 6. Final feature extraction
% ==========================================================
% Spike threshold after convolution with the STA kernel
analysis_config.features.spike_SNR_threshold = 6;
% Minimum separation between peaks
analysis_config.features.min_peak_distance_frames = 5;
% ---------------------------------------------------------
% Refine detected peak
% ---------------------------------------------------------
% Search +/- this many frames around the preliminary peak
analysis_config.features.peak_refine_radius_frames = 2;
% ---------------------------------------------------------
% Edge exclusion
% ---------------------------------------------------------
analysis_config.features.spike_edge_pre_frames = 50;
analysis_config.features.spike_edge_post_frames = 101;
% ---------------------------------------------------------
% AP waveform extraction
% ---------------------------------------------------------
analysis_config.features.spike_extract_pre_frames = 20;
analysis_config.features.spike_extract_post_frames = 30;
% ---------------------------------------------------------
% Spike-onset localization
% ---------------------------------------------------------
% Uses the shared analysis_config.onset settings above.
% No fixed 13:15 second-derivative search window is used.
% ---------------------------------------------------------
% Onset-aligned AP waveform
% ---------------------------------------------------------
analysis_config.features.onset_pre_frames = 10;
analysis_config.features.onset_post_frames = 30;
% ---------------------------------------------------------
% Mixed-trace QC
% ---------------------------------------------------------
analysis_config.features.close_spike_fail_frames = 5;
% ---------------------------------------------------------
% Post-stimulus decay fitting
% ---------------------------------------------------------
analysis_config.features.decay_trim_start_frames = 10;
analysis_config.features.decay_trim_end_frames = 30;

%% =========================================================
% 7. Progress reporting
% ==========================================================
analysis_config.progress.enable = true;
analysis_config.progress.percent_interval = 5;

%% =========================================================
% 8. Data storage
% ==========================================================
analysis_config.storage.save_processed_traces = true;
analysis_config.storage.save_exp_corrected_traces = false;
analysis_config.storage.save_normalized_traces = true;



%% =========================================================
% 9. Validate configuration
% ==========================================================
assert(analysis_config.exp_fit.min_fitted_baseline_fraction > 0 && ...
    analysis_config.exp_fit.min_fitted_baseline_fraction < 1, ...
    'exp_fit min_fitted_baseline_fraction must be between 0 and 1.');
assert(analysis_config.exp_fit.max_corrected_fold > 1, ...
    'exp_fit max_corrected_fold must be > 1.');
assert(analysis_config.drift.max_baseline_range > 0, ...
    'Residual drift max_baseline_range must be > 0.');
assert(analysis_config.drift.min_anchor_periods >= 2, ...
    'Residual drift min_anchor_periods must be >= 2.');
assert(ismember(string(analysis_config.drift.method), "linear"), ...
    'Residual drift method must currently be "linear".');
assert(analysis_config.snr.step_threshold > 0, ...
    'SNR step threshold must be > 0.');
assert(analysis_config.snr.spike_threshold > 0, ...
    'SNR spike threshold must be > 0.');
assert(analysis_config.snr.extreme_step_threshold > 0, ...
    'Extreme step SNR threshold must be > 0.');
assert(analysis_config.snr.extreme_spike_threshold > 0, ...
    'Extreme spike SNR threshold must be > 0.');
assert(analysis_config.snr.extreme_noise_std_threshold > 0, ...
    'Extreme noise std threshold must be > 0.');
assert(analysis_config.snr.edge_exclusion_frames >= 0, ...
    'SNR edge exclusion must be >= 0 frames.');
assert(analysis_config.sta.spike_SNR_threshold > 0, ...
    'STA spike SNR threshold must be > 0.');
assert(analysis_config.features.spike_SNR_threshold > 0, ...
    'Feature spike SNR threshold must be > 0.');
assert(analysis_config.progress.percent_interval > 0 && ...
    analysis_config.progress.percent_interval <= 100, ...
    'Progress interval must be between 0 and 100.');
assert(analysis_config.features.spike_extract_pre_frames >= ...
    analysis_config.features.onset_pre_frames, ...
    'AP extraction window is too short before the spike.');
assert(analysis_config.features.spike_extract_post_frames >= ...
    analysis_config.features.onset_post_frames, ...
    'AP extraction window is too short after the spike.');
assert(analysis_config.features.spike_extract_pre_frames >= ...
    analysis_config.onset.search_prepeak_frames + ...
    analysis_config.features.onset_pre_frames, ...
    ['AP extraction window is too short for the requested onset search ' ...
     'plus pre-onset waveform window.']);
assert(analysis_config.sta.extract_pre_frames >= ...
    analysis_config.sta.kernel_half_width_frames, ...
    'STA extraction needs enough pre-peak frames for the kernel.');
assert(analysis_config.sta.extract_post_frames >= ...
    analysis_config.sta.kernel_half_width_frames, ...
    'STA extraction needs enough post-peak frames for the kernel.');
assert(analysis_config.onset.smooth_window_frames >= 1, ...
    'Onset smoothing window must be >= 1 frame.');
assert(analysis_config.onset.dvdt_fraction > 0 && ...
    analysis_config.onset.dvdt_fraction < 1, ...
    'Onset dV/dt fraction must be between 0 and 1.');
end