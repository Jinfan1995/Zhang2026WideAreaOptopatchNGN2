function snr_info = traces_SNR_filter( ...
    traces_nm, ...
    baseline_frames, ...
    median_filter_size, ...
    step_SNR_th, ...
    spike_SNR_th, ...
    snr_config)
% TRACES_SNR_FILTER
% Calculate step SNR, spike SNR, baseline noise, and SNR-related QC for
% normalized Optopatch traces.
%
% INPUT
% -------------------------------------------------------------------------
% traces_nm
%       nCells x nFrames normalized traces.
%
% baseline_frames
%       Frames used to estimate high-frequency baseline noise.
%
% median_filter_size
%       Window used for low-pass filtering.
%
% step_SNR_th
%       Minimum acceptable step-response SNR.
%
% spike_SNR_th
%       Minimum acceptable spike SNR.
%
% snr_config   (optional)
%       analysis_config.snr structure. Supported fields:
%
%       extreme_step_threshold
%       extreme_spike_threshold
%       extreme_noise_std_threshold
%       edge_exclusion_frames
%
%       If snr_config is omitted, defaults are:
%
%       extreme_step_threshold      = 100
%       extreme_spike_threshold     = 100
%       extreme_noise_std_threshold = 1
%       edge_exclusion_frames       = 400
%
% OUTPUT
% -------------------------------------------------------------------------
% snr_info.step
% snr_info.spike
% snr_info.noise_std
%
% Detailed QC flags:
%
% snr_info.fail_low_step
% snr_info.fail_low_spike
% snr_info.fail_step_gt_spike
% snr_info.fail_extreme_step
% snr_info.fail_extreme_spike
% snr_info.fail_extreme_noise
% snr_info.fail_invalid
%
% snr_info.fail
%       Combined SNR-QC failure flag. This still includes low-SNR cells
%       for reporting purposes.
%
% IMPORTANT PIPELINE CONVENTION
% -------------------------------------------------------------------------
% For Layer-2 exclusion, it is recommended to build qc_abnormal_SNR from:
%
%   fail_step_gt_spike
%   fail_extreme_step
%   fail_extreme_spike
%   fail_extreme_noise
%   fail_invalid
%
% Low SNR alone does not have to exclude a biologically silent cell from
% Layer 2.
% -------------------------------------------------------------------------

%% ========================================================================
% 1. Input information
% =========================================================================

[nCells,nFrames] = size(traces_nm);

baseline_frames = ...
    unique(baseline_frames(:));

baseline_frames = ...
    baseline_frames( ...
        baseline_frames >= 1 & ...
        baseline_frames <= nFrames);

if isempty(baseline_frames)
    error( ...
        'traces_SNR_filter:NoBaselineFrames', ...
        'No valid baseline frames were supplied.');
end

if median_filter_size < 1
    error( ...
        'traces_SNR_filter:InvalidMedianFilterSize', ...
        'median_filter_size must be >= 1.');
end


%% ========================================================================
% 2. Read optional SNR configuration
% =========================================================================

extreme_step_SNR_th = 100;
extreme_spike_SNR_th = 100;
extreme_noise_std_th = 1;
edge_exclusion = 400;

if nargin >= 6 && ~isempty(snr_config)

    if ~isstruct(snr_config)
        error( ...
            'traces_SNR_filter:InvalidSNRConfig', ...
            'snr_config must be a structure.');
    end

    if isfield(snr_config,'extreme_step_threshold')
        extreme_step_SNR_th = ...
            double(snr_config.extreme_step_threshold);
    end

    if isfield(snr_config,'extreme_spike_threshold')
        extreme_spike_SNR_th = ...
            double(snr_config.extreme_spike_threshold);
    end

    if isfield(snr_config,'extreme_noise_std_threshold')
        extreme_noise_std_th = ...
            double(snr_config.extreme_noise_std_threshold);
    end

    if isfield(snr_config,'edge_exclusion_frames')
        edge_exclusion = ...
            double(snr_config.edge_exclusion_frames);
    end

end

if ~isfinite(extreme_step_SNR_th) || extreme_step_SNR_th <= 0
    error( ...
        'traces_SNR_filter:InvalidExtremeStepThreshold', ...
        'Extreme step SNR threshold must be > 0.');
end

if ~isfinite(extreme_spike_SNR_th) || extreme_spike_SNR_th <= 0
    error( ...
        'traces_SNR_filter:InvalidExtremeSpikeThreshold', ...
        'Extreme spike SNR threshold must be > 0.');
end

if ~isfinite(extreme_noise_std_th) || extreme_noise_std_th <= 0
    error( ...
        'traces_SNR_filter:InvalidExtremeNoiseThreshold', ...
        'Extreme noise std threshold must be > 0.');
end

if ~isfinite(edge_exclusion) || edge_exclusion < 0
    error( ...
        'traces_SNR_filter:InvalidEdgeExclusion', ...
        'edge_exclusion_frames must be >= 0.');
end

edge_exclusion = floor(edge_exclusion);


%% ========================================================================
% 3. Low-pass and high-pass traces
% =========================================================================

lowpass_traces = ...
    movmedian( ...
        traces_nm, ...
        median_filter_size, ...
        2, ...
        'omitnan');

highpass_traces = ...
    traces_nm - lowpass_traces;


%% ========================================================================
% 4. Estimate baseline noise
% =========================================================================

noise_std = ...
    std( ...
        highpass_traces(:,baseline_frames), ...
        0, ...
        2, ...
        'omitnan');


%% ========================================================================
% 5. Define frames used for SNR signal measurement
% =========================================================================

signal_start = ...
    1 + edge_exclusion;

signal_end = ...
    nFrames - edge_exclusion;

if signal_start > signal_end
    error( ...
        'traces_SNR_filter:MovieTooShort', ...
        ['Movie is too short for a %d-frame edge exclusion ' ...
         'at both ends.'], ...
        edge_exclusion);
end

signal_frames = ...
    signal_start:signal_end;


%% ========================================================================
% 6. Calculate SNR
% =========================================================================

step_signal = ...
    max( ...
        lowpass_traces(:,signal_frames), ...
        [], ...
        2, ...
        'omitnan');

spike_signal = ...
    max( ...
        highpass_traces(:,signal_frames), ...
        [], ...
        2, ...
        'omitnan');

SNR_step = ...
    step_signal ./ noise_std;

SNR_spike = ...
    spike_signal ./ noise_std;


%% ========================================================================
% 7. Identify invalid traces
% =========================================================================

invalid_trace = ...
    ~any(isfinite(traces_nm),2);

invalid_noise = ...
    ~isfinite(noise_std) | ...
    noise_std <= 0;

invalid_SNR = ...
    ~isfinite(SNR_step) | ...
    ~isfinite(SNR_spike);


%% ========================================================================
% 8. Apply SNR QC criteria
% =========================================================================

low_step_SNR = ...
    SNR_step < step_SNR_th;

low_spike_SNR = ...
    SNR_spike < spike_SNR_th;

step_larger_than_spike = ...
    SNR_step > 1.5 .* SNR_spike;

extreme_step_SNR = ...
    SNR_step > extreme_step_SNR_th;

extreme_spike_SNR = ...
    SNR_spike > extreme_spike_SNR_th;

extreme_noise_std = ...
    noise_std > extreme_noise_std_th;


% Broad SNR-QC report flag.
% Low SNR remains included here for reporting, but your final qc_pass can
% exclude only the abnormal/artifact criteria.
fail = ...
    invalid_trace | ...
    invalid_noise | ...
    invalid_SNR | ...
    low_step_SNR | ...
    low_spike_SNR | ...
    step_larger_than_spike | ...
    extreme_step_SNR | ...
    extreme_spike_SNR | ...
    extreme_noise_std;


%% ========================================================================
% 9. Store measurements
% =========================================================================

snr_info = struct;

snr_info.fail = ...
    logical(fail);

snr_info.step = ...
    SNR_step;

snr_info.spike = ...
    SNR_spike;

snr_info.noise_std = ...
    noise_std;


%% ========================================================================
% 10. Store detailed QC reasons
% =========================================================================

snr_info.fail_low_step = ...
    logical(low_step_SNR);

snr_info.fail_low_spike = ...
    logical(low_spike_SNR);

snr_info.fail_step_gt_spike = ...
    logical(step_larger_than_spike);

snr_info.fail_extreme_step = ...
    logical(extreme_step_SNR);

snr_info.fail_extreme_spike = ...
    logical(extreme_spike_SNR);

snr_info.fail_extreme_noise = ...
    logical(extreme_noise_std);

snr_info.fail_invalid = ...
    logical( ...
        invalid_trace | ...
        invalid_noise | ...
        invalid_SNR);


%% ========================================================================
% 11. Save thresholds used for reproducibility
% =========================================================================

snr_info.thresholds = struct;

snr_info.thresholds.step_min = ...
    double(step_SNR_th);

snr_info.thresholds.spike_min = ...
    double(spike_SNR_th);

snr_info.thresholds.extreme_step = ...
    extreme_step_SNR_th;

snr_info.thresholds.extreme_spike = ...
    extreme_spike_SNR_th;

snr_info.thresholds.extreme_noise_std = ...
    extreme_noise_std_th;

snr_info.thresholds.edge_exclusion_frames = ...
    edge_exclusion;


%% ========================================================================
% 12. Final report
% =========================================================================

nFailed = ...
    nnz(snr_info.fail);

fprintf( ...
    'SNR filter failed: %d / %d cells (%.2f%%)\n', ...
    nFailed, ...
    nCells, ...
    100 .* nFailed ./ nCells);

fprintf( ...
    ['Extreme SNR/noise QC: step=%d, spike=%d, noise=%d cells\n'], ...
    nnz(snr_info.fail_extreme_step), ...
    nnz(snr_info.fail_extreme_spike), ...
    nnz(snr_info.fail_extreme_noise));

end
