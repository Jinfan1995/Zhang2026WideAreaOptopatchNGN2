function [traces_corrected, drift_info] = ...
    traces_residual_drift_corr( ...
        traces_nm, ...
        protocol, ...
        drift_config)
% TRACES_RESIDUAL_DRIFT_CORR
%
% Quantify and correct residual baseline drift AFTER normalization.
%
% Biological / QC definition used by this pipeline:
%
%   baseline_drift_range =
%       max(median value of each baseline period) -
%       min(median value of each baseline period)
%
% A trace is considered usable only when:
%
%       baseline_drift_range < drift_config.max_baseline_range
%
% Current default threshold:
%
%       max_baseline_range = 0.15
%
% For traces that pass this QC, fit:
%
%       baseline(t) = slope * t_sec + intercept
%
% using ONLY the baseline-period medians, then subtract that line from the
% entire normalized trace:
%
%       corrected_trace = traces_nm - fitted_baseline
%
% Traces that fail QC are returned as NaN when:
%
%       drift_config.nan_failed_traces = true
%
% INPUT
% -------------------------------------------------------------------------
% traces_nm
%       nCells x nFrames normalized traces.
%
% protocol
%       Must contain protocol.fps and baseline information.
%
%       Preferred:
%           protocol.baseline.period(i).frames
%
%       Fallback:
%           protocol.baseline.frames
%       In fallback mode, contiguous baseline-frame blocks are automatically
%       converted into baseline periods.
%
% drift_config
%       analysis_config.drift
%
% OUTPUT
% -------------------------------------------------------------------------
% traces_corrected
%       Residual-drift-corrected normalized traces.
%
% drift_info.anchor_frames
%       1 x nPeriods anchor frame for each baseline period.
%
% drift_info.anchor_values
%       nCells x nPeriods median normalized baseline values.
%
% drift_info.baseline_range
%       Pre-correction max(anchor)-min(anchor).
%
% drift_info.slope_per_sec
% drift_info.intercept
% drift_info.r2
% drift_info.post_correction_range
% drift_info.n_valid_anchors
%
% drift_info.fail
% drift_info.correction_applied
% drift_info.fail_reason
%
% -------------------------------------------------------------------------

%% ========================================================================
% 1. Inputs and configuration
% =========================================================================

[nCells,nFrames] = ...
    size(traces_nm);

if ~isfield(protocol,'fps') || ...
        ~isfinite(protocol.fps) || ...
        protocol.fps <= 0

    error( ...
        'traces_residual_drift_corr:InvalidFPS', ...
        'protocol.fps must be finite and > 0.');

end

max_baseline_range = 0.15;
min_anchor_periods = 3;
method = "linear";
nan_failed_traces = true;

if nargin >= 3 && ~isempty(drift_config)

    if ~isstruct(drift_config)
        error( ...
            'traces_residual_drift_corr:InvalidConfig', ...
            'drift_config must be a structure.');
    end

    if isfield(drift_config,'max_baseline_range')
        max_baseline_range = ...
            double(drift_config.max_baseline_range);
    end

    if isfield(drift_config,'min_anchor_periods')
        min_anchor_periods = ...
            double(drift_config.min_anchor_periods);
    end

    if isfield(drift_config,'method')
        method = ...
            string(drift_config.method);
    end

    if isfield(drift_config,'nan_failed_traces')
        nan_failed_traces = ...
            logical(drift_config.nan_failed_traces);
    end

end

if ~isfinite(max_baseline_range) || ...
        max_baseline_range <= 0

    error( ...
        'traces_residual_drift_corr:InvalidThreshold', ...
        'max_baseline_range must be > 0.');

end

if ~isfinite(min_anchor_periods) || ...
        min_anchor_periods < 2

    error( ...
        'traces_residual_drift_corr:InvalidMinAnchors', ...
        'min_anchor_periods must be >= 2.');

end

if method ~= "linear"

    error( ...
        'traces_residual_drift_corr:UnsupportedMethod', ...
        'Only method="linear" is currently supported.');

end


%% ========================================================================
% 2. Build baseline periods
% =========================================================================

baseline_periods = {};

if isfield(protocol,'baseline') && ...
        isfield(protocol.baseline,'period') && ...
        ~isempty(protocol.baseline.period)

    for iPeriod = 1:numel(protocol.baseline.period)

        frames = ...
            double(protocol.baseline.period(iPeriod).frames(:)');

        frames = ...
            unique(frames( ...
                frames >= 1 & ...
                frames <= nFrames));

        if ~isempty(frames)
            baseline_periods{end+1} = frames; %#ok<AGROW>
        end

    end

elseif isfield(protocol,'baseline') && ...
        isfield(protocol.baseline,'frames') && ...
        ~isempty(protocol.baseline.frames)

    % Fallback:
    % split protocol.baseline.frames into contiguous runs.

    frames = ...
        unique(double(protocol.baseline.frames(:)'));

    frames = ...
        frames(frames >= 1 & frames <= nFrames);

    if ~isempty(frames)

        run_start = 1;

        breaks = ...
            find(diff(frames) > 1);

        run_end = ...
            [breaks, numel(frames)];

        for iRun = 1:numel(run_end)

            idx = ...
                run_start:run_end(iRun);

            baseline_periods{end+1} = ...
                frames(idx); %#ok<AGROW>

            run_start = ...
                run_end(iRun) + 1;

        end

    end

else

    error( ...
        'traces_residual_drift_corr:NoBaselineDefinition', ...
        'Protocol contains no usable baseline periods/frames.');

end

nPeriods = ...
    numel(baseline_periods);

if nPeriods < min_anchor_periods

    error( ...
        'traces_residual_drift_corr:TooFewBaselinePeriods', ...
        ['Protocol provides only %d baseline periods, but at least %d ' ...
         'are required.'], ...
        nPeriods, ...
        min_anchor_periods);

end


%% ========================================================================
% 3. Baseline anchor frames
% =========================================================================

anchor_frames = ...
    nan(1,nPeriods);

for iPeriod = 1:nPeriods

    anchor_frames(iPeriod) = ...
        median(baseline_periods{iPeriod});

end

anchor_times_sec = ...
    (anchor_frames - 1) ./ ...
    double(protocol.fps);

all_times_sec = ...
    ((1:nFrames) - 1) ./ ...
    double(protocol.fps);


%% ========================================================================
% 4. Calculate baseline-period medians
% =========================================================================

anchor_values = ...
    nan(nCells,nPeriods);

for iPeriod = 1:nPeriods

    frames = ...
        baseline_periods{iPeriod};

    anchor_values(:,iPeriod) = ...
        median( ...
            double(traces_nm(:,frames)), ...
            2, ...
            'omitnan');

end


%% ========================================================================
% 5. Preallocate outputs
% =========================================================================

traces_corrected = ...
    nan(size(traces_nm),'like',traces_nm);

drift_info = struct;

drift_info.anchor_frames = ...
    anchor_frames;

drift_info.anchor_times_sec = ...
    anchor_times_sec;

drift_info.anchor_values = ...
    anchor_values;

drift_info.baseline_range = ...
    nan(nCells,1);

drift_info.slope_per_sec = ...
    nan(nCells,1);

drift_info.intercept = ...
    nan(nCells,1);

drift_info.r2 = ...
    nan(nCells,1);

drift_info.post_correction_range = ...
    nan(nCells,1);

drift_info.n_valid_anchors = ...
    zeros(nCells,1,'uint16');

drift_info.fail = ...
    false(nCells,1);

drift_info.correction_applied = ...
    false(nCells,1);

drift_info.fail_reason = ...
    strings(nCells,1);

drift_info.threshold = ...
    max_baseline_range;


%% ========================================================================
% 6. Quantify drift and correct passing cells
% =========================================================================

for iCell = 1:nCells

    y = ...
        anchor_values(iCell,:);

    valid = ...
        isfinite(y) & ...
        isfinite(anchor_times_sec);

    nValid = ...
        nnz(valid);

    drift_info.n_valid_anchors(iCell) = ...
        uint16(nValid);


    %% --------------------------------------------------------------------
    % 6A. Enough valid baseline anchors?
    % ---------------------------------------------------------------------

    if nValid < min_anchor_periods

        drift_info.fail(iCell) = true;
        drift_info.fail_reason(iCell) = ...
            "too_few_valid_baseline_anchors";

        if ~nan_failed_traces
            traces_corrected(iCell,:) = ...
                traces_nm(iCell,:);
        end

        continue

    end


    %% --------------------------------------------------------------------
    % 6B. Baseline-drift range
    % ---------------------------------------------------------------------

    valid_y = ...
        y(valid);

    baseline_range = ...
        max(valid_y) - min(valid_y);

    drift_info.baseline_range(iCell) = ...
        baseline_range;


    % User-defined QC:
    % meaningful baseline must have drift STRICTLY LESS than 0.15.
    if ~isfinite(baseline_range) || ...
            baseline_range >= max_baseline_range

        drift_info.fail(iCell) = true;
        drift_info.fail_reason(iCell) = ...
            "baseline_drift_range_too_large";

        if ~nan_failed_traces
            traces_corrected(iCell,:) = ...
                traces_nm(iCell,:);
        end

        continue

    end


    %% --------------------------------------------------------------------
    % 6C. Fit linear drift using only baseline-period anchors
    % ---------------------------------------------------------------------

    x = ...
        anchor_times_sec(valid);

    yfit = ...
        valid_y;

    p = ...
        polyfit( ...
            x, ...
            yfit, ...
            1);

    slope_per_sec = ...
        p(1);

    intercept = ...
        p(2);

    fitted_anchor = ...
        polyval(p,x);

    residual_anchor = ...
        yfit - fitted_anchor;


    %% --------------------------------------------------------------------
    % 6D. R^2 diagnostic
    % ---------------------------------------------------------------------

    ss_res = ...
        sum(residual_anchor.^2,'omitnan');

    y_centered = ...
        yfit - mean(yfit,'omitnan');

    ss_tot = ...
        sum(y_centered.^2,'omitnan');

    if isfinite(ss_tot) && ss_tot > 0

        r2 = ...
            1 - ss_res ./ ss_tot;

    else

        % Flat baseline: slope is effectively zero and a linear correction
        % is still valid.
        r2 = ...
            NaN;

    end


    %% --------------------------------------------------------------------
    % 6E. Subtract fitted baseline drift from entire normalized trace
    % ---------------------------------------------------------------------

    fitted_full_baseline = ...
        polyval( ...
            p, ...
            all_times_sec);

    corrected_trace = ...
        double(traces_nm(iCell,:)) - ...
        fitted_full_baseline;


    if any(~isfinite(corrected_trace))

        drift_info.fail(iCell) = true;
        drift_info.fail_reason(iCell) = ...
            "nonfinite_corrected_trace";

        if ~nan_failed_traces
            traces_corrected(iCell,:) = ...
                traces_nm(iCell,:);
        end

        continue

    end


    %% --------------------------------------------------------------------
    % 6F. Measure post-correction baseline range
    % ---------------------------------------------------------------------

    corrected_anchor = ...
        nan(1,nPeriods);

    for iPeriod = 1:nPeriods

        corrected_anchor(iPeriod) = ...
            median( ...
                corrected_trace(baseline_periods{iPeriod}), ...
                'omitnan');

    end

    valid_post = ...
        isfinite(corrected_anchor);

    if nnz(valid_post) >= 2

        post_range = ...
            max(corrected_anchor(valid_post)) - ...
            min(corrected_anchor(valid_post));

    else

        post_range = ...
            NaN;

    end


    %% --------------------------------------------------------------------
    % 6G. Store successful correction
    % ---------------------------------------------------------------------

    traces_corrected(iCell,:) = ...
        cast( ...
            corrected_trace, ...
            'like', ...
            traces_nm);

    drift_info.slope_per_sec(iCell) = ...
        slope_per_sec;

    drift_info.intercept(iCell) = ...
        intercept;

    drift_info.r2(iCell) = ...
        r2;

    drift_info.post_correction_range(iCell) = ...
        post_range;

    drift_info.correction_applied(iCell) = ...
        true;

end


%% ========================================================================
% 7. Summary
% =========================================================================

nFailed = ...
    nnz(drift_info.fail);

nCorrected = ...
    nnz(drift_info.correction_applied);

fprintf( ...
    ['Residual baseline drift: %d / %d cells failed range QC ' ...
     '(threshold < %.3f); %d cells linearly corrected.\n'], ...
    nFailed, ...
    nCells, ...
    max_baseline_range, ...
    nCorrected);

end
