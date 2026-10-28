function [traces_corrected, exp_fit_info] = ...
    traces_exp_corr( ...
        raw_traces, ...
        baseline_frames, ...
        median_filter_size, ...
        exp_config)
% TRACES_EXP_CORR
% Robust exponential photobleaching / slow-drift correction.
%
% The function fits only the median-filtered baseline portions of each raw
% trace, then divides the full raw trace by the fitted baseline.
%
% PRIMARY MODEL
%   exp2:
%       y = a*exp(b*x) + c*exp(d*x)
%
% FALLBACK
%   If exp2 is numerically or physically unstable, exp1 can be attempted:
%       y = a*exp(b*x)
%
% IMPORTANT FAILURE GUARDS
% -------------------------------------------------------------------------
% A mathematically successful exp2 fit can still be unusable. In particular,
% opposite-signed exponential components can nearly cancel at an intermediate
% frame. Dividing by that very small fitted baseline creates a huge artificial
% spike in the corrected trace.
%
% This function therefore rejects fits that:
%   1. contain non-finite fitted values
%   2. change sign relative to the observed baseline fluorescence
%   3. approach too close to zero relative to the observed baseline scale
%   4. produce non-finite corrected values
%   5. produce extreme corrected-trace amplification
%
% Failed cells remain in their original row but:
%
%       traces_corrected(iCell,:) = NaN
%       exp_fit_info.fail(iCell)  = true
%
% INPUT
% -------------------------------------------------------------------------
% raw_traces
%       nCells x nFrames raw fluorescence matrix.
%
% baseline_frames
%       Frames used to fit the slow fluorescence baseline.
%
% median_filter_size
%       Median-filter window in frames.
%
% exp_config (optional)
%       analysis_config.exp_fit structure.
%
%       Supported fields:
%           model
%           use_exp1_fallback
%           min_fitted_baseline_fraction
%           max_corrected_fold
%
% OUTPUT
% -------------------------------------------------------------------------
% traces_corrected
%
% exp_fit_info.fail
% exp_fit_info.rmse
% exp_fit_info.relative_rmse
% exp_fit_info.model_used
% exp_fit_info.used_fallback
% exp_fit_info.fail_reason
% exp_fit_info.min_fitted_baseline_fraction
% exp_fit_info.max_corrected_fold
%
% -------------------------------------------------------------------------

%% ========================================================================
% 1. Input information
% =========================================================================

[nCells,nFrames] = size(raw_traces);

baseline_frames = ...
    unique(baseline_frames(:));

baseline_frames = ...
    baseline_frames( ...
        baseline_frames >= 1 & ...
        baseline_frames <= nFrames);

if isempty(baseline_frames)
    error( ...
        'traces_exp_corr:NoBaselineFrames', ...
        'No valid baseline frames were supplied.');
end

if median_filter_size < 1
    error( ...
        'traces_exp_corr:InvalidMedianFilterSize', ...
        'median_filter_size must be >= 1.');
end


%% ========================================================================
% 2. Configuration
% =========================================================================

primary_model = "exp2";
use_exp1_fallback = true;
min_fitted_baseline_fraction = 0.05;
max_corrected_fold = 50;

if nargin >= 4 && ~isempty(exp_config)

    if ~isstruct(exp_config)
        error( ...
            'traces_exp_corr:InvalidConfig', ...
            'exp_config must be a structure.');
    end

    if isfield(exp_config,'model')
        primary_model = ...
            string(exp_config.model);
    end

    if isfield(exp_config,'use_exp1_fallback')
        use_exp1_fallback = ...
            logical(exp_config.use_exp1_fallback);
    end

    if isfield(exp_config,'min_fitted_baseline_fraction')
        min_fitted_baseline_fraction = ...
            double(exp_config.min_fitted_baseline_fraction);
    end

    if isfield(exp_config,'max_corrected_fold')
        max_corrected_fold = ...
            double(exp_config.max_corrected_fold);
    end

end

if ~ismember(primary_model,["exp1","exp2"])
    error( ...
        'traces_exp_corr:InvalidModel', ...
        'exp_fit.model must be "exp1" or "exp2".');
end

if ~isfinite(min_fitted_baseline_fraction) || ...
        min_fitted_baseline_fraction <= 0 || ...
        min_fitted_baseline_fraction >= 1

    error( ...
        'traces_exp_corr:InvalidBaselineFraction', ...
        'min_fitted_baseline_fraction must be between 0 and 1.');
end

if ~isfinite(max_corrected_fold) || ...
        max_corrected_fold <= 1

    error( ...
        'traces_exp_corr:InvalidCorrectedFold', ...
        'max_corrected_fold must be > 1.');
end


%% ========================================================================
% 3. Preallocate outputs
% =========================================================================

traces_corrected = ...
    nan(size(raw_traces),'like',raw_traces);

exp_fit_info = struct;

exp_fit_info.fail = ...
    false(nCells,1);

exp_fit_info.rmse = ...
    nan(nCells,1);

exp_fit_info.relative_rmse = ...
    nan(nCells,1);

exp_fit_info.model_used = ...
    strings(nCells,1);

exp_fit_info.used_fallback = ...
    false(nCells,1);

exp_fit_info.fail_reason = ...
    strings(nCells,1);

exp_fit_info.min_fitted_baseline_fraction = ...
    nan(nCells,1);

exp_fit_info.max_corrected_fold = ...
    nan(nCells,1);


%% ========================================================================
% 4. Frame coordinates
% =========================================================================
%
% Normalize the frame axis to approximately 0..1 before fitting.
% This improves numerical conditioning relative to fitting exponents directly
% on frame numbers as large as ~8000.

all_frames = ...
    (1:nFrames)';

if nFrames > 1

    all_x = ...
        (double(all_frames) - 1) ./ ...
        double(nFrames - 1);

else

    all_x = ...
        zeros(nFrames,1);

end

baseline_x = ...
    all_x(baseline_frames);


%% ========================================================================
% 5. Progress
% =========================================================================

progress_step = ...
    max(1,round(nCells/20));

fprintf( ...
    'Starting robust exponential correction for %d cells...\n', ...
    nCells);


%% ========================================================================
% 6. Process cells independently
% =========================================================================

for iCell = 1:nCells

    raw_trace = ...
        double(raw_traces(iCell,:));

    try

        %% ----------------------------------------------------------------
        % 6A. Median-filter raw trace
        % -----------------------------------------------------------------

        filtered_trace = ...
            medfilt1( ...
                raw_trace, ...
                median_filter_size, ...
                'omitnan', ...
                'truncate');


        %% ----------------------------------------------------------------
        % 6B. Baseline fitting data
        % -----------------------------------------------------------------

        y_baseline = ...
            filtered_trace(baseline_frames)';

        valid_idx = ...
            isfinite(baseline_x) & ...
            isfinite(y_baseline);

        x_fit = ...
            baseline_x(valid_idx);

        y_fit = ...
            y_baseline(valid_idx);

        if numel(y_fit) < 5
            error( ...
                'traces_exp_corr:TooFewPoints', ...
                'Too few valid baseline points.');
        end


        %% ----------------------------------------------------------------
        % 6C. Robust observed baseline scale/sign
        % -----------------------------------------------------------------

        baseline_center = ...
            median(y_fit,'omitnan');

        baseline_scale = ...
            median(abs(y_fit),'omitnan');

        if ~isfinite(baseline_center) || ...
                ~isfinite(baseline_scale) || ...
                baseline_scale <= 0 || ...
                baseline_center == 0

            error( ...
                'traces_exp_corr:InvalidObservedBaseline', ...
                'Observed baseline has invalid scale/sign.');

        end


        %% ----------------------------------------------------------------
        % 6D. Try primary model
        % -----------------------------------------------------------------

        [ok, fit_result, fitted_baseline, gof, ...
            fit_reason, fit_stats] = ...
            local_try_exponential_fit( ...
                x_fit, ...
                y_fit, ...
                all_x, ...
                primary_model, ...
                baseline_center, ...
                baseline_scale, ...
                min_fitted_baseline_fraction);


        %% ----------------------------------------------------------------
        % 6E. Optional exp1 fallback
        % -----------------------------------------------------------------

        used_fallback = false;

        if ~ok && ...
                use_exp1_fallback && ...
                primary_model == "exp2"

            [ok2, fit_result2, fitted_baseline2, gof2, ...
                fit_reason2, fit_stats2] = ...
                local_try_exponential_fit( ...
                    x_fit, ...
                    y_fit, ...
                    all_x, ...
                    "exp1", ...
                    baseline_center, ...
                    baseline_scale, ...
                    min_fitted_baseline_fraction);

            if ok2

                ok = true;

                fit_result = ...
                    fit_result2; %#ok<NASGU>

                fitted_baseline = ...
                    fitted_baseline2;

                gof = ...
                    gof2;

                fit_reason = ...
                    "";

                fit_stats = ...
                    fit_stats2;

                model_used = ...
                    "exp1";

                used_fallback = ...
                    true;

            else

                fit_reason = ...
                    "exp2: " + fit_reason + ...
                    " | exp1: " + fit_reason2;

            end

        else

            model_used = ...
                primary_model;

        end


        %% ----------------------------------------------------------------
        % 6F. If both models failed validation -> QC failure
        % -----------------------------------------------------------------

        if ~ok

            error( ...
                'traces_exp_corr:UnstableFit', ...
                '%s', ...
                char(fit_reason));

        end


        %% ----------------------------------------------------------------
        % 6G. Correct raw trace
        % -----------------------------------------------------------------

        corrected_trace = ...
            raw_trace(:) ./ fitted_baseline;


        %% ----------------------------------------------------------------
        % 6H. Validate corrected trace
        % -----------------------------------------------------------------

        if any(~isfinite(corrected_trace))

            error( ...
                'traces_exp_corr:InvalidCorrectedTrace', ...
                'Corrected trace contains non-finite values.');

        end

        corrected_baseline_scale = ...
            median( ...
                abs(corrected_trace(baseline_frames)), ...
                'omitnan');

        if ~isfinite(corrected_baseline_scale) || ...
                corrected_baseline_scale <= 0

            error( ...
                'traces_exp_corr:InvalidCorrectedBaseline', ...
                'Corrected baseline has invalid scale.');

        end

        corrected_max_fold = ...
            max(abs(corrected_trace),[],'omitnan') ./ ...
            corrected_baseline_scale;

        exp_fit_info.max_corrected_fold(iCell) = ...
            corrected_max_fold;

        if ~isfinite(corrected_max_fold) || ...
                corrected_max_fold > max_corrected_fold

            error( ...
                'traces_exp_corr:ExtremeCorrectedTrace', ...
                ['Corrected trace amplification %.3g exceeds allowed ' ...
                 'fold %.3g.'], ...
                corrected_max_fold, ...
                max_corrected_fold);

        end


        %% ----------------------------------------------------------------
        % 6I. Store successful result
        % -----------------------------------------------------------------

        traces_corrected(iCell,:) = ...
            cast( ...
                corrected_trace', ...
                'like', ...
                raw_traces);

        exp_fit_info.rmse(iCell) = ...
            gof.rmse;

        exp_fit_info.relative_rmse(iCell) = ...
            gof.rmse ./ baseline_scale;

        exp_fit_info.model_used(iCell) = ...
            model_used;

        exp_fit_info.used_fallback(iCell) = ...
            used_fallback;

        exp_fit_info.min_fitted_baseline_fraction(iCell) = ...
            fit_stats.min_abs_fitted ./ baseline_scale;


    catch ME

        %% ---------------------------------------------------------------
        % Fit / correction failure
        % ----------------------------------------------------------------

        traces_corrected(iCell,:) = ...
            NaN;

        exp_fit_info.fail(iCell) = ...
            true;

        exp_fit_info.rmse(iCell) = ...
            NaN;

        exp_fit_info.relative_rmse(iCell) = ...
            NaN;

        if strlength(exp_fit_info.model_used(iCell)) == 0
            exp_fit_info.model_used(iCell) = ...
                "";
        end

        exp_fit_info.fail_reason(iCell) = ...
            string(ME.identifier) + ": " + string(ME.message);

    end


    %% --------------------------------------------------------------------
    % Progress
    % ---------------------------------------------------------------------

    if mod(iCell,progress_step) == 0 || ...
            iCell == 1 || ...
            iCell == nCells

        fprintf( ...
            ['Exponential correction: %d/%d (%.1f%%) | ' ...
             'failed: %d | exp1 fallback: %d\n'], ...
            iCell, ...
            nCells, ...
            100*iCell/nCells, ...
            nnz(exp_fit_info.fail(1:iCell)), ...
            nnz(exp_fit_info.used_fallback(1:iCell)));

    end

end


%% ========================================================================
% 7. Summary
% =========================================================================

nFailed = ...
    nnz(exp_fit_info.fail);

nFallback = ...
    nnz(exp_fit_info.used_fallback);

fprintf( ...
    ['Exponential correction: %d / %d cells failed (%.2f%%); ' ...
     '%d used exp1 fallback.\n'], ...
    nFailed, ...
    nCells, ...
    100*nFailed/nCells, ...
    nFallback);

end


%% ========================================================================
% Local helper
% =========================================================================

function [ok,fit_result,fitted_baseline,gof,reason,stats] = ...
    local_try_exponential_fit( ...
        x_fit, ...
        y_fit, ...
        all_x, ...
        model_name, ...
        baseline_center, ...
        baseline_scale, ...
        min_baseline_fraction)

ok = false;

fit_result = [];
fitted_baseline = [];
gof = struct('rmse',NaN);
reason = "";
stats = struct('min_abs_fitted',NaN);

try

    [fit_result,gof] = ...
        fit( ...
            x_fit, ...
            y_fit, ...
            char(model_name));

    fitted_baseline = ...
        feval( ...
            fit_result, ...
            double(all_x));

catch ME

    reason = ...
        "fit exception: " + string(ME.message);

    return

end


% Non-finite baseline
if any(~isfinite(fitted_baseline))

    reason = ...
        "fitted baseline contains non-finite values";

    return

end


% The fitted fluorescence baseline should keep the same sign as the
% robustly observed raw baseline. A sign change implies a zero crossing.
baseline_sign = ...
    sign(baseline_center);

if any(fitted_baseline .* baseline_sign <= 0)

    reason = ...
        "fitted baseline changes sign / crosses zero";

    return

end


% Relative near-zero guard.
% Using eps alone is much too permissive for fluorescence data.
min_abs_fitted = ...
    min(abs(fitted_baseline),[],'omitnan');

stats.min_abs_fitted = ...
    min_abs_fitted;

minimum_allowed = ...
    min_baseline_fraction .* baseline_scale;

if min_abs_fitted < minimum_allowed

    reason = ...
        sprintf( ...
            ['fitted baseline approaches zero: min=%.4g, ' ...
             'required >= %.4g'], ...
            min_abs_fitted, ...
            minimum_allowed);

    return

end

ok = true;

end
