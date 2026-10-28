function [traces_nm, norm_info] = traces_flip_norm(exp_corrected_traces,baseline_frames)

% 1. Input information
[nCells, nFrames] = size(exp_corrected_traces);
baseline_frames = baseline_frames(:);

% Remove duplicate frames
baseline_frames = unique(baseline_frames);

% Keep only valid frame indices
baseline_frames = baseline_frames( ...
    baseline_frames >= 1 & ...
    baseline_frames <= nFrames);

if isempty(baseline_frames)
    error('traces_flip_norm_v2:NoBaselineFrames', ...
        'No valid baseline frames were supplied.');
end

% 2. Preallocate outputs
traces_nm = nan(size(exp_corrected_traces), 'like', exp_corrected_traces);
norm_info.fail = false(nCells,1);
norm_info.baseline = nan(nCells,1);

% 3. Normalize each cell
for iCell = 1:nCells
    trace = exp_corrected_traces(iCell,:);
    try
        baseline_data = double(trace(baseline_frames));
        valid_baseline = isfinite(baseline_data);
        if ~any(valid_baseline)
            error('traces_flip_norm_v2:NoValidBaseline', 'No valid baseline values are available.');
        end
        baseline = mean(baseline_data(valid_baseline), 'omitnan');
        if ~isfinite(baseline)
           error('traces_flip_norm_v2:InvalidBaseline', 'Baseline is not finite.');
        end
        if abs(baseline) < eps
           error('traces_flip_norm_v2:ZeroBaseline', 'Baseline is zero or extremely close to zero.');
        end

        normalized_trace = (baseline - double(trace)) ./ baseline;
        if all(~isfinite(normalized_trace))
           error('traces_flip_norm_v2:InvalidNormalizedTrace', 'Normalized trace contains no finite values.');
        end

        traces_nm(iCell,:) = cast(normalized_trace, 'like', exp_corrected_traces);
        norm_info.baseline(iCell) = baseline;
  
    catch
        % Normalization failure
        % Do NOT remove this cell.
        traces_nm(iCell,:) = NaN;
        norm_info.fail(iCell) = true;
        norm_info.baseline(iCell) = NaN;
    end
end
