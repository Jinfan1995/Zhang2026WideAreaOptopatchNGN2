function [peak_time, corrected_trace] = mini_ali_patch_detect_spike_times( ...
    patch_movie, ...
    patch_coordinate, ...
    patch_mask, ...
    dfdx, ...
    dfdy, ...
    shift, ...
    median_filter_size, ...
    threshold_factor, ...
    minimal_interval)

%% =========================================================
% Basic information
% ==========================================================

[nr,nc,nFrames] = size(patch_movie);

mask_idx = find(patch_mask);

%% =========================================================
% 1. Generate raw average trace
% ==========================================================

movie_flat = reshape( ...
    patch_movie, ...
    nr*nc, ...
    nFrames);

% Use only pixels inside the cell mask
raw_trace = sum(single(movie_flat(mask_idx,:)), 1)./ numel(patch_mask);


%% =========================================================
% 2. Calculate average motion coefficients
% ==========================================================

rowRange = ...
    patch_coordinate(1):...
    patch_coordinate(1)+nr-1;

colRange = ...
    patch_coordinate(2):...
    patch_coordinate(2)+nc-1;

dfdx_sub = single( ...
    dfdx(rowRange,colRange));

dfdy_sub = single( ...
    dfdy(rowRange,colRange));


gx = sum( ...
    dfdx_sub(mask_idx), ...
    'all') ...
    ./ numel(patch_mask);

gy = sum( ...
    dfdy_sub(mask_idx), ...
    'all') ...
    ./ numel(patch_mask);


%% =========================================================
% 3. Motion-correct ONLY the trace
% ==========================================================

corrected_trace = ...
    raw_trace ...
    - gx .* single(shift(1,:)) ...
    - gy .* single(shift(2,:));


%% =========================================================
% 4. Temporal spike detection
% ==========================================================

dF_trace = ...
    medfilt1( ...
        corrected_trace, ...
        median_filter_size*4) ...
    - corrected_trace;


%% Remove edges
edgeFrames = round(median_filter_size/2)*4;

dF_trace(1:min(edgeFrames,nFrames)) = 0;

dF_trace( ...
    max(1,nFrames-edgeFrames+1):nFrames) = 0;


%% Noise estimate
std_noise = mad(dF_trace);

threshold = ...
    threshold_factor .* std_noise;

dF_trace(dF_trace < threshold) = 0;


%% Detect peaks
[~,peak_time] = findpeaks(dF_trace);


%% =========================================================
% 5. Remove peaks too close together
% ==========================================================

% minimal_interval = 15;

if ~isempty(peak_time)

    keep = [true, ...
        diff(peak_time) >= minimal_interval];

    peak_time = peak_time(keep);

end


peak_time = peak_time(:);

end