function [BW, info] = mini_ali_brightness_threshold(I_bgsub, validMask, params)

if nargin < 3
    params = struct;
end

if ~isfield(params,'seg_smooth_sigma')
    params.seg_smooth_sigma = 1.5;
end

if ~isfield(params,'seg_high_sigma')
    params.seg_high_sigma = 5;
end

if ~isfield(params,'seg_low_sigma')
    params.seg_low_sigma = 2.5;
end

if ~isfield(params,'seg_min_area')
    params.seg_min_area = 20;
end


%% ------------------------------------------------------------
% 1. Slight smoothing
% ------------------------------------------------------------

Iseg = imgaussfilt( ...
    single(I_bgsub), ...
    params.seg_smooth_sigma);


%% ------------------------------------------------------------
% 2. Get only valid finite pixels
% ------------------------------------------------------------

valid = ...
    validMask & ...
    isfinite(Iseg);

x = Iseg(valid);


%% ------------------------------------------------------------
% 3. Estimate background center from histogram MODE
%
% Bright cells generate a positive tail, so median can shift upward.
% Histogram mode usually corresponds better to the background population.
% ------------------------------------------------------------

[counts,edges] = histcounts(x,256);

[~,idx] = max(counts);

bg_center = ...
    (edges(idx) + edges(idx+1))/2;


%% ------------------------------------------------------------
% 4. Estimate background noise from LEFT side only
%
% This avoids bright-cell contamination.
% ------------------------------------------------------------

left_values = ...
    x(x <= bg_center);

if numel(left_values) > 20

    sigma_bg = ...
        1.4826 * ...
        median(abs(left_values - bg_center));

else

    % fallback
    sigma_bg = ...
        1.4826 * mad(x,1);

end


sigma_bg = max(sigma_bg,eps('single'));


%% ------------------------------------------------------------
% 5. Two thresholds
% ------------------------------------------------------------

th_high = ...
    bg_center + ...
    params.seg_high_sigma * sigma_bg;

th_low = ...
    bg_center + ...
    params.seg_low_sigma * sigma_bg;


%% ------------------------------------------------------------
% 6. Strong bright-pixel seeds
% ------------------------------------------------------------

BW_seed = ...
    validMask & ...
    Iseg >= th_high;


%% ------------------------------------------------------------
% 7. Permissive candidate region
% ------------------------------------------------------------

BW_low = ...
    validMask & ...
    Iseg >= th_low;


%% ------------------------------------------------------------
% 8. Hysteresis growth
%
% Keep low-threshold pixels ONLY if they are connected to
% a strong high-threshold seed.
% ------------------------------------------------------------

BW = ...
    imreconstruct(BW_seed,BW_low);


%% ------------------------------------------------------------
% 9. Remove tiny noise
% ------------------------------------------------------------

BW = ...
    bwareaopen( ...
        BW, ...
        params.seg_min_area, ...
        8);


%% QC
info.bg_center = bg_center;
info.sigma_bg = sigma_bg;

info.th_high = th_high;
info.th_low = th_low;

info.BW_seed = BW_seed;
info.BW_low = BW_low;

end