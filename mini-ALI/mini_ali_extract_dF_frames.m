function dF_frames = mini_ali_extract_dF_frames( ...
    movie, ...
    rowRange, ...
    colRange, ...
    frame_idx, ...
    median_filter, ...
    dfdx, ...
    dfdy, ...
    shift)

%% =========================================================
% Extract motion-corrected dF images only at selected frames
% ==========================================================

frame_idx = round(frame_idx(:));

nRequested = length(frame_idx);

nr = length(rowRange);
nc = length(colRange);

nFramesTotal = size(movie,3);

dF_frames = zeros( ...
    nr, ...
    nc, ...
    nRequested, ...
    'single');


%% Amount of temporal padding
%
% Using median_filter frames on either side is conservative
% and ensures movmedian has enough neighboring frames.
%
pad = double(median_filter);


for k = 1:nRequested

    t = frame_idx(k);

    %% Temporal window around requested frame

    t0 = max(1, t-pad);
    t1 = min(nFramesTotal, t+pad);


    %% Only load this tiny temporal window

    local_movie = movie(rowRange, colRange, t0:t1);


    %% Motion correct only these frames

    local_movie = ...
        mini_ali_regenerateMotionCorrectedSubmovie( ...
        local_movie, ...
        [rowRange(1),colRange(1)], ...
        dfdx, ...
        dfdy, ...
        shift(:,t0:t1));


    %% Calculate temporal baseline only within this window

    local_dF = -(local_movie - movmedian( local_movie, ...
            median_filter, ...
            3));


    %% Locate requested frame within local movie

    local_t = t - t0 + 1;

    dF_frames(:,:,k) = ...
        local_dF(:,:,local_t);

end

end