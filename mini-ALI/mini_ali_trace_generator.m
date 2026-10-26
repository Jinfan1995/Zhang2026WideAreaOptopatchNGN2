function cell_traces = mini_ali_trace_generator(footprint,movie,Background,dfdx, dfdy, shift)

[nrow,ncol,frame_count] = size(movie);

nCells = numel(footprint);

%% =========================================================
% Build global summed footprint ONCE
% ==========================================================

footprint_sum = zeros( ...
    nrow*ncol, ...
    1, ...
    'single');

for i = 1:nCells

    idx = double(footprint{i}{1});
    values = single(footprint{i}{2});

    footprint_sum(idx) = ...
        footprint_sum(idx) + values;

end


%% =========================================================
% Preallocate output
% ==========================================================

cell_traces = zeros( ...
    nCells, ...
    frame_count, ...
    'single');


%% =========================================================
% Extract trace for each footprint
% ==========================================================

for i = 1:nCells

    if mod(i,50) == 0
        fprintf( ...
            'Trace extraction %d / %d\n', ...
            i,nCells);
    end


    %% Get sparse footprint

    idx = double(footprint{i}{1});
    values = single(footprint{i}{2});


    %% Convert global indices to row / column

    [row,col] = ind2sub( ...
        [nrow,ncol], ...
        idx);


    %% Bounding box

    patch_rows = min(row):max(row);
    patch_cols = min(col):max(col);

    nr = length(patch_rows);
    nc = length(patch_cols);


    %% Convert footprint coordinates to local patch coordinates

    local_row = row - patch_rows(1) + 1;
    local_col = col - patch_cols(1) + 1;

    local_idx = sub2ind( ...
        [nr,nc], ...
        local_row, ...
        local_col);


    %% Build SMALL local footprint only

    footprint_local = zeros( ...
        nr*nc,1,'single');

    footprint_local(local_idx) = values;


    %% Global indices corresponding to this local patch

    [R,C] = ndgrid( ...
        patch_rows, ...
        patch_cols);

    global_patch_idx = sub2ind( ...
        [nrow,ncol], ...
        R(:), ...
        C(:));


    %% Other overlapping footprints

    footprint_background = ...
        footprint_sum(global_patch_idx) - ...
        footprint_local;


    %% Two regressors

    B = [ ...
        footprint_local, ...
        footprint_background];


    %% Load movie patch

    footprint_movie = single( ...
        movie( ...
        patch_rows, ...
        patch_cols, ...
        :));


    %% Motion correction

    footprint_movie = ...
        mini_ali_regenerateMotionCorrectedSubmovie( ...
        footprint_movie, ...
        [patch_rows(1),patch_cols(1)], ...
        dfdx, ...
        dfdy, ...
        shift);


    %% Background subtraction

    footprint_movie = ...
        footprint_movie - ...
        single(Background( ...
        patch_rows, ...
        patch_cols));


    flattened_movie = reshape( ...
        footprint_movie, ...
        [], ...
        frame_count);


    %% Linear regression

    BtB = B' * B;
    BtA = B' * flattened_movie;

    lambda = ...
        1e-6 * trace(BtB) / size(B,2);

    region_traces = ...
        (BtB + lambda*eye(size(B,2),'single')) ...
        \ BtA;


    cell_traces(i,:) = ...
        region_traces(1,:);

end
end