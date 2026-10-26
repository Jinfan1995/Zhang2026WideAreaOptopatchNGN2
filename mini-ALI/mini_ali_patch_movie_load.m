function [patch_coordinate,patch_movie,patch_mask] = mini_ali_patch_movie_load(movie,regions,patch_size)
% mini_ali_patch_movie_load load patch movies based on region crops
% 
% input--
%   movie:   nrow x ncol x nframe, stack of raw movie,
%   regions:    regions from image segmation based on binary map. a sturct that contains regions.Centroid, regions.BoundingBox,
%   regions.PixelIdxList, regions.PixelList,
%   patch_size:   a vector for patch size, [height, width];
% output--
%   patch_coordinate:    the up left coordinate of patch in raw movie
%   patch_movie:   cropped sub movie
%   patch_mask:   binary mask for selected region pixels
% Created 04/21/2026 Jinfan Zhang 

sz = size(movie);

% 1) Decide the patch size for each footprint, based on patch size, or boudingbox. also consider the edge of movie. 
min_X = min(floor(regions.Centroid(1,2)-patch_size(2)*0.5),floor(regions.BoundingBox(1,2)+0.5));
min_Y = min(floor(regions.Centroid(1,1)-patch_size(1)*0.5),floor(regions.BoundingBox(1,1)+0.5));
max_X = max(floor(regions.Centroid(1,2)+patch_size(2)*0.5),floor(regions.BoundingBox(1,2)+0.5+regions.BoundingBox(1,4)));
max_Y = max(floor(regions.Centroid(1,1)+patch_size(1)*0.5),floor(regions.BoundingBox(1,1)+0.5+regions.BoundingBox(1,3)));
patch_X = max(min_X,1):min(max_X,sz(1));
patch_Y = max(min_Y,1):min(max_Y,sz(2));

% 2) Record the up left coordinate of each patch
patch_coordinate = [patch_X(1) patch_Y(1)];

% 3) Spatial crop the movie
patch_movie = movie(patch_X,patch_Y,:); 

% 4) Generate binary region mask

patch_mask = false( ...
    length(patch_X), ...
    length(patch_Y));
% Convert region pixel indices to global row/column
[row,col] = ind2sub( ...
    sz(1:2), ...
    regions.PixelIdxList);
% Keep pixels inside this patch
valid = ...
    row >= patch_X(1) & ...
    row <= patch_X(end) & ...
    col >= patch_Y(1) & ...
    col <= patch_Y(end);
row = row(valid);
col = col(valid);
% Convert to local coordinates
local_row = row - patch_X(1) + 1;
local_col = col - patch_Y(1) + 1;
local_idx = sub2ind( ...
    size(patch_mask), ...
    local_row, ...
    local_col);
patch_mask(local_idx) = true;

end

