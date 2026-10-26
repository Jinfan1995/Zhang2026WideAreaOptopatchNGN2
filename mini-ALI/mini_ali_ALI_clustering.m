function [clust_idx,clust_cen,alimap] = mini_ali_ALI_clustering(COMs,factor,nrow,ncol,threshold)
% mini_ali_ALI_clustering Generate spike density map (or ALI map) and then
% cluster spikes in space
%
% input--
%   COMs: information of each spike, nspikes*4 matrix [x,y,time,dF_max]
%   factor: ALI map contains factor times more pixel in both row,col dimensions
%   nrow: rows of original movie
%   ncol: columns of original movie
%   threshold: density threshold that defines the cluster
% output--
%   clust_idx: cluster index for each spike, non-clustered spikes are 0
%   alimap: spike density map (or ALI map) using 2D histogram
%
% Created 05/23/2026 Jinfan zhang

% generate ALI map
% factor=4; % specify the resolution of the ALI map (4x higher than the original pixel resolution)
[cnt,cen]=mini_ali_density_map(COMs(:,1:2),[nrow,ncol],factor); %Count the number of APs in each spatial bin
% smooth the map a bit
alimap = imgaussfilt(cnt, 2);

% Automatically detect ALI clusters from the ALI map
% threshold=0.20; % detect peaks in the ALI map that are brighter than the value 'th'
map_th=alimap;
map_th(alimap<threshold)=0;
BW=imregionalmax(map_th,8);
[I,J]=find(BW);
clust_cen=[cen{1}(I)',cen{2}(J)']; % nclust x 2 center locations of detected clusters

% use watershed to seperate neghbouring cluster
map_th_neg = -map_th;
map_th_neg(map_th_neg == 0) = -Inf;
cluster_label = watershed(imimposemin(map_th_neg, BW));
cluster_label(alimap<threshold) = 0;
% generate cluster index for each single spike. index = 0 for spikes out of
% cluster
row_idx = round((COMs(:,1)-0.5) .* factor + 0.5);
col_idx = round((COMs(:,2)-0.5) .* factor + 0.5);
% protect boundaries
row_idx = max(1, min(size(cluster_label,1), row_idx));
col_idx = max(1, min(size(cluster_label,2), col_idx));
linear_idx = sub2ind(size(cluster_label), row_idx, col_idx);
clust_idx = cluster_label(linear_idx);