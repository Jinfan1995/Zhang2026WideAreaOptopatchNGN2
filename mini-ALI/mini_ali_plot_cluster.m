function fig_handle = mini_ali_plot_cluster(reference_image, COMs, clust_idx, ax)

% mini_ali_plot_cluster Plot spike locations, use different color for spikes
% beloning to different clusters
%
% input--
%   reference_image: image showing cell structure
%   COMs: nspike x 2, (row,col) coordinates of each spike
%   clust_idx: nspike x 1, cluster index of each spike 
%   ax, optional axes handle
% output--
%   fig_handle, figure handle
%   
% Template: Chen, Tsai-Wen et al.(2025)  
% Created 04/21/2026 Jinfan Zhang
%% =========================================================
% Create figure if axes not provided
% ==========================================================

if nargin < 4 || isempty(ax)

    fig_handle = figure( ...
        'Color','w', ...
        'Position',[100 100 1000 850]);

    ax = axes(fig_handle);

else

    fig_handle = ancestor(ax,'figure');

end


%% =========================================================
% Reference image
% ==========================================================

imshow( ...
    reference_image, ...
    [], ...
    'Parent',ax, ...
    'InitialMagnification','fit');

hold(ax,'on');


%% =========================================================
% Handle empty input
% ==========================================================

if isempty(COMs) || isempty(clust_idx)

    title(ax,'No localized spikes');

    axis(ax,'image');
    axis(ax,'off');

    return

end


%% =========================================================
% Cluster information
% ==========================================================

if isempty(clust_idx)
    nclust = 0;
else
    nclust = double(max(clust_idx));
end


%% =========================================================
% Colors
% ==========================================================

if nclust > 0

    cmap = hsv(nclust);

    % Reproducible shuffled colors
    rng(1);
    perm = randperm(nclust);
    cmap = cmap(perm,:);

else

    cmap = zeros(0,3);

end
%% =========================================================
% Plot spikes belonging to clusters
% ==========================================================

for i = 1:nclust

    idx = clust_idx == i;

    if ~any(idx)
        continue
    end


    %% Spike positions

    plot( ...
        ax, ...
        COMs(idx,2), ...
        COMs(idx,1), ...
        '.', ...
        'Color',cmap(i,:), ...
        'MarkerSize',4);


    %% Cluster center

    center_row = mean(COMs(idx,1));
    center_col = mean(COMs(idx,2));


    %% Cluster ID

    text( ...
        ax, ...
        center_col + 3, ...
        center_row + 3, ...
        num2str(i), ...
        'Color',cmap(i,:), ...
        'FontSize',8, ...
        'FontWeight','bold', ...
        'Clipping','on');

end


%% =========================================================
% Optional: show unclustered spikes
% ==========================================================

unclustered = clust_idx == 0;

if any(unclustered)

    plot( ...
        ax, ...
        COMs(unclustered,2), ...
        COMs(unclustered,1), ...
        '.', ...
        'Color',[0.6 0.6 0.6], ...
        'MarkerSize',2);

end


%% =========================================================
% Formatting
% ==========================================================

title( ...
    ax, ...
    sprintf( ...
        'ALI spike clusters: %d clusters | %d spikes', ...
        nclust, ...
        size(COMs,1)));

axis(ax,'image');
axis(ax,'off');

hold(ax,'off');

end