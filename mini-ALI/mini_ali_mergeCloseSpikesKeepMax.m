function merged_spikes = mini_ali_mergeCloseSpikesKeepMax(spikes, distThresh)
% mergeCloseSpikesKeepMax
%
% INPUT:
%   spikes      N x 4 matrix: [row, col, time, spike_value]
%   distThresh  distance threshold, e.g. 8 pixels
%
% OUTPUT:
%   merged_spikes M x 4 matrix
%   For each same-time spatial cluster
%     1. Group spatially close spikes.
%     2. Identify connected spatial clusters.
%     3. Keep only the spike with the largest spike_value within each cluster.

    if nargin < 2
        distThresh = 8;
    end

    if isempty(spikes)
        merged_spikes = spikes;
        return
    end

    % Sort by time
    spikes = sortrows(spikes, 3);

    % Find contiguous blocks having the same time
    change_idx = [1; find(diff(spikes(:,3)) ~= 0) + 1; size(spikes,1) + 1];
    nTimes = numel(change_idx) - 1;
    merged_cell = cell(nTimes, 1);

    % Process each time block
    for k = 1:nTimes

        idx1 = change_idx(k);
        idx2 = change_idx(k+1) - 1;
    
        % All spikes in S occur at exactly the same time
        S = spikes(idx1:idx2, :);
    
        % Only one spike at this time
        if size(S,1) == 1
           merged_cell{k} = S;
           continue
        end

        % Spatial grouping

        xy = S(:,1:2);
        n = size(xy,1);

        if n == 1       
            group_id = 1;   
        else
            % Find nearby spikes only
            idxNeighbors = rangesearch(xy, xy, distThresh);    
            edge_i = cell(n,1);
            edge_j = cell(n,1);
        
            for ii = 1:n
        
                neighbors = idxNeighbors{ii};
        
                % j > ii simultaneously:
                % 1. removes self
                % 2. avoids duplicate edges
                neighbors = neighbors(neighbors > ii);
        
                edge_i{ii} = repmat(ii, numel(neighbors), 1);
                edge_j{ii} = neighbors(:);
        
            end
        
            edge_i = vertcat(edge_i{:});
            edge_j = vertcat(edge_j{:});
        
            if isempty(edge_i)
        
                % No spikes are spatially connected
                group_id = 1:n;
        
            else
        
                G = graph(edge_i, edge_j, [], n);
        
                group_id = conncomp(G);
        
            end
    
        end
        % Keep strongest spike from each spatial group
        nGroup = max(group_id);
        S_merged = zeros(nGroup, size(S,2));

        for g = 1:nGroup
            members = find(group_id == g);
            % Find maximum spike_value (column 4)
            [~, localMaxIdx] = max(S(members,4));
            bestIdx = members(localMaxIdx);
            % Keep original spike information
            S_merged(g,:) = S(bestIdx,:);
        end
        merged_cell{k} = S_merged;
    end

 % Combine results
 merged_spikes = vertcat(merged_cell{:});

 % Technically already ordered by time, but keep this for safety
 merged_spikes = sortrows(merged_spikes, 3);

end













   
