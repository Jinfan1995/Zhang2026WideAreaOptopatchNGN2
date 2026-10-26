function [spike_fine]=mini_ali_spk_fine(dF_AP,npix,radius,spk_coarse)
% mini_ali_spk_fine Localize AP with sub-pixel precision through weighted averagaing
%
% input--
%   dF_AP: nrow x ncol x nspike, stack of df_AP images, contain positively going spikes
%   npix: scalar, number of pixels around the peak for calculating the average
%   radius: scalar
%   spk_coarse: nspike x 3, coarse coordinates [row,col,frame] of action potentials, integer pixel values
% output--
%   spike_fine: nspike x 3 array
%         sloc(i,1) stores the row location of the i th spike;
%         sloc(i,2) stores the col location of the i th spike;
%         sloc(i,3) stores the frame number of the i th spike;
%   spike_BW: nrow x ncol x nspike, BW images of the region included in the average
% Template: Chen, Tsai-Wen et al.(2025)  
% Created 04/21/2026 Jinfan Zhang

sz = size(dF_AP);

spike_fine = spk_coarse;
for i=1:size(dF_AP,3)
    current_AP = dF_AP(:,:,i);
    picked_pix = zeros(npix,2);
    % initial pixels form spike detection
    [~, idx] = max(current_AP,[],'all','linear');
    [initial_x, initial_y] = ind2sub(size(current_AP), idx);
    picked_pix(1,:) = [initial_x,initial_y];
    % Search npix brightest pixels within the radous around the initial pixel      
    pix_count = 2;
    pixel_qeue = [];
    while pix_count <= npix
        % add 4 neighbouring pixels from last picked pixel for searching
        pix = picked_pix(pix_count-1,:)';
        search_coordinate = [pix + [[1 0]',[-1 0]', [0 1]',[0 -1]']];search_coordinate = search_coordinate';
        % delete the pixels that is beyond the searching radious
        delta_vector = search_coordinate - picked_pix(1,:);
        distance = delta_vector(:,1).^2 + delta_vector(:,2).^2;
        search_coordinate(distance>radius^2,:)=[];
        % delete pixels over the patch boundary
        delete_x = search_coordinate(:,1)<=0 | search_coordinate(:,1)> sz(1);
        delete_y = search_coordinate(:,2)<=0 | search_coordinate(:,2)> sz(2);
        search_coordinate(delete_x|delete_y,:) = [];
        % delete pixels from 
        % add new pixels to search qeue
        pixel_qeue = [pixel_qeue;search_coordinate];
        % Delete pixels already picked
        overlap_pix = [];
        for j = 1:size(pixel_qeue,1)
            if ismember(pixel_qeue(j,:),picked_pix,'rows')
                overlap_pix = [overlap_pix,j];
            end
        end
        pixel_qeue(overlap_pix,:)=[];
        if isempty(pixel_qeue)
            break
        end

        % find the pixel in the qeue with max value
        queue_idx = sub2ind(size(current_AP), pixel_qeue(:,1), pixel_qeue(:,2));
        [~,max_idx] = max(current_AP(queue_idx));
        picked_pix(pix_count,:) = pixel_qeue(max_idx,:);

        pix_count = pix_count + 1;
    end
    
    delete_index = picked_pix(:,1)==0 | picked_pix(:,2)==0;
    picked_pix(delete_index,:) = [];
    
    picked_pix_x = picked_pix(:,1);
    picked_pix_y = picked_pix(:,2);

    spk_idx = sub2ind(size(current_AP), picked_pix_x, picked_pix_y);  % NOTE: (row=y, col=x)

    values = double(current_AP(spk_idx));
    values(values < 0) = 0;
    weights = values.^2;
    total_mass = sum(weights);

    if total_mass > 0
        x_center = sum(weights .* double(picked_pix_x)) / total_mass;
        y_center = sum(weights .* double(picked_pix_y)) / total_mass;
    else
        x_center = NaN;
        y_center = NaN;
    end

    spike_fine(i,[1 2]) = [x_center,y_center];
end

        


    
    
    
    
    
    
    
    
    
    
    
    
    
 