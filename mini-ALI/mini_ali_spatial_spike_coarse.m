function spk = ...
    mini_ali_spatial_spike_coarse( ...
    dF_AP, ...
    patch_mask, ...
    peak_time)

nSpikes = length(peak_time);

if nSpikes == 0

    spk = zeros(0,4);
    return

end


%% =========================================================
% Mask pixels outside segmented region
% ==========================================================

dF_flat = reshape( ...
    dF_AP, ...
    [], ...
    nSpikes);

mask_idx = find(patch_mask);

%% Rather than generate a full replicated mask,
% just inspect mask pixels
%
masked_dF = ...
    dF_flat(mask_idx,:);


%% =========================================================
% Find brightest spatial pixel for every spike simultaneously
% ==========================================================

[peak_value,local_idx] = ...
    max(masked_dF,[],1);

peak_linear_idx = ...
    mask_idx(local_idx);

[peak_row,peak_col] = ...
    ind2sub( ...
        size(patch_mask), ...
        peak_linear_idx);


%% Output
spk = [ ...
    peak_row(:), ...
    peak_col(:), ...
    peak_time(:), ...
    peak_value(:)];

end