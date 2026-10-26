function [validMask, bgMap, Icorr, qc] = mini_ali_well_mask_and_background(I, params)
% WELL_MASK_AND_BACKGROUND
%
% Conservative well-edge detection and background estimation for MiniALI.
%
% This version is designed to strongly prefer FALSE NEGATIVES over
% FALSE POSITIVES: if the algorithm is not convinced that a true
% outside-well compartment exists, the entire FOV is kept.
%
% Key safeguards against false well-edge detection:
%   1. Only the largest border-connected dark component is considered.
%   2. A small candidate is accepted only with extremely strong optical
%      evidence that it is truly outside the well.
%   3. True outside-well regions must be substantially darker than the
%      in-well region.
%   4. True outside-well regions must also have substantially lower
%      high-frequency texture than the neuronal region.
%   5. The boundary must be a meaningful macroscopic structure.
%   6. If confidence is insufficient, validMask = true(size(I)).
%
%
% INPUT
% -----
% I
%   2-D reference image. Prefer a temporal mean/median image.
%
% params
%   Optional parameter structure.
%
%
% OUTPUT
% ------
% validMask
%   true  = safe inside-well pixels
%   false = outside well / rim / guard region
%
% bgMap
%   Smooth background estimate from mask-aware morphological opening.
%   Invalid pixels are NaN.
%
% Icorr
%   Background-corrected image. Invalid pixels are NaN.
%
% qc
%   Detailed QC information.
%
%
% Recommended starting parameters:
%
%   params.well_sigma = 30;
%   params.dark_fraction = 0.10;
%   params.inner_edge_exclusion_px = 100;  % reject 100 px inside detected well edge
%
%   params.bg_tile_px = 96;
%   params.bg_percentile = 20;
%   params.bg_smooth_sigma = 48;
%
%
% -------------------------------------------------------------------------


%% ========================================================================
% 0. Parameters
% ========================================================================

if nargin < 2 || isempty(params)
    params = struct;
end


% -------------------------------------------------------------------------
% Large-scale well detection
% -------------------------------------------------------------------------

params = set_default(params, 'well_sigma', 30);

params = set_default(params, 'dark_fraction', 0.10);

params = set_default(params, 'intensity_low_percentile', 5);
params = set_default(params, 'intensity_high_percentile', 90);


% -------------------------------------------------------------------------
% Candidate area requirements
%
% min_outside_fraction:
%   candidates >= this size can be accepted using the normal confidence
%   tests.
%
% min_small_outside_fraction:
%   smaller candidates can ONLY be accepted if they satisfy the stricter
%   "veryStrongOpticalEvidence" test below.
% -------------------------------------------------------------------------

params = set_default(params, 'min_outside_fraction', 0.08);
params = set_default(params, 'min_small_outside_fraction', 0.025);
params = set_default(params, 'max_outside_fraction', 0.90);


% -------------------------------------------------------------------------
% Border / geometry requirements
% -------------------------------------------------------------------------

params = set_default(params, 'min_border_outside_fraction', 0.03);

params = set_default(params, 'min_outside_depth_px', 40);

params = set_default(params, 'min_boundary_length_fraction', 0.05);

params = set_default(params, 'min_boundary_span_fraction', 0.30);


% -------------------------------------------------------------------------
% Global darkness
%
% outside_inside_intensity_ratio =
%
%       median(outside) / median(inside)
%
% Lower means outside is much darker.
% -------------------------------------------------------------------------

params = set_default(params, ...
    'max_outside_inside_intensity_ratio', 0.70);

% Small candidate requires MUCH stronger darkness.
params = set_default(params, ...
    'small_candidate_max_intensity_ratio', 0.50);


% -------------------------------------------------------------------------
% Texture test
%
% Real outside-well regions should normally be relatively featureless.
% False dark regions inside the well still contain neurons/neurites.
%
% texture_ratio =
%
%       median(texture outside) / median(texture inside)
% -------------------------------------------------------------------------

params = set_default(params, 'texture_sigma_px', 3);

params = set_default(params, ...
    'max_outside_inside_texture_ratio', 0.60);

% Small candidate requires much stronger texture suppression.
params = set_default(params, ...
    'small_candidate_max_texture_ratio', 0.40);


% -------------------------------------------------------------------------
% Local boundary confirmation
% -------------------------------------------------------------------------

params = set_default(params, 'edge_band_px', 15);

params = set_default(params, 'min_contrast_fraction', 0.08);

params = set_default(params, 'edge_gradient_ratio', 2.0);

params = set_default(params, 'edge_gradient_percentile', 70);


% Extra rejection distance INSIDE the detected well edge.
%
% Once a true well edge is detected, all pixels outside the well are
% rejected, AND this many additional pixels on the INNER side of the edge
% are rejected.  Set to 0 to keep pixels right up to the detected edge.
%
% Default here is 100 px for the first test requested.
%
% Backward compatibility: if an older caller supplies well_guard_px but not
% inner_edge_exclusion_px, use that value.
if ~isfield(params,'inner_edge_exclusion_px') || ...
        isempty(params.inner_edge_exclusion_px)

    if isfield(params,'well_guard_px') && ...
            ~isempty(params.well_guard_px)

        params.inner_edge_exclusion_px = ...
            params.well_guard_px;

    else

        params.inner_edge_exclusion_px = 100;

    end
end


% -------------------------------------------------------------------------
% Background estimation
%
% Background is estimated using the user's original MiniALI strategy:
%
%   1. grayscale morphological opening with a disk
%   2. Gaussian smoothing
%
% Both operations are MASK-AWARE.  Near the well boundary, only pixels
% inside validMask are allowed to contribute.  Pixels outside the well are
% never inserted as zeros and therefore cannot bias the background.
% -------------------------------------------------------------------------

% Radius of the disk used for morphological opening.
% This should generally be larger than the bright cellular structures that
% you want removed from the background estimate.
params = set_default(params, 'background_disk_radius', 20);

% Equivalent to:
%       Background = imgaussfilt(Background,20)
% in the original MiniALI workflow.
params = set_default(params, 'background_gaussian_sigma', 20);

% Small denominator threshold used by normalized masked Gaussian smoothing.
params = set_default(params, 'background_min_weight', 1e-6);



%% ========================================================================
% 1. Validate image
% ========================================================================

if ndims(I) ~= 2
    error('I must be a 2-D image.');
end

I = single(I);

[H,W] = size(I);



%% ========================================================================
% 2. Large-scale image for well detection
% ========================================================================

Ilarge = imgaussfilt( ...
    I, ...
    params.well_sigma, ...
    'Padding','replicate');



%% ========================================================================
% 3. Robust intensity range
% ========================================================================

pLow = prctile( ...
    Ilarge(:), ...
    params.intensity_low_percentile);

pHigh = prctile( ...
    Ilarge(:), ...
    params.intensity_high_percentile);

dynamicRange = ...
    max(pHigh - pLow, eps('single'));



%% ========================================================================
% 4. Initial dark-pixel mask
%
% This mask is ONLY a candidate generator.
% It does not by itself mean "outside well".
% ========================================================================

darkThreshold = ...
    pLow + ...
    params.dark_fraction * dynamicRange;

darkPixels = ...
    Ilarge <= darkThreshold;



%% ========================================================================
% 5. Border-connected dark pixels
% ========================================================================

borderMask = false(H,W);

borderMask(1,:)   = true;
borderMask(end,:) = true;
borderMask(:,1)   = true;
borderMask(:,end) = true;

marker = ...
    darkPixels & borderMask;

if any(marker(:))

    outsideAll = ...
        imreconstruct(marker, darkPixels);

else

    outsideAll = ...
        false(H,W);

end



%% ========================================================================
% 6. Keep only the LARGEST connected border-dark component
% ========================================================================

outsideCandidate = ...
    false(H,W);

CC = ...
    bwconncomp(outsideAll,8);

componentArea = [];

largestArea = 0;
largestIdx = NaN;

if CC.NumObjects > 0

    componentArea = ...
        cellfun(@numel, CC.PixelIdxList);

    [largestArea,largestIdx] = ...
        max(componentArea);

    outsideCandidate( ...
        CC.PixelIdxList{largestIdx}) = true;

end



%% ========================================================================
% 7. Candidate size / border connectivity
% ========================================================================

outsideFraction = ...
    nnz(outsideCandidate) / numel(outsideCandidate);

borderOutsideFraction = ...
    nnz(outsideCandidate & borderMask) / ...
    max(nnz(borderMask),1);

normalSizeCandidate = ...
    outsideFraction >= params.min_outside_fraction && ...
    outsideFraction <= params.max_outside_fraction;

smallSizeCandidate = ...
    outsideFraction >= params.min_small_outside_fraction && ...
    outsideFraction < params.min_outside_fraction;

strongBorderConnection = ...
    borderOutsideFraction >= ...
    params.min_border_outside_fraction;



%% ========================================================================
% 8. Candidate depth
%
% Small shallow edge patches are common false positives.
% ========================================================================

if any(outsideCandidate(:))

    outsideDepthMap = ...
        bwdist(~outsideCandidate);

    maxOutsideDepth = ...
        max(outsideDepthMap(:));

else

    outsideDepthMap = ...
        zeros(H,W,'single');

    maxOutsideDepth = 0;

end

deepOutsideRegion = ...
    maxOutsideDepth >= ...
    params.min_outside_depth_px;



%% ========================================================================
% 9. Define an INSIDE reference region
%
% Do not compare the outside candidate directly against pixels immediately
% adjacent to the boundary because the optical rim can contaminate them.
% ========================================================================

distanceFromOutside = ...
    bwdist(outsideCandidate);

insideReference = ...
    ~outsideCandidate & ...
    distanceFromOutside > ...
    max(2*params.edge_band_px, params.inner_edge_exclusion_px);

% Fallback for unusually small images
if nnz(insideReference) < 0.05*numel(I)

    insideReference = ...
        ~outsideCandidate;

end



%% ========================================================================
% 10. Global intensity-ratio test
%
% This is more intuitive than only using a local gradient:
% a true outside-well compartment should be substantially darker overall.
% ========================================================================

outsideMedian = NaN;
insideMedian = NaN;
outsideInsideIntensityRatio = NaN;

strongGlobalDarkness = false;
veryStrongGlobalDarkness = false;

if any(outsideCandidate(:)) && any(insideReference(:))

    outsideMedian = ...
        median(Ilarge(outsideCandidate));

    insideMedian = ...
        median(Ilarge(insideReference));

    if insideMedian > 0

        outsideInsideIntensityRatio = ...
            outsideMedian / insideMedian;

    else

        % If intensity values are not positive, fall back to a
        % dynamic-range-normalized contrast.
        outsideInsideIntensityRatio = ...
            1 - (insideMedian - outsideMedian) / dynamicRange;

    end

    strongGlobalDarkness = ...
        outsideInsideIntensityRatio <= ...
        params.max_outside_inside_intensity_ratio;

    veryStrongGlobalDarkness = ...
        outsideInsideIntensityRatio <= ...
        params.small_candidate_max_intensity_ratio;

end



%% ========================================================================
% 11. Texture / feature-content test
%
% This is the major safeguard against the false-positive case where a
% slightly darker corner STILL CONTAINS NEURONS.
%
% Build a high-pass image:
%
%       Ihigh = I - lowpass(I)
%
% Then use local absolute high-frequency energy as a texture map.
% ========================================================================

IlowTexture = ...
    imgaussfilt( ...
        I, ...
        params.texture_sigma_px, ...
        'Padding','replicate');

Ihigh = ...
    I - IlowTexture;

textureMap = ...
    imgaussfilt( ...
        abs(Ihigh), ...
        params.texture_sigma_px, ...
        'Padding','replicate');

outsideTexture = NaN;
insideTexture = NaN;
outsideInsideTextureRatio = NaN;

lowOutsideTexture = false;
veryLowOutsideTexture = false;

if any(outsideCandidate(:)) && any(insideReference(:))

    outsideTexture = ...
        median(textureMap(outsideCandidate));

    insideTexture = ...
        median(textureMap(insideReference));

    outsideInsideTextureRatio = ...
        outsideTexture / ...
        max(insideTexture, eps('single'));

    lowOutsideTexture = ...
        outsideInsideTextureRatio <= ...
        params.max_outside_inside_texture_ratio;

    veryLowOutsideTexture = ...
        outsideInsideTextureRatio <= ...
        params.small_candidate_max_texture_ratio;

end



%% ========================================================================
% 12. Internal boundary only
%
% Excludes the image-border portion of the candidate perimeter.
% ========================================================================

candidateBoundary = ...
    outsideCandidate & ...
    imdilate( ...
        ~outsideCandidate, ...
        strel('disk',1,0));

candidateBoundary(1,:)   = false;
candidateBoundary(end,:) = false;
candidateBoundary(:,1)   = false;
candidateBoundary(:,end) = false;


boundaryLength = ...
    nnz(candidateBoundary);

minimumBoundaryLength = ...
    params.min_boundary_length_fraction * ...
    min(H,W);

reasonableBoundaryLength = ...
    boundaryLength >= minimumBoundaryLength;



%% ========================================================================
% 13. Boundary span
%
% A true well boundary should normally be a macroscopic feature rather than
% a tiny local indentation.
% ========================================================================

boundaryXSpan = 0;
boundaryYSpan = 0;

boundaryXSpanFraction = 0;
boundaryYSpanFraction = 0;
maxBoundarySpanFraction = 0;

if any(candidateBoundary(:))

    [yb,xb] = ...
        find(candidateBoundary);

    boundaryXSpan = ...
        max(xb) - min(xb) + 1;

    boundaryYSpan = ...
        max(yb) - min(yb) + 1;

    boundaryXSpanFraction = ...
        boundaryXSpan / W;

    boundaryYSpanFraction = ...
        boundaryYSpan / H;

    maxBoundarySpanFraction = ...
        max( ...
            boundaryXSpanFraction, ...
            boundaryYSpanFraction);

end

reasonableBoundarySpan = ...
    maxBoundarySpanFraction >= ...
    params.min_boundary_span_fraction;



%% ========================================================================
% 14. Local contrast across boundary
% ========================================================================

contrast = NaN;

insideValue = NaN;
outsideValue = NaN;

insideBand = false(H,W);
outsideBand = false(H,W);

if any(candidateBoundary(:))

    distanceFromInside = ...
        bwdist(~outsideCandidate);

    insideBand = ...
        ~outsideCandidate & ...
        distanceFromOutside <= ...
        params.edge_band_px;

    outsideBand = ...
        outsideCandidate & ...
        distanceFromInside <= ...
        params.edge_band_px;

    if any(insideBand(:)) && ...
            any(outsideBand(:))

        insideValue = ...
            median(Ilarge(insideBand));

        outsideValue = ...
            median(Ilarge(outsideBand));

        contrast = ...
            insideValue - outsideValue;

    end

end

strongContrast = false;

if isfinite(contrast)

    strongContrast = ...
        contrast >= ...
        params.min_contrast_fraction * ...
        dynamicRange;

end



%% ========================================================================
% 15. Gradient evidence near candidate boundary
% ========================================================================

[Gmag,~] = ...
    imgradient(Ilarge);

gradientReference = ...
    prctile(Gmag(:),75) + ...
    eps('single');

edgeGradient = NaN;

edgeNeighbourhood = ...
    false(H,W);

if any(candidateBoundary(:))

    edgeNeighbourhood = ...
        imdilate( ...
            candidateBoundary, ...
            strel('disk',3,0));

    edgeGradientValues = ...
        Gmag(edgeNeighbourhood);

    if ~isempty(edgeGradientValues)

        edgeGradient = ...
            prctile( ...
                edgeGradientValues, ...
                params.edge_gradient_percentile);

    end

end

strongGradient = false;

if isfinite(edgeGradient)

    strongGradient = ...
        edgeGradient >= ...
        params.edge_gradient_ratio * ...
        gradientReference;

end



%% ========================================================================
% 16. Size-dependent optical confidence
%
% Normal-sized candidate:
%   darkness + low texture are required.
%
% Small candidate:
%   only accepted if it is VERY dark AND VERY textureless.
%
% This is the key behavior that prevents a dark neuronal corner from being
% interpreted as outside well.
% ========================================================================

normalOpticalEvidence = ...
    strongGlobalDarkness && ...
    lowOutsideTexture;

veryStrongOpticalEvidence = ...
    veryStrongGlobalDarkness && ...
    veryLowOutsideTexture;

sizeAndOpticsPass = ...
    (normalSizeCandidate && normalOpticalEvidence) || ...
    (smallSizeCandidate && veryStrongOpticalEvidence);



%% ========================================================================
% 17. Final edge decision
%
% Intentionally conservative.
% ========================================================================

edgeDetected = ...
    sizeAndOpticsPass && ...
    strongBorderConnection && ...
    deepOutsideRegion && ...
    reasonableBoundaryLength && ...
    reasonableBoundarySpan && ...
    (strongContrast || strongGradient);



%% ========================================================================
% 18. Construct final valid mask + inward edge exclusion
% ========================================================================

if edgeDetected

    % For pixels INSIDE the well, bwdist(outsideCandidate) is effectively
    % the distance to the detected well boundary.  Therefore this creates
    % an inward exclusion band of exactly the requested thickness.
    distanceToOutside = ...
        bwdist(outsideCandidate);

    innerEdgeExclusionMask = ...
        ~outsideCandidate & ...
        distanceToOutside <= ...
        params.inner_edge_exclusion_px;

    % Reject both the true outside-well region and the additional inner
    % safety band.
    rejectMask = ...
        outsideCandidate | ...
        innerEdgeExclusionMask;

    validMask = ...
        ~rejectMask;

else

    % If we are not confident, do NOT invent a boundary and do NOT apply
    % an inner-edge exclusion band.
    innerEdgeExclusionMask = ...
        false(H,W);

    rejectMask = ...
        false(H,W);

    validMask = ...
        true(H,W);

end



%% ========================================================================
% 19. Mask-aware morphological-opening background estimate
%
% Original MiniALI logic:
%
%       se = strel('disk',params.background_disk_radius);
%       Background = imopen(reference_image,se);
%
% The important modification here is MASK AWARENESS.
%
% A normal imopen() does not know that pixels outside the well should be
% ignored.  If the outside region is represented by zero/dark values, the
% erosion step can pull those values into the valid region near the well
% boundary and distort the estimated background.
%
% Here:
%
%   - erosion sees invalid pixels as +Inf, so invalid pixels cannot become
%     the local minimum;
%
%   - dilation sees invalid pixels as -Inf, so invalid pixels cannot become
%     the local maximum;
%
% Therefore, when the structuring element overlaps the well boundary, the
% background is estimated only from the AVAILABLE valid pixels inside the
% well.
% ========================================================================

backgroundMask = ...
    validMask & isfinite(I);

seBackground = ...
    strel( ...
        'disk', ...
        double(params.background_disk_radius), ...
        0);

Background_opened = ...
    masked_gray_open( ...
        I, ...
        backgroundMask, ...
        seBackground);



%% ========================================================================
% 20. Mask-aware Gaussian smoothing
%
% Original MiniALI logic:
%
%       Background = imgaussfilt(Background,20);
%
% A direct Gaussian filter is problematic at a mask boundary because:
%
%   - filling invalid pixels with 0 biases the result downward;
%   - filling invalid pixels with NaN propagates NaNs;
%   - ordinary filtering has no concept of "available pixels only".
%
% Instead, use normalized convolution:
%
%       smooth(background .* mask)
%       -------------------------
%             smooth(mask)
%
% Thus, each output pixel is averaged only from valid background pixels
% that actually exist within the Gaussian kernel.
% ========================================================================

if params.background_gaussian_sigma > 0

    bgForSmooth = ...
        Background_opened;

    % Invalid pixels contribute ZERO to the numerator.
    bgForSmooth(~backgroundMask) = 0;

    backgroundWeight = ...
        single(backgroundMask);

    bgNumerator = ...
        imgaussfilt( ...
            bgForSmooth, ...
            params.background_gaussian_sigma, ...
            'Padding','replicate');

    bgDenominator = ...
        imgaussfilt( ...
            backgroundWeight, ...
            params.background_gaussian_sigma, ...
            'Padding','replicate');

    bgMap = ...
        bgNumerator ./ ...
        max( ...
            bgDenominator, ...
            single(params.background_min_weight));

else

    bgDenominator = ...
        single(backgroundMask);

    bgMap = ...
        Background_opened;

end


% Only the valid well area is meaningful for downstream analysis.
bgMap(~validMask) = NaN;



%% ========================================================================
% 21. Background-estimation QC
% ========================================================================

% Difference between the raw reference image and the morphological-opening
% estimate before Gaussian smoothing.  Bright cells should largely appear
% here rather than in Background_opened.
backgroundOpeningResidual = ...
    I - Background_opened;

backgroundOpeningResidual(~validMask) = NaN;


%% ========================================================================
% 22. Background subtraction
% ========================================================================

Icorr = ...
    I - bgMap;

Icorr(~validMask) = NaN;

bgMap(~validMask) = NaN;



%% ========================================================================
% 23. QC
% ========================================================================

qc = struct;

% Final result
qc.edge_detected = edgeDetected;


% Candidate generation
qc.dark_threshold = darkThreshold;
qc.intensity_low = pLow;
qc.intensity_high = pHigh;
qc.dynamic_range = dynamicRange;

qc.n_outside_components = CC.NumObjects;
qc.component_area = componentArea;

qc.largest_component_index = largestIdx;
qc.largest_component_area = largestArea;

qc.outside_fraction = outsideFraction;
qc.border_outside_fraction = borderOutsideFraction;


% Size classification
qc.normal_size_candidate = normalSizeCandidate;
qc.small_size_candidate = smallSizeCandidate;


% Geometry
qc.max_outside_depth = maxOutsideDepth;
qc.min_required_outside_depth = ...
    params.min_outside_depth_px;

qc.boundary_length = boundaryLength;
qc.minimum_boundary_length = minimumBoundaryLength;

qc.boundary_x_span_fraction = ...
    boundaryXSpanFraction;

qc.boundary_y_span_fraction = ...
    boundaryYSpanFraction;

qc.max_boundary_span_fraction = ...
    maxBoundarySpanFraction;

qc.min_boundary_span_fraction = ...
    params.min_boundary_span_fraction;


% Global intensity
qc.outside_median = outsideMedian;
qc.inside_median = insideMedian;

qc.outside_inside_intensity_ratio = ...
    outsideInsideIntensityRatio;

qc.max_outside_inside_intensity_ratio = ...
    params.max_outside_inside_intensity_ratio;

qc.small_candidate_max_intensity_ratio = ...
    params.small_candidate_max_intensity_ratio;


% Texture
qc.outside_texture = outsideTexture;
qc.inside_texture = insideTexture;

qc.outside_inside_texture_ratio = ...
    outsideInsideTextureRatio;

qc.max_outside_inside_texture_ratio = ...
    params.max_outside_inside_texture_ratio;

qc.small_candidate_max_texture_ratio = ...
    params.small_candidate_max_texture_ratio;


% Local edge evidence
qc.inside_value = insideValue;
qc.outside_value = outsideValue;

qc.edge_contrast = contrast;

qc.min_required_contrast = ...
    params.min_contrast_fraction * ...
    dynamicRange;

qc.edge_gradient = edgeGradient;
qc.gradient_reference = gradientReference;

qc.min_required_gradient = ...
    params.edge_gradient_ratio * ...
    gradientReference;


% Boolean confidence components
qc.strong_border_connection = ...
    strongBorderConnection;

qc.deep_outside_region = ...
    deepOutsideRegion;

qc.reasonable_boundary_length = ...
    reasonableBoundaryLength;

qc.reasonable_boundary_span = ...
    reasonableBoundarySpan;

qc.strong_global_darkness = ...
    strongGlobalDarkness;

qc.very_strong_global_darkness = ...
    veryStrongGlobalDarkness;

qc.low_outside_texture = ...
    lowOutsideTexture;

qc.very_low_outside_texture = ...
    veryLowOutsideTexture;

qc.normal_optical_evidence = ...
    normalOpticalEvidence;

qc.very_strong_optical_evidence = ...
    veryStrongOpticalEvidence;

qc.size_and_optics_pass = ...
    sizeAndOpticsPass;

qc.strong_contrast = ...
    strongContrast;

qc.strong_gradient = ...
    strongGradient;


% Debug images
qc.large_scale_image = Ilarge;

qc.dark_pixels = darkPixels;
qc.outside_all = outsideAll;
qc.outside_candidate = outsideCandidate;

qc.outside_depth_map = outsideDepthMap;

qc.inside_reference = insideReference;

qc.texture_map = textureMap;

qc.candidate_boundary = candidateBoundary;

qc.inside_band = insideBand;
qc.outside_band = outsideBand;

qc.edge_neighbourhood = edgeNeighbourhood;

qc.reject_mask = rejectMask;
qc.inner_edge_exclusion_mask = innerEdgeExclusionMask;
qc.inner_edge_exclusion_px = params.inner_edge_exclusion_px;


% Background diagnostics
qc.background_mask = backgroundMask;
qc.background_opened = Background_opened;
qc.background_opening_residual = backgroundOpeningResidual;
qc.background_gaussian_weight = bgDenominator;

qc.background_disk_radius = ...
    params.background_disk_radius;

qc.background_gaussian_sigma = ...
    params.background_gaussian_sigma;


end



%% ========================================================================
% Helper: mask-aware grayscale morphological opening
% ========================================================================

function opened = masked_gray_open(I, mask, se)
% MASKED_GRAY_OPEN
%
% Performs a flat grayscale opening while allowing ONLY mask==true pixels
% to contribute to the morphology.
%
% Standard opening:
%
%       opened = imdilate(imerode(I,se),se)
%
% For erosion, invalid pixels are replaced with +Inf so they cannot become
% the neighborhood minimum.
%
% For dilation, invalid pixels are replaced with -Inf so they cannot become
% the neighborhood maximum.
%
% This is especially important when the structuring element touches the
% well boundary.

I = single(I);
mask = logical(mask);

opened = nan(size(I),'single');

if ~any(mask(:))
    return;
end


% -------------------------------------------------------------------------
% Erosion using only available valid pixels
% -------------------------------------------------------------------------

erosionInput = ...
    I;

erosionInput(~mask) = ...
    Inf;

eroded = ...
    imerode(erosionInput,se);


% -------------------------------------------------------------------------
% Dilation using only valid eroded pixels
% -------------------------------------------------------------------------

dilationInput = ...
    eroded;

dilationInput(~mask) = ...
    -Inf;

opened = ...
    imdilate(dilationInput,se);


% Outside the selected region is undefined.
opened(~mask) = ...
    NaN;


% Safety fallback for any pathological non-finite values that remain
% inside the mask.
bad = ...
    mask & ~isfinite(opened);

opened(bad) = ...
    I(bad);

end



%% ========================================================================
% Helper
% ========================================================================

function s = set_default(s,fieldName,value)

if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = value;
end

end
