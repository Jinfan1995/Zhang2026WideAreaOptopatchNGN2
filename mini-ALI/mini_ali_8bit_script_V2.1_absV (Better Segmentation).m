%% =========================================================
% Experiment setup
% ==========================================================

% Data needs to be organized as EXP_folder/Plate_1/Well_1/movie_folder1
% root_folder is your EXP_folder
root_path = uigetdir();

experiment_id = "EXP001_iPSC_NGN2_2061_CSMV2_screening"; % change accordingly
pipeline_version = "miniALI_v2.1_absV";

% Find all FOVs
fov_table = find_optopatch_fovs(root_path);
fprintf('Found %d FOVs\n', height(fov_table));

% Analysis directory
analysis_root = fullfile(root_path, "Analysis");

if ~exist(analysis_root,'dir')
    mkdir(analysis_root);
end

% save the FOV index
parquetwrite( ...
    fullfile(analysis_root,"fov_index.parquet"), ...
    fov_table);

%% =========================================================
% Loop through FOVs
% ==========================================================

%% Process each FOV
for iFOV = 1:height(fov_table)

    plate_name = fov_table.plate(iFOV);
    well_name = fov_table.well(iFOV);
    fov_number = fov_table.fov(iFOV);
    fov_name = fov_table.fov_name(iFOV);
    fov_path = fov_table.fov_path(iFOV);

    plate_number = str2double( ...
    regexp(plate_name, '\d+', 'match', 'once'));

    % Unique FOV identity
    % Global sequential FOV index within this experiment
    fov_index = uint32(iFOV);
    % Unique FOV ID
    fov_id = sprintf( ...
        '%s_P%02d_%s_F%03d', ...
        experiment_id, ...
        plate_number, ...
        well_name, ...
        fov_number);    
    fov_id = string(fov_id);

    fprintf( ...
        '\nProcessing plate %d / well %s / FOV %d (%d/%d)\n', ...
        plate_number, ...
        well_name, ...
        fov_number, ...
        iFOV, ...
        height(fov_table));


    %% Create corresponding analysis folder

    sub_analysis_folder = fullfile( ...
        analysis_root, ...
        plate_name, ...
        "Well_" + well_name, ...
        fov_name);

    if ~exist(sub_analysis_folder,'dir')
        mkdir(sub_analysis_folder);
    end
    %% Load video info
    % For Luminos
    cameradata = load(fullfile(fov_path,"output_data.mat"));
    movie_path = fullfile(fov_path,"frames1.bin");
    nrow = double(cameradata.Device_Data{1,2}.ROI(4));  
    ncol = double(cameradata.Device_Data{1,2}.ROI(2));   

    bit = 8;
    frame_count = cameradata.Device_Data{1, 2}.frames_requested; % Number of frames in the movie
    fps = 1 / 1.25 * 1000; % Frames per second
    time = (0:frame_count-1) / fps; % Time vector

    %% FOV metadata
    fov_metadata = struct;
    
    % Identity
    fov_metadata.experiment_id = experiment_id;
    
    fov_metadata.plate = uint16(plate_number);
    fov_metadata.well = well_name;
    fov_metadata.fov = uint16(fov_number);
    
    fov_metadata.fov_id = fov_id;
    fov_metadata.fov_index = fov_index;
    
    % Source
    fov_metadata.source_folder = fov_name;
    fov_metadata.fov_path = fov_path;
    fov_metadata.movie_path = movie_path;
    
    % Imaging information
    fov_metadata.height_px = uint16(nrow);
    fov_metadata.width_px = uint16(ncol);
    
    fov_metadata.frames_requested = uint32(frame_count);
    
    fov_metadata.fps = single(fps);
    fov_metadata.bit_depth = uint8(bit);
    % information for actual camera exposure time and DMD on time
    fov_metadata.CamExp = 0.00123; % sec
    fov_metadata.Dmd_on_time = 0.0005; % sec
    
    % Analysis
    fov_metadata.pipeline_version = pipeline_version;

    %% miniALI analysis parameters
    
    params = struct;

    % reference image
    params.reference_frames = uint16([1 300]);
    params.Optopatch_end_frame = 8399;
    params.DMD_off_frames = uint16([8400 8410]);
    params.pulse_frames = uint16([8411 8428]);

    
    % segmentation
    % params.intensity_percentile = single(92);% try 97–99
    params.watershed_percentile = single(95);
    params.background_disk_radius = uint16(50);% radius > neuron size
    params.background_gaussian_sigma = 40;
    params.inner_edge_exclusion_px = 125;
    params.seg_smooth_sigma = 2;
    params.seg_high_sigma = 6; % 4-5 for loosing to much cells, 7 for too much false positive
    params.seg_low_sigma = 3;% 2-2.5 for loosing to much cells, 3.5 for too much false positive
    params.seg_min_area = 15;

    % motion correction
    params.motion_reference_frames = uint16([1 800]);
    
    % trace processing
    params.median_filter_frames = uint16(8);
    
    % spike detection
    params.threshold_factor = single(4);
    params.nsvd = uint16(10);

    % Fine spike localization
    params.fine_npix = uint16(25);
    params.fine_radius_px = single(6);
    
    % Temporal spike detection
    params.min_spike_interval_frames = uint16(25);
    
    % spatial
    params.patch_size = uint16([27 27]);
    params.spike_merge_distance_px = single(8);
    
    % clustering
    params.ALI_cluster_factor = single(4);
    params.ALI_cluster_threshold = single(0.15);

    % Footprint
    params.footprint_radius_px = single(6);
    
    % store all parameters to metadata
    fov_metadata.analysis_parameters = params;
    
    median_filter = params.median_filter_frames;
    threshold_factor = params.threshold_factor;
    nsvd = params.nsvd;
    patch_size = double(params.patch_size);

    %% Load whole movie
    % load the raw 8-bit movie
    movie = readBinMov(movie_path, nrow, ncol,bit);
    % No end frames truncation comapred to raw V2.0 version
    % movie(:,:,end-9:end) = [];
    
    
    %% Load DMD Off background and substract
    reference_image_DMD_off = median (movie(:,:,params.DMD_off_frames(1):params.DMD_off_frames(2)),3);
    reference_image_DMD_off = uint8(reference_image_DMD_off);
    movie = movie -reference_image_DMD_off;

    %% Load reference image for cell segmentation
    % mean image of first 300 frames of movie for reference
    reference_image = mean (movie(:,:,params.reference_frames(1):params.reference_frames(2)),3);
    % Extract background
    se = strel('disk', double(params.background_disk_radius));   
    [validMask, Background, reference_corrected, qc] = mini_ali_well_mask_and_background(reference_image, params);
    % correct for negative valuee
    reference_seg = reference_corrected;
    reference_seg(reference_seg < 0) = 0;
    % reference_corrected(reference_corrected<0)=0;
    % Background = reference_image -reference_corrected;

    %% movie truncation, exclude dark frames and 
    pulse_start = params.pulse_frames(1);
    pulse_end   = params.pulse_frames(2);
    idx1 = pulse_start:2:pulse_end;
    idx2 = pulse_start+1:2:pulse_end;

    scale_factor = single(fov_metadata.CamExp / fov_metadata.Dmd_on_time);
    movie_pulses_combined = (single(movie(:,:,idx1)) + single(movie(:,:,idx2))) .* scale_factor;
    movie_pulses_combined = uint8(movie_pulses_combined);

    movie = cat(3, movie(:,:,1:params.DMD_off_frames(1)-1), movie_pulses_combined);

    clear movie_pulses_combined

    frame_count = size(movie,3);
    fov_metadata.frames_analyzed = uint32(frame_count);

    %% regionprops
    % Binarize the image based on intensity
    [BW_intensity,info] = mini_ali_brightness_threshold(reference_seg,validMask,params);
    
    % Binarize the image based on watershed of neuron soma
    I = reference_image;
    I = imgaussfilt(I, 1);  % reduce noise
    I_soma = imtophat(I, strel('disk', 3));% highlight the soma
    % background substraction
    I_bg = imopen(I_soma, strel('disk',double(params.background_disk_radius)));  % bigger scale
    I_bg = imgaussfilt(I_bg, 20); 
    I_corr = I_soma - I_bg;
    I_corr = mat2gray(I_corr);
    % generate binary mask
    th = prctile(I_corr(:), params.watershed_percentile);   % try 97–99
    BW = I_corr > th;
    BW = bwareaopen(BW, 8);     % remove small junk
    BW = imclose(BW, strel('disk',2));
    D = bwdist(~BW);
    D_smooth = imgaussfilt(D, 1);
    % Find regional maxima (cell centers)
    markers = imextendedmax(I_corr, 0.03);
    D_neg = -D_smooth;
    D_neg(~BW) = -Inf;
    L = watershed(imimposemin(D_neg, markers));
    L(~BW) = 0;   % remove background
    BW_watershed = bwareaopen(L>0, 2);     % remove small junk
    
    % Combine two binary image
    Binary_mask = BW_intensity & BW_watershed;
    % label all regions
   regions = regionprops(Binary_mask, 'BoundingBox', 'Centroid', 'PixelIdxList');
    
     % Display the binary mask with labeling
    % figure;
    % imshow(Binary_mask);
    % hold on;
    % for b = 1:length(regions)
    % 
    %     plot(regions(b).Centroid(1), regions(b).Centroid(2), 'r*');
    %     bbox = regions(b).BoundingBox;
    %     rectangle('Position', [bbox(1), bbox(2), bbox(3), bbox(4)], 'EdgeColor', 'g', 'LineWidth', 1);
    %     text(regions(b).Centroid(1)+5, regions(b).Centroid(2)+5,num2str(b),'Color', [1 0 0])
    % end
    % hold off;
    % title('Detected Regions with Centroids and Bounding Boxes');



    %%  Motion correction for whole movie
    [shift,dfdx,dfdy, ~] = mini_ali_motion_corr_info(movie, params.motion_reference_frames(1):params.motion_reference_frames(2));
    shift_movie = shift(:,1:params.DMD_off_frames(1)-1);
    
    %% Spike Detection for each region
    disp('Spike Detection for each Region')
    
    spk_COM = cell(size(regions,1),1);% record all spike COMs in x,y,frame
    non_spk_region_id = cell(size(regions,1),1);% record region id for non-spiking regions
    non_spk_COM = cell(size(regions,1),1);% record all COMs for non-spiking regions in x,y,frame
    non_spk_trace = cell(size(regions,1),1);% record traces for non-spiking regions
    non_spk_mask = cell(size(regions,1),1);% record masks for non-spiking regions
    
    
    for c = 1:size(regions,1)
        % Report the progress
        disp(['Spike Detection in no.' num2str(c) ' region of ' num2str(size(regions,1)) ' regions']);
       
        % load patch movie
        %[patch_coordinate, patch_movie, patch_mask] = mini_ali_patch_movie_load(movie(:,:,1:params.DMD_off_frames(1)-1), regions(c), patch_size);
        [patch_coordinate, patch_movie, patch_mask] = mini_ali_patch_movie_load(movie, regions(c), patch_size);
        patch_movie = patch_movie(:,:,1:params.DMD_off_frames(1)-1);
       
        % Detect spike TIMES using only 1-D trace
        [peak_time,region_trace] = mini_ali_patch_detect_spike_times(patch_movie, patch_coordinate, patch_mask, dfdx, dfdy, shift_movie, median_filter, threshold_factor, params.min_spike_interval_frames);
    
        % skip the loop if no spike detected
        if isempty(peak_time)
            non_spk_region_id{c} = uint32(c);
            non_spk_COM{c} = regions(c).Centroid;
            non_spk_trace{c} = region_trace;
            non_spk_mask{c} = uint32(regions(c).PixelIdxList);
        continue
        end
    
        % Convert patch coordinates to global ranges
        patch_rows = patch_coordinate(1):patch_coordinate(1)+size(patch_movie,1)-1;
        patch_cols = patch_coordinate(2):patch_coordinate(2)+size(patch_movie,2)-1;
    
        % Generate dF ONLY for candidate spike frames
        dF_AP = mini_ali_extract_dF_frames(movie, patch_rows, patch_cols, peak_time, median_filter, dfdx, dfdy, shift);
    
        % Coarse spatial position
        spk_coarse = mini_ali_spatial_spike_coarse(dF_AP, patch_mask, peak_time);
    
        % SVD denoising
        dF_AP = mini_ali_SVD_denoise(dF_AP, nsvd);
    
        % Fine localization
        spk_fine = mini_ali_spk_fine(dF_AP .* patch_mask, params.fine_npix, params.fine_radius_px, spk_coarse);
        spk_fine(any(isnan(spk_fine),2),:) = [];
    
        % exclude spk_fine points that are localized outside of region pixels
        spk_fine_pixel_idx = round(spk_fine(:,1:2)); 
        spk_fine_pixel_idx = sub2ind(size(patch_mask),spk_fine_pixel_idx(:,1),spk_fine_pixel_idx(:,2));
        spk_fine_pixel_mask = patch_mask (spk_fine_pixel_idx);
        spk_exclude = spk_fine_pixel_mask == 0;
        spk_fine(spk_exclude,:) = [];
    
        % project back to the whole movie coordinate
        spk_fine (:,[1 2]) = spk_fine (:,[1 2]) + patch_coordinate - [1 1];
    
        spk_COM{c} = spk_fine;
    end
    
    % Merge spikes if close in space and same in time
    COMs = vertcat(spk_COM{:});
    COMs = mini_ali_mergeCloseSpikesKeepMax(COMs, params.spike_merge_distance_px);
    
    % Compact non-spiking candidates
    nonspk_valid = ~cellfun(@isempty, non_spk_region_id);    
    non_spk_region_id = vertcat(non_spk_region_id{nonspk_valid});   
    non_spk_COM = vertcat(non_spk_COM{nonspk_valid});   
    non_spk_trace = vertcat(non_spk_trace{nonspk_valid});  
    non_spk_mask = non_spk_mask(nonspk_valid);
    %% ALI spike cluster;
    [clust_idx,clust_cen,alimap] = mini_ali_ALI_clustering(COMs,params.ALI_cluster_factor,nrow,ncol,params.ALI_cluster_threshold);
    % visualize the spike clustering 
    %mini_ali_plot_cluster(reference_image, COMs, clust_idx);

    %% Extract footprint for each spike-cluster
    nclust = size(clust_cen,1);
    footprint_cell =  cell(1,nclust);% record all spike COMs in x,y,frame
    footprint_center_cell =  cell(1,nclust);% record all footprint center in x,y
    valid_count = 0;
    
    for i= 1:nclust
         disp(['Footprint extraction for No.' num2str(i) ' cluster of ' num2str(nclust) ' regions']);    
        % extract cluster spikes
        clust_spk = COMs(clust_idx == i,:);
        if isempty(clust_spk)
            continue
        end
        valid_count = valid_count + 1;
        clust_cent = mean(clust_spk(:,1:2),1);
        clust_cent_round = round(clust_cent);
        % Generate a footprint in neighbouring 21x21 pixels of peaks
        selection_x = ...
        max(clust_cent_round(1,1)-params.footprint_radius_px,1):...
        min(clust_cent_round(1,1)+params.footprint_radius_px,nrow);
        selection_y = ...
        max(clust_cent_round(1,2)-params.footprint_radius_px,1):...
        min(clust_cent_round(1,2)+params.footprint_radius_px,ncol);
        % Spike frames belonging to this neuron
        footprint_frame = round(clust_spk(:,3));
        % Generate ONLY required dF frames
        footprint_stack = mini_ali_extract_dF_frames(movie, selection_x, selection_y, footprint_frame, median_filter, dfdx, dfdy, shift);
        % Average spike images
        footprint = mean(footprint_stack, 3);
        footprint(footprint < 0) = 0;
        % Save footprint as index and value
        [row_idx, col_idx] = ndgrid(selection_x, selection_y);
        idx = sub2ind([nrow,ncol], row_idx(:), col_idx(:));
        values = footprint(:);
        footprint_info = {uint32(idx), single(values)};
        footprint_cell{valid_count} =footprint_info;
        footprint_center_cell{valid_count} = clust_cent;
    end
    footprint_cell = footprint_cell(1:valid_count);
    footprint_center_cell = footprint_center_cell(1:valid_count);
    footprint = footprint_cell;
    footprint_center = vertcat(footprint_center_cell{:});
    %% Find pseudoinverse and left multiply to get traces for each footprint
    cell_traces = mini_ali_trace_generator(footprint,movie,Background,dfdx, dfdy, shift);

    %% Build master cell table for this FOV
    
    % Generate Cell ID
    nCells = size(cell_traces, 1);
    local_cell_id = uint16((1:nCells)');
    cell_id = strings(nCells,1);
    for c = 1:nCells
        cell_id(c) = sprintf('%s_C%04d', fov_id, c);
    end
    
    % Coordinates
    y_px = single(footprint_center(:,1));
    x_px = single(footprint_center(:,2));
    trace_row = uint32((1:nCells)');
    
    % Metadata columns
    experiment_col = repmat(experiment_id,nCells,1);
    plate_col = repmat(uint16(plate_number),nCells,1);
    well_col = repmat(well_name,nCells,1);
    fov_col = repmat(uint16(fov_number),nCells,1);
    fov_id_col = repmat(fov_id,nCells,1);
    fov_index_col = repmat(fov_index,nCells,1);
    
    % Cell table
    
    cells = table( ...
        cell_id, ...
        experiment_col, ...
        plate_col, ...
        well_col, ...
        fov_col, ...
        fov_id_col, ...
        fov_index_col, ...
        local_cell_id, ...
        x_px, ...
        y_px, ...
        trace_row);
    
    cells.Properties.VariableNames = { ...
        'cell_id', ...
        'experiment_id', ...
        'plate', ...
        'well', ...
        'fov', ...
        'fov_id', ...
        'fov_index', ...
        'local_cell_id', ...
        'x_px', ...
        'y_px', ...
        'trace_row'};

    %% save raw traces in "signals.h5"
    signals_file = fullfile(sub_analysis_folder, "signals.h5");
    
    if exist(signals_file,'file')
        delete(signals_file);
    end
    
    if nCells > 0
        trace_chunk_cells = min(256,nCells);
        trace_chunk_frames = min(2048,size(cell_traces,2));
        h5create( ...
            signals_file, ...
            '/traces/raw', ...
            size(cell_traces), ...
            'Datatype','single', ...
            'ChunkSize',[ ...
                trace_chunk_cells ...
                trace_chunk_frames], ...
            'Deflate',4);
    
        h5write( ...
            signals_file, ...
            '/traces/raw', ...
            cell_traces);
    
    else
        warning( ...
            'No spiking cells found in %s', ...
            fov_id);
    end
    
    %% save the cells table separately
    parquetwrite( ...
        fullfile(sub_analysis_folder,"cells.parquet"), ...
        cells);
    
    save( ...
        fullfile(sub_analysis_folder,"footprints.mat"), ...
        "footprint", ...
        "-v7.3");
    
    reference_image = single(reference_image);
    save( ...
        fullfile(sub_analysis_folder,"reference_image.mat"), ...
        "reference_image", ...
        "-v7.3");
    save( ...
        fullfile(sub_analysis_folder,"fov_metadata.mat"), ...
        "fov_metadata");
    
    %% save all non-spiking candidiate
    nNonSpk = size(non_spk_COM,1);
    
    if nNonSpk > 0
    
        candidate_id = strings(nNonSpk,1);
    
        for c = 1:nNonSpk
            candidate_id(c) = sprintf('%s_NS%04d', fov_id, c);
        end
    
        nonspk_y = single(non_spk_COM(:,2));
        nonspk_x = single(non_spk_COM(:,1));
    
        nonspiking_candidates = table( ...
            candidate_id, ...
            non_spk_region_id, ...
            repmat(uint16(plate_number),nNonSpk,1), ...
            repmat(well_name,nNonSpk,1), ...
            repmat(uint16(fov_number),nNonSpk,1), ...
            repmat(fov_id,nNonSpk,1), ...
            nonspk_x, ...
            nonspk_y, ...
        'VariableNames', { ...
            'candidate_id', ...
            'region_id', ...
            'plate', ...
            'well', ...
            'fov', ...
            'fov_id', ...
            'x_px', ...
            'y_px'});
    
        parquetwrite( ...
            fullfile(sub_analysis_folder,...
            "nonspiking_candidates.parquet"), ...
            nonspiking_candidates);
    
    end
    % save non-spiking traces
    if nNonSpk > 0
    
        non_spk_trace = single(non_spk_trace);
    
        save( ...
            fullfile( ...
                sub_analysis_folder, ...
                "nonspiking_data.mat"), ...
            "non_spk_trace", ...
            "non_spk_mask", ...
            "-v7.3");
    end

    %% QC Report Generation
    % Generate QC report sub folder
    qc_folder = fullfile(sub_analysis_folder, "QC");
    if ~exist(qc_folder,'dir')
        mkdir(qc_folder);
    end
    % Generate QC report sub folder
    figs = mini_ali_generate_QC_plots( ...
    reference_image, ...
    Binary_mask, ...
    regions, ...
    COMs, ...
    clust_idx, ...
    footprint_center, ...
    cell_traces(:,1:end-20), ...
    fps, ...
    qc_folder, ...
    fov_id, ...
    "qc");

    % Plot individual QC matlab figure
    % figs = mini_ali_generate_QC_plots( ...
    % reference_image, ...
    % Binary_mask, ...
    % regions, ...
    % COMs, ...
    % clust_idx, ...
    % footprint_center, ...
    % cell_traces, ...
    % fps, ...
    % qc_folder, ...
    % fov_id, ...
    % "traces", ...
    % 'Visible',"on", ...
    % 'Save',false);
end



