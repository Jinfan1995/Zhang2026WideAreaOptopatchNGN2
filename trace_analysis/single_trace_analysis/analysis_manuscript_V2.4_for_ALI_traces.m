%% =========================================================
% Optopatch Trace Analysis V2.4 - Legacy Trace Input + residual drift QC
% Input:
%   legacy trace matrix stored in a MAT file
% Output:
%   one row per cell with permanent cell_id
% ==========================================================

%% =========================================================
% 1. Experiment
% ==========================================================
root_path = uigetdir();
if root_path == 0
    error("No experiment folder selected.");
end
analysis_root = fullfile( ...
    root_path, ...
    "Analysis", ...
    "TraceAnalysis");
if ~exist(analysis_root,'dir')
    mkdir(analysis_root);
end

%% =========================================================
% 2. Load experimental protocol
% ==========================================================
protocol = make_protocol_double_ramp_800Hz();

%% =========================================================
% 3. Load analysis configuration
% ==========================================================
analysis_config = make_trace_analysis_config_v2();

%% =========================================================
% 4. Load / initialize legacy trace dataset
% ==========================================================
%
% USER SETTINGS FOR LEGACY DATA
%
% The selected root_path should contain this MAT file.
legacy_file = fullfile(root_path, "cell_profile.mat");

% Name of the trace matrix variable inside the MAT file.
legacy_trace_variable = "ALI_traces";

% Permanent dataset identifier used when creating cell_id values.
% Change this for each independent legacy dataset.
legacy_dataset_id = "JZ043";

% Allow legacy matrices whose frame count does not exactly match
% protocol.nFrames. Protocol validation below still requires all
% explicitly referenced baseline/stimulus frames to exist.
legacy_strict_frame_count = false;

if ~isfile(legacy_file)
    error("Legacy MAT file not found: %s", legacy_file);
end

legacy_data = load(legacy_file);

if ~isfield(legacy_data, legacy_trace_variable)
    error( ...
        'Variable "%s" was not found in %s.', ...
        legacy_trace_variable, ...
        legacy_file);
end

legacy_raw_traces = legacy_data.(legacy_trace_variable);

if ~isnumeric(legacy_raw_traces) || ndims(legacy_raw_traces) ~= 2
    error( ...
        'Legacy trace variable "%s" must be a numeric 2-D matrix.', ...
        legacy_trace_variable);
end

data = ...
    mini_ali_initialize_legacy_dataset( ...
        legacy_raw_traces, ...
        protocol, ...
        'DatasetID', legacy_dataset_id, ...
        'StrictFrameCount', legacy_strict_frame_count);

cells = data.cells;
raw_traces = data.traces;

fprintf("Loaded legacy dataset: %s\n", legacy_dataset_id);
fprintf("Loaded %d cells\n", height(cells));
fprintf( ...
    "Trace size: %d cells x %d frames\n", ...
    size(raw_traces,1), ...
    size(raw_traces,2));
%% =========================================================
% 5. Validate protocol
% ==========================================================
nCells = height(cells);
nFrames = size(raw_traces,2);
if size(raw_traces,1) ~= nCells
    error("Number of cells and traces do not match.");
end

if max(protocol.baseline.frames) > nFrames
    error("Baseline frame exceeds trace length.");
end
fprintf("Protocol validation passed.\n");

%% =========================================================
% 6. Initialize QC table
% ==========================================================
qc = table;
qc.cell_id = cells.cell_id;
% Preprocessing
qc.exp_fit_fail = false(nCells,1);
qc.exp_fit_rmse = nan(nCells,1);
qc.normalization_fail = false(nCells,1);
qc.normalization_baseline = nan(nCells,1);

% Residual baseline drift
qc.residual_drift_fail = false(nCells,1);
qc.residual_drift_range = nan(nCells,1);
qc.residual_drift_slope_per_sec = nan(nCells,1);
qc.residual_drift_r2 = nan(nCells,1);

% SNR
qc.snr_fail = false(nCells,1);
qc.low_step_SNR = false(nCells,1);
qc.low_spike_SNR = false(nCells,1);
qc.abnormal_SNR = false(nCells,1);
qc.snr_low_step_fail = false(nCells,1);
qc.snr_low_spike_fail = false(nCells,1);
qc.snr_step_gt_spike_fail = false(nCells,1);
qc.snr_extreme_step_fail = false(nCells,1);
qc.snr_extreme_spike_fail = false(nCells,1);
qc.snr_extreme_noise_fail = false(nCells,1);
qc.snr_invalid_fail = false(nCells,1);
qc.snr_step = nan(nCells,1);
qc.snr_spike = nan(nCells,1);
qc.noise_std = nan(nCells,1);
% STA
qc.no_STA_spike = false(nCells,1);
qc.STA_kernel_available = false(nCells,1);
qc.STA_initial_spike_count = zeros(nCells,1,'uint16');
qc.STA_valid_spike_count = zeros(nCells,1,'uint16');
qc.STA_onset_available = false(nCells,1);
qc.STA_onset_to_peak_frames = nan(nCells,1);
qc.STA_onset_used_fallback = false(nCells,1);
% Feature extraction
qc.analysis_fail = false(nCells,1);
qc.mixed_trace = false(nCells,1);
qc.no_detected_spike = false(nCells,1);
qc.spike_shape_unavailable = false(nCells,1);
qc.STA_onset_fail = false(nCells,1);
qc.individual_onset_fail_count = zeros(nCells,1);

%% =========================================================
% 7. Exponential correction
% ==========================================================
[traces_corrected, exp_info] = ...
    traces_exp_corr( ...
        raw_traces, ...
        protocol.baseline.frames, ...
        protocol.analysis.median_filter_size, ...
        analysis_config.exp_fit);

qc.exp_fit_fail = exp_info.fail;
qc.exp_fit_rmse = exp_info.rmse;
%% =========================================================
% 8. Normalize traces
% ==========================================================
[traces_nm_pre_drift, norm_info] = ...
    traces_flip_norm( ...
        traces_corrected, ...
        protocol.baseline.frames);

qc.normalization_fail = norm_info.fail;
qc.normalization_baseline = norm_info.baseline;

%% =========================================================
% 9. Residual baseline-drift QC and correction
% ==========================================================

[traces_nm, drift_info] = ...
    traces_residual_drift_corr( ...
        traces_nm_pre_drift, ...
        protocol, ...
        analysis_config.drift);

qc.residual_drift_fail = drift_info.fail;
qc.residual_drift_range = drift_info.baseline_range;
qc.residual_drift_slope_per_sec = drift_info.slope_per_sec;
qc.residual_drift_r2 = drift_info.r2;

%% =========================================================
% 10. SNR QC
% ==========================================================
snr_info = ...
    traces_SNR_filter( ...
        traces_nm, ...
        protocol.baseline.frames, ...
        protocol.analysis.median_filter_size, ...
        analysis_config.snr.step_threshold, ...
        analysis_config.snr.spike_threshold, ...
        analysis_config.snr);

% Overall SNR QC
qc.snr_fail = snr_info.fail;
% SNR measurements
qc.snr_step = snr_info.step;
qc.snr_spike = snr_info.spike;
qc.noise_std = snr_info.noise_std;
% Individual SNR failure reasons
qc.snr_low_step_fail  = snr_info.fail_low_step;
qc.snr_low_spike_fail  = snr_info.fail_low_spike;
qc.snr_step_gt_spike_fail = snr_info.fail_step_gt_spike;
qc.snr_extreme_step_fail = snr_info.fail_extreme_step;
qc.snr_extreme_spike_fail = snr_info.fail_extreme_spike;
qc.snr_extreme_noise_fail = snr_info.fail_extreme_noise;
qc.snr_invalid_fail = snr_info.fail_invalid;

qc.abnormal_SNR = ...
    snr_info.fail_step_gt_spike | ...
    snr_info.fail_extreme_step | ...
    snr_info.fail_extreme_spike | ...
    snr_info.fail_extreme_noise | ...
    snr_info.fail_invalid;
%% =========================================================
% 11. Determine which traces can enter spike analysis
% ==========================================================
valid_preprocessing = ...
    ~qc.exp_fit_fail & ...
    ~qc.normalization_fail & ...
    ~qc.residual_drift_fail & ...
    ~qc.abnormal_SNR;
fprintf( ...
    "%d / %d cells have valid preprocessing.\n", ...
    nnz(valid_preprocessing), ...
    nCells);
%% =========================================================
% 12. Build STA kernels
% ==========================================================
[kernel, sta_info] = traces_STA_kernel(traces_nm, protocol, analysis_config, valid_preprocessing);
% Store QC
qc.no_STA_spike = sta_info.no_spike;
qc.STA_kernel_available = sta_info.kernel_available;
qc.STA_initial_spike_count = sta_info.initial_spike_count;
qc.STA_valid_spike_count = sta_info.valid_spike_count;
qc.STA_onset_available = sta_info.onset_available;
qc.STA_onset_to_peak_frames = sta_info.onset_to_peak_frames;
qc.STA_onset_used_fallback = sta_info.onset_used_fallback;

%% =========================================================
% 13. Feature extraction
% ==========================================================

[features, spike_data, feature_qc] = ...
    traces_extract_features( ...
        traces_nm, ...
        kernel, ...
        protocol, ...
        analysis_config);
%% =========================================================
% 14. Update QC after feature extraction
% ==========================================================

qc.analysis_fail = feature_qc.invalid_input;
qc.no_detected_spike = feature_qc.no_spike;
qc.mixed_trace = feature_qc.close_spikes;
qc.spike_shape_unavailable = feature_qc.spike_shape_unavailable;
qc.STA_onset_fail = feature_qc.STA_onset_fail;
qc.individual_onset_fail_count = feature_qc.individual_onset_fail_count;
qc.low_step_SNR = snr_info.fail_low_step;
qc.low_spike_SNR = snr_info.fail_low_spike;

%% =========================================================
% 15. Build cells_analysis
% ==========================================================
nCells = height(cells);
assert(size(traces_nm,1) == nCells, ...
    'traces_nm and cells have different numbers of cells.');
assert(height(features) == nCells, ...
    'features and cells have different numbers of cells.');
assert(numel(qc.exp_fit_fail) == nCells, ...
    'QC vectors are not aligned with cells.');

% Start from original cell identity / metadata table
cells_analysis = cells;
% ----------------------------------------------------------
% SNR measurements
% ----------------------------------------------------------
cells_analysis.SNR_step = single(snr_info.step);
cells_analysis.SNR_spike = single(snr_info.spike);
cells_analysis.noise_std = single(snr_info.noise_std);
% ----------------------------------------------------------
% Drift / normalization metrics
% ----------------------------------------------------------
cells_analysis.exp_fit_rmse = single(qc.exp_fit_rmse);
cells_analysis.normalization_baseline = single(qc.normalization_baseline);
cells_analysis.residual_drift_range = single(qc.residual_drift_range);
cells_analysis.residual_drift_slope_per_sec = single(qc.residual_drift_slope_per_sec);
cells_analysis.residual_drift_r2 = single(qc.residual_drift_r2);
cells_analysis.residual_drift_post_range = ...
    single(drift_info.post_correction_range);
% ----------------------------------------------------------
% Preprocessing QC
% ----------------------------------------------------------
cells_analysis.qc_exp_fit_fail = logical(qc.exp_fit_fail);
cells_analysis.qc_normalization_fail = logical(qc.normalization_fail);
cells_analysis.qc_residual_drift_fail = logical(qc.residual_drift_fail);
% ----------------------------------------------------------
% SNR QC
% ----------------------------------------------------------
cells_analysis.qc_low_step_SNR = logical(qc.low_step_SNR);
cells_analysis.qc_low_spike_SNR = logical(qc.low_spike_SNR);
cells_analysis.qc_extreme_step_SNR = logical(qc.snr_extreme_step_fail);
cells_analysis.qc_extreme_spike_SNR = logical(qc.snr_extreme_spike_fail);
cells_analysis.qc_extreme_noise_std = logical(qc.snr_extreme_noise_fail);
cells_analysis.qc_abnormal_SNR = logical(qc.abnormal_SNR);
% ----------------------------------------------------------
% STA QC
% ----------------------------------------------------------
cells_analysis.qc_no_STA_spike = logical(qc.no_STA_spike);
cells_analysis.STA_kernel_available = logical(qc.STA_kernel_available);
cells_analysis.STA_initial_spike_count = uint16(qc.STA_initial_spike_count);
cells_analysis.STA_valid_spike_count = uint16(qc.STA_valid_spike_count);
cells_analysis.STA_onset_available = logical(qc.STA_onset_available);
cells_analysis.STA_onset_to_peak_frames = single(qc.STA_onset_to_peak_frames);
cells_analysis.STA_onset_to_peak_ms = ...
    single(qc.STA_onset_to_peak_frames ./ protocol.fps * 1000);
cells_analysis.STA_onset_used_fallback = logical(qc.STA_onset_used_fallback);
% ---------------------------------------------------------
% Feature extraction QC
% ----------------------------------------------------------
cells_analysis.qc_analysis_fail = logical(qc.analysis_fail);
cells_analysis.qc_mixed_trace = logical(qc.mixed_trace);
cells_analysis.qc_no_detected_spike = logical(qc.no_detected_spike);
cells_analysis.qc_spike_shape_unavailable = logical(qc.spike_shape_unavailable);
cells_analysis.qc_STA_onset_fail = logical(qc.STA_onset_fail);
cells_analysis.individual_onset_fail_count = uint16(qc.individual_onset_fail_count);
% ----------------------------------------------------------
% Validate QC vector lengths
% ----------------------------------------------------------
nCells = height(cells_analysis);
assert(numel(qc.exp_fit_rmse) == nCells);
assert(numel(qc.normalization_baseline) == nCells);
assert(numel(qc.residual_drift_fail) == nCells);
assert(numel(qc.residual_drift_range) == nCells);
assert(numel(qc.STA_kernel_available) == nCells);
assert(numel(qc.STA_initial_spike_count) == nCells);
assert(numel(qc.STA_valid_spike_count) == nCells);
assert(isequal(double(cells_analysis.analysis_row(:)), (1:nCells)'), ...
    'cells_analysis.analysis_row is not correctly aligned.');
%% =========================================================
% 16. Define overall QC pass
% ==========================================================
% In this case, a biologically non-spiking neurons still pass qc
cells_analysis.qc_pass = ...
    ~cells_analysis.qc_exp_fit_fail & ...
    ~cells_analysis.qc_normalization_fail & ...
    ~cells_analysis.qc_residual_drift_fail & ...
    ~cells_analysis.qc_abnormal_SNR & ...
    ~cells_analysis.qc_analysis_fail;

%% =========================================================
% 17. Attach cell IDs to feature table
% ==========================================================
assert(height(features) == height(cells), ...
    'Feature table lost cell alignment.');
features.cell_id = cells.cell_id;
if ismember('cell_index', features.Properties.VariableNames)
    features.cell_index = [];
end
features = movevars(features, 'cell_id', 'Before', 1);

%% =========================================================
% 18. Finalize spike-level data
% ==========================================================
for iCell = 1:nCells
    spike_data(iCell).cell_id = cells.cell_id(iCell);
end

%% =========================================================
% 19. Save final analysis
% ==========================================================
% save feature matrix
save( ...
    fullfile(analysis_root, 'analysis_v2.mat'), ...
    'cells_analysis', ...
    'features', ...
    'qc', ...
    'protocol', ...
    'analysis_config', ...
    'sta_info', ...
    'drift_info', ...
    '-v7.3');
parquetwrite( ...
    fullfile( ...
        analysis_root, ...
        "cells_analysis.parquet"), ...
    cells_analysis);
parquetwrite( ...
    fullfile( ...
        analysis_root, ...
        "features.parquet"), ...
    features);
% save spike information structure
save( ...
    fullfile(analysis_root, 'spike_data.mat'), ...
    'spike_data', ...
    '-v7.3');

analysis_metadata = struct;
analysis_metadata.protocol = protocol;
analysis_metadata.analysis_config = analysis_config;
analysis_metadata.residual_drift = struct( ...
    'method',string(analysis_config.drift.method), ...
    'max_baseline_range',double(analysis_config.drift.max_baseline_range), ...
    'min_anchor_periods',double(analysis_config.drift.min_anchor_periods));
analysis_metadata.created_datetime = datetime("now");
save(fullfile(analysis_root, "analysis_metadata.mat"), "analysis_metadata");

%% =========================================================
% 20. Save processed traces to HDF5
% ==========================================================
if analysis_config.storage.save_processed_traces
    processed_file = fullfile(analysis_root, "processed_signals.h5");
    % Start with a clean file
    if isfile(processed_file)
        delete(processed_file);
    end    
    % Matrix dimensions
    trace_size = size(traces_nm);
    % HDF5 chunk size
    % This works well for reading subsets of cells and time ranges.
    chunk_size = [min(64, trace_size(1)), min(1024, trace_size(2))];
    % ----------------------------------------------------------
    % Exponential-corrected traces
    % ----------------------------------------------------------
    if analysis_config.storage.save_exp_corrected_traces
        h5create( ...
            processed_file, ...
            "/traces/exp_corrected", ...
            trace_size, ...
            'Datatype','single', ...
            'ChunkSize',chunk_size, ...
            'Deflate',4);
        h5write( ...
            processed_file, ...
            "/traces/exp_corrected", ...
            single(traces_corrected));
    end    
    % ----------------------------------------------------------
    % Final normalized traces used for analysis
    % (after residual linear baseline-drift correction)
    % ----------------------------------------------------------
    if analysis_config.storage.save_normalized_traces
        h5create( ...
            processed_file, ...
            "/traces/normalized", ...
            trace_size, ...
            'Datatype','single', ...
            'ChunkSize',chunk_size, ...
            'Deflate',4);
        h5write( ...
            processed_file, ...
            "/traces/normalized", ...
            single(traces_nm));
    end
    % -------------------------------------------------------
    % Metadata
    % -------------------------------------------------------
    h5writeatt(processed_file, "/", "analysis_version", char(analysis_config.version));
    h5writeatt(processed_file, "/", "protocol_name", char(protocol.name));
    h5writeatt(processed_file, "/", "fps", double(protocol.fps));
    h5writeatt(processed_file, "/", "n_cells", uint32(nCells));
    h5writeatt(processed_file, "/", "n_frames", uint32(nFrames));
    h5writeatt(processed_file, "/", "cell_metadata_file", "cells_analysis.parquet");
    h5writeatt(processed_file, "/", "feature_file", "features.parquet");
    h5writeatt(processed_file, "/", "residual_drift_corrected", 'true');
    h5writeatt(processed_file, "/", "residual_drift_max_range", ...
        double(analysis_config.drift.max_baseline_range));
    h5writeatt(processed_file, "/traces", "row_mapping", "Row i corresponds to cells_analysis.analysis_row == i"); 
end
