function [group_features, group_STA, group_STA_features, group_info, cell_STA_features, population_qc] = ...
    build_population_features( ...
        features, spike_data, protocol, population_definition, ...
        population_config, cells_analysis)
% BUILD_POPULATION_FEATURES
% Build Layer-2 population/perturbation phenotypes from:
%   Layer 1: features table (one row per cell)
%   Layer 0: spike_data struct (individual AP events)
%
% Two population-definition modes are supported:
%
%   1) population_definition.mode = "cell_id"
%      Use a mapping table, e.g.
%          population_definition.mapping = sequencing_mapping;
%          population_definition.group_variable = "gene";
%      Required mapping columns:
%          cell_id
%          <group_variable>
%
%   2) population_definition.mode = "all_cells"
%      Treat the entire input dataset as one arrayed perturbation, e.g.
%          population_definition.group_id = "KCNQ2_KD";
%
% OUTPUT
% -------------------------------------------------------------------------
% group_features
%   Long-format Layer-2 summaries of every numeric/logical Layer-1 feature:
%       group_id, feature_name, mean, std, sem, median, q25, q75, n
%
% group_STA
%   Struct containing 100-frame-bin population STA waveforms. Individual APs
%   are NOT pooled directly. The hierarchy is:
%       APs -> median STA within cell/bin -> mean across cell STAs.
%
% group_STA_features
%   Long-format summary of STA-shape features for each group/bin. Includes
%   both:
%       population_STA_value  = feature(mean population waveform)
%       mean/std/sem/...      = distribution across per-cell STAs
%
% group_info
%   Group-level counts.
%
% cell_STA_features
%   Layer-1.5 table: one row per cell/bin with a valid cell STA. This is
%   useful for downstream statistics and preserves equal cell weighting.
%
% population_qc
%   Struct summarizing population assignment, QC filtering, spike-data
%   availability, and STA record availability for this Layer-2 run.
%
% IMPORTANT
% -------------------------------------------------------------------------
% - cell_id is the join key between features, spike_data, mapping, and
%   optional cells_analysis.
% - 0 remains a valid biological zero; NaN means unavailable measurement.
% - n is feature-specific and counts finite values only.
% -------------------------------------------------------------------------

if nargin < 5 || isempty(population_config)
    population_config = make_population_analysis_config();
end

if nargin < 6
    cells_analysis = table();
end

%% ========================================================================
% 1. Validate core inputs
% =========================================================================
if ~istable(features)
    error('features must be a MATLAB table.');
end

if ~ismember('cell_id',features.Properties.VariableNames)
    error('features must contain cell_id.');
end

if ~isstruct(spike_data)
    error('spike_data must be a struct array.');
end

if ~isfield(protocol,'fps') || ~isfinite(protocol.fps) || protocol.fps <= 0
    error('protocol.fps is missing or invalid.');
end

if isfield(protocol,'nFrames') && isfinite(protocol.nFrames)
    nFrames = double(protocol.nFrames);
else
    % Fall back to latest epoch endpoint if nFrames is not explicitly saved.
    if ~isfield(protocol,'epoch') || isempty(protocol.epoch)
        error('Need protocol.nFrames or protocol.epoch to determine trace length.');
    end
    nFrames = max(double([protocol.epoch.end_frame]));
end

fps = double(protocol.fps);
nFeatureCells = height(features);
feature_cell_keys = local_id_keys(features.cell_id);

% Require unique Layer-1 cell identity.
if numel(unique(feature_cell_keys)) ~= nFeatureCells
    error('features.cell_id contains duplicate cell IDs.');
end

% Map spike_data by its permanent cell_id.
nSpikeCells = numel(spike_data);
spike_cell_keys = strings(nSpikeCells,1);
for iCell = 1:nSpikeCells
    if isfield(spike_data(iCell),'cell_id') && ~isempty(spike_data(iCell).cell_id)
        spike_cell_keys(iCell) = local_scalar_id_key(spike_data(iCell).cell_id);
    end
end

nonempty_spike_keys = spike_cell_keys ~= "";
if numel(unique(spike_cell_keys(nonempty_spike_keys))) ~= nnz(nonempty_spike_keys)
    error('spike_data contains duplicate non-empty cell_id values.');
end

[has_spike_data, spike_row_for_feature] = ...
    ismember(feature_cell_keys,spike_cell_keys);

%% ========================================================================
% 2. Resolve population assignment
% =========================================================================
if ~isfield(population_definition,'mode')
    error('population_definition.mode is required.');
end

mode = lower(string(population_definition.mode));
selected = false(nFeatureCells,1);
group_label = strings(nFeatureCells,1);
input_group_label = strings(nFeatureCells,1);

switch mode
    case "all_cells"
        if ~isfield(population_definition,'group_id') || ...
                strlength(string(population_definition.group_id)) == 0
            error('all_cells mode requires population_definition.group_id.');
        end

        selected(:) = true;
        group_label(:) = string(population_definition.group_id);
        input_group_label = group_label;

    case "cell_id"
        if ~isfield(population_definition,'mapping') || ...
                ~istable(population_definition.mapping)
            error('cell_id mode requires population_definition.mapping table.');
        end

        mapping = population_definition.mapping;
        if ~ismember('cell_id',mapping.Properties.VariableNames)
            error('population_definition.mapping must contain cell_id.');
        end

        if isfield(population_definition,'group_variable')
            group_variable = char(string(population_definition.group_variable));
        else
            group_variable = 'group_id';
        end

        if ~ismember(group_variable,mapping.Properties.VariableNames)
            error('Mapping table does not contain group variable "%s".',group_variable);
        end

        mapping_keys = local_id_keys(mapping.cell_id);
        if numel(unique(mapping_keys)) ~= numel(mapping_keys)
            error(['population_definition.mapping contains duplicate cell_id values. ' ...
                   'Each cell should map to one population group.']);
        end

        mapping_groups = string(mapping.(group_variable));
        valid_mapping = mapping_keys ~= "" & ~ismissing(mapping_groups) & mapping_groups ~= "";

        mapping_keys = mapping_keys(valid_mapping);
        mapping_groups = mapping_groups(valid_mapping);

        [selected,mapping_row] = ismember(feature_cell_keys,mapping_keys);
        group_label(selected) = mapping_groups(mapping_row(selected));
        input_group_label = group_label;

    otherwise
        error('Unknown population_definition.mode: %s',mode);
end

%% ========================================================================
% 3. Optional qc_pass filtering
% =========================================================================
use_qc_pass = false;
if isfield(population_config,'use_qc_pass')
    use_qc_pass = logical(population_config.use_qc_pass);
end

qc_pass_for_feature = true(nFeatureCells,1);

if use_qc_pass
    if isempty(cells_analysis) || ~istable(cells_analysis)
        error('population_config.use_qc_pass=true requires cells_analysis table.');
    end

    required_qc_vars = {'cell_id','qc_pass'};
    if ~all(ismember(required_qc_vars,cells_analysis.Properties.VariableNames))
        error('cells_analysis must contain cell_id and qc_pass.');
    end

    qc_keys = local_id_keys(cells_analysis.cell_id);
    if numel(unique(qc_keys)) ~= numel(qc_keys)
        error('cells_analysis.cell_id contains duplicate values.');
    end

    [has_qc,qc_row] = ismember(feature_cell_keys,qc_keys);
    qc_pass_for_feature(:) = false;
    qc_pass_for_feature(has_qc) = logical(cells_analysis.qc_pass(qc_row(has_qc)));

    selected = selected & qc_pass_for_feature;
end

selected = selected & group_label ~= "" & ~ismissing(group_label);

if ~any(selected)
    error('No cells remain after population assignment/QC filtering.');
end

groups = unique(group_label(selected),'stable');
nGroups = numel(groups);

%% ========================================================================
% 4. Group info
% =========================================================================
group_info = table( ...
    groups(:), ...
    zeros(nGroups,1), ...
    zeros(nGroups,1), ...
    zeros(nGroups,1), ...
    'VariableNames', ...
    {'group_id','n_cells_input','n_cells_analyzed','n_cells_with_spike_data'});

for iGroup = 1:nGroups
    g = groups(iGroup);
    group_info.n_cells_input(iGroup) = nnz(input_group_label == g);
    group_mask = selected & group_label == g;
    group_info.n_cells_analyzed(iGroup) = nnz(group_mask);
    group_info.n_cells_with_spike_data(iGroup) = nnz(group_mask & has_spike_data);
end

%% ========================================================================
% 5. Layer-1 feature population summaries
% =========================================================================
feature_names = string(features.Properties.VariableNames);
exclude_names = ["cell_id","cell_index","analysis_row"];
if isfield(population_config,'summary') && ...
        isfield(population_config.summary,'exclude_feature_names')
    exclude_names = unique([ ...
        exclude_names, ...
        string(population_config.summary.exclude_feature_names)]);
end

is_summary_feature = false(size(feature_names));
for iVar = 1:numel(feature_names)
    var_name = char(feature_names(iVar));
    if any(feature_names(iVar) == exclude_names)
        continue
    end
    x = features.(var_name);
    if isnumeric(x) || islogical(x)
        if isvector(x) && numel(x) == nFeatureCells
            is_summary_feature(iVar) = true;
        end
    end
end

summary_feature_names = feature_names(is_summary_feature);
nSummaryFeatures = numel(summary_feature_names);
nSummaryRows = nGroups * nSummaryFeatures;

group_id_col = strings(nSummaryRows,1);
feature_name_col = strings(nSummaryRows,1);
mean_col = nan(nSummaryRows,1);
std_col = nan(nSummaryRows,1);
sem_col = nan(nSummaryRows,1);
median_col = nan(nSummaryRows,1);
q25_col = nan(nSummaryRows,1);
q75_col = nan(nSummaryRows,1);
n_col = zeros(nSummaryRows,1);

row = 0;
for iGroup = 1:nGroups
    group_mask = selected & group_label == groups(iGroup);

    for iVar = 1:nSummaryFeatures
        row = row + 1;
        var_name = char(summary_feature_names(iVar));
        feature_vector = features.(var_name);
        x = double(feature_vector(group_mask));
        S = local_summary_stats(x);

        group_id_col(row) = groups(iGroup);
        feature_name_col(row) = summary_feature_names(iVar);
        mean_col(row) = S.mean;
        std_col(row) = S.std;
        sem_col(row) = S.sem;
        median_col(row) = S.median;
        q25_col(row) = S.q25;
        q75_col(row) = S.q75;
        n_col(row) = S.n;
    end
end

group_features = table( ...
    group_id_col,feature_name_col, ...
    mean_col,std_col,sem_col,median_col,q25_col,q75_col,n_col, ...
    'VariableNames', ...
    {'group_id','feature_name','mean','std','sem','median','q25','q75','n'});

%% ========================================================================
% 6. Define fixed 100-frame STA bins
% =========================================================================
default_bin_frames = 100;
if isfield(population_config,'STA') && ...
        isfield(population_config.STA,'bin_frames')
    bin_frames = double(population_config.STA.bin_frames);
else
    bin_frames = default_bin_frames;
end

if ~isscalar(bin_frames) || ~isfinite(bin_frames) || bin_frames < 1
    error('population_config.STA.bin_frames must be a positive scalar.');
end
bin_frames = round(bin_frames);

bin_start = (1:bin_frames:nFrames)';
bin_end = min(bin_start + bin_frames - 1,nFrames);
nBins = numel(bin_start);

min_spikes_per_cell_bin = 1;
if isfield(population_config,'STA') && ...
        isfield(population_config.STA,'min_spikes_per_cell_bin')
    min_spikes_per_cell_bin = ...
        max(1,round(double(population_config.STA.min_spikes_per_cell_bin)));
end

cell_STA_method = "median";
if isfield(population_config,'STA') && ...
        isfield(population_config.STA,'cell_STA_method')
    cell_STA_method = lower(string(population_config.STA.cell_STA_method));
end

population_STA_method = "mean";
if isfield(population_config,'STA') && ...
        isfield(population_config.STA,'population_STA_method')
    population_STA_method = lower(string(population_config.STA.population_STA_method));
end

min_cells_for_population_STA = 1;
if isfield(population_config,'STA') && ...
        isfield(population_config.STA,'min_cells_for_population_STA')
    min_cells_for_population_STA = ...
        max(1,round(double(population_config.STA.min_cells_for_population_STA)));
end

store_cell_STA_trace = false;
if isfield(population_config,'STA') && ...
        isfield(population_config.STA,'store_cell_STA_trace_in_table')
    store_cell_STA_trace = ...
        logical(population_config.STA.store_cell_STA_trace_in_table);
end

% Peak position in each Layer-0 peak-aligned waveform.
if isfield(population_config,'STA') && isfield(population_config.STA,'peak_idx')
    peak_idx = round(double(population_config.STA.peak_idx));
elseif isfield(protocol,'analysis') && isfield(protocol.analysis,'spike_extract_pre_frames')
    peak_idx = round(double(protocol.analysis.spike_extract_pre_frames)) + 1;
elseif isfield(population_definition,'peak_idx')
    peak_idx = round(double(population_definition.peak_idx));
else
    % Normal/current pipeline stores this in analysis_config, not protocol.
    % We infer it from the waveform width later if needed; current default
    % extraction uses 20 pre-peak frames.
    peak_idx = 21;
end

% Prefer explicitly supplied onset config.
if isfield(population_config,'onset')
    onset_cfg = population_config.onset;
elseif isfield(population_definition,'onset')
    onset_cfg = population_definition.onset;
else
    error(['Layer-2 STA characterization requires onset settings. ' ...
           'Set population_config.onset = analysis_config.onset before calling.']);
end

%% ========================================================================
% 7. Build one median STA per cell / 100-frame bin (Layer 1.5)
% =========================================================================
selected_rows = find(selected);
nSelected = numel(selected_rows);

progress_enable = false;
progress_interval = inf;
if isfield(population_config,'progress') && ...
        isfield(population_config.progress,'enable') && ...
        population_config.progress.enable
    progress_enable = true;
    percent_interval = 10;
    if isfield(population_config.progress,'percent_interval')
        percent_interval = double(population_config.progress.percent_interval);
    end
    progress_interval = max(1,round(nSelected * percent_interval / 100));
end

% Internal dynamic records. Number of populated cell/bin combinations is
% not known in advance and is usually much smaller than nCells*nBins.
record_cell_id = strings(0,1);
record_group_id = strings(0,1);
record_bin = zeros(0,1);
record_start = zeros(0,1);
record_end = zeros(0,1);
record_n_spikes = zeros(0,1);
record_trace = cell(0,1);
record_available = false(0,1);
record_shape = nan(0,20);

% STA shape features are ordered by threshold level:
%   20% -> 50% -> 80%
% Within each threshold block:
%   total width (ms), rise duration (ms), fall duration (ms),
%   rise duration (frames), fall duration (frames).
% Non-threshold waveform features follow afterward.
shape_feature_names = [ ...
    "STA_width20_ms", ...
    "STA_width20_rise_ms", ...
    "STA_width20_fall_ms", ...
    "STA_width20_rise_frames", ...
    "STA_width20_fall_frames", ...
    "STA_width50_ms", ...
    "STA_width50_rise_ms", ...
    "STA_width50_fall_ms", ...
    "STA_width50_rise_frames", ...
    "STA_width50_fall_frames", ...
    "STA_width80_ms", ...
    "STA_width80_rise_ms", ...
    "STA_width80_fall_ms", ...
    "STA_width80_rise_frames", ...
    "STA_width80_fall_frames", ...
    "STA_height", ...
    "STA_AHP", ...
    "STA_AHP_delay_ms", ...
    "STA_upstroke", ...
    "STA_downstroke"];

fprintf('\n============================================================\n');
fprintf('Building Layer-2 population features\n');
fprintf('============================================================\n');
fprintf('Selected cells: %d\n',nSelected);
fprintf('Groups:         %d\n',nGroups);
fprintf('STA bins:       %d (%d frames/bin)\n',nBins,bin_frames);
fprintf('------------------------------------------------------------\n');

t_start = tic;

for iSelected = 1:nSelected
    feature_row = selected_rows(iSelected);

    if ~has_spike_data(feature_row)
        continue
    end

    iSpikeCell = spike_row_for_feature(feature_row);

    if ~isfield(spike_data(iSpikeCell),'spike_frames') || ...
            ~isfield(spike_data(iSpikeCell),'peak_aligned_traces')
        continue
    end

    spike_frames = double(spike_data(iSpikeCell).spike_frames(:));
    waveforms = double(spike_data(iSpikeCell).peak_aligned_traces);

    if isempty(spike_frames) || isempty(waveforms)
        continue
    end

    if size(waveforms,1) ~= numel(spike_frames)
        warning('build_population_features:SpikeWaveformMismatch', ...
            'Skipping cell_id %s: spike_frames and peak_aligned_traces are not aligned.', ...
            feature_cell_keys(feature_row));
        continue
    end

    waveform_len = size(waveforms,2);
    if peak_idx < 1 || peak_idx > waveform_len
        error(['Configured peak_idx=%d is incompatible with Layer-0 waveform length=%d. ' ...
               'Set population_config.STA.peak_idx explicitly.'], ...
               peak_idx,waveform_len);
    end

    valid_event = ...
        isfinite(spike_frames) & ...
        spike_frames >= 1 & ...
        spike_frames <= nFrames & ...
        all(isfinite(waveforms),2);

    spike_frames = spike_frames(valid_event);
    waveforms = waveforms(valid_event,:);

    if isempty(spike_frames)
        continue
    end

    event_bin = floor((spike_frames-1) ./ bin_frames) + 1;
    event_bin = min(max(event_bin,1),nBins);
    populated_bins = unique(event_bin(:)');

    for iBinLocal = 1:numel(populated_bins)
        iBin = populated_bins(iBinLocal);
        in_bin = event_bin == iBin;
        nSpikesBin = nnz(in_bin);

        if nSpikesBin < min_spikes_per_cell_bin
            continue
        end

        W = waveforms(in_bin,:);

        switch cell_STA_method
            case "median"
                cell_STA = median(W,1,'omitnan');
            case "mean"
                cell_STA = mean(W,1,'omitnan');
            otherwise
                error('Unknown population_config.STA.cell_STA_method: %s',cell_STA_method);
        end

        F = characterize_STA_waveform(cell_STA,fps,peak_idx,onset_cfg);

        record_cell_id(end+1,1) = feature_cell_keys(feature_row); %#ok<AGROW>
        record_group_id(end+1,1) = group_label(feature_row); %#ok<AGROW>
        record_bin(end+1,1) = iBin; %#ok<AGROW>
        record_start(end+1,1) = bin_start(iBin); %#ok<AGROW>
        record_end(end+1,1) = bin_end(iBin); %#ok<AGROW>
        record_n_spikes(end+1,1) = nSpikesBin; %#ok<AGROW>
        record_trace{end+1,1} = single(cell_STA); %#ok<AGROW>
        record_available(end+1,1) = logical(F.available); %#ok<AGROW>
        record_shape(end+1,:) = [ ... %#ok<AGROW>
            F.width20_ms, ...
            F.width20_rise_ms, ...
            F.width20_fall_ms, ...
            F.width20_rise_frames, ...
            F.width20_fall_frames, ...
            F.width50_ms, ...
            F.width50_rise_ms, ...
            F.width50_fall_ms, ...
            F.width50_rise_frames, ...
            F.width50_fall_frames, ...
            F.width80_ms, ...
            F.width80_rise_ms, ...
            F.width80_fall_ms, ...
            F.width80_rise_frames, ...
            F.width80_fall_frames, ...
            F.height, ...
            F.AHP, ...
            F.AHP_delay_ms, ...
            F.upstroke, ...
            F.downstroke];
    end

    if progress_enable && ...
            (mod(iSelected,progress_interval) == 0 || iSelected == nSelected)
        fprintf('Cell STA extraction: %5.1f%% (%d/%d)\n', ...
            100*iSelected/nSelected,iSelected,nSelected);
    end
end

%% ========================================================================
% 8. Build Layer-1.5 cell_STA_features table
% =========================================================================
nRecords = numel(record_bin);

if nRecords == 0
    cell_STA_features = table();
else
    % cell_id is emitted as a stable string join key in this Layer-1.5 table.
    cell_STA_features = table( ...
        record_cell_id, ...
        record_group_id, ...
        record_bin, ...
        record_start, ...
        record_end, ...
        record_n_spikes, ...
        record_available, ...
        record_shape(:,1), ...
        record_shape(:,2), ...
        record_shape(:,3), ...
        record_shape(:,4), ...
        record_shape(:,5), ...
        record_shape(:,6), ...
        record_shape(:,7), ...
        record_shape(:,8), ...
        record_shape(:,9), ...
        record_shape(:,10), ...
        record_shape(:,11), ...
        record_shape(:,12), ...
        record_shape(:,13), ...
        record_shape(:,14), ...
        record_shape(:,15), ...
        record_shape(:,16), ...
        record_shape(:,17), ...
        record_shape(:,18), ...
        record_shape(:,19), ...
        record_shape(:,20), ...
        'VariableNames',{ ...
            'cell_id','group_id','bin_number','start_frame','end_frame', ...
            'n_spikes','STA_available', ...
            'STA_width20_ms', ...
            'STA_width20_rise_ms','STA_width20_fall_ms', ...
            'STA_width20_rise_frames','STA_width20_fall_frames', ...
            'STA_width50_ms', ...
            'STA_width50_rise_ms','STA_width50_fall_ms', ...
            'STA_width50_rise_frames','STA_width50_fall_frames', ...
            'STA_width80_ms', ...
            'STA_width80_rise_ms','STA_width80_fall_ms', ...
            'STA_width80_rise_frames','STA_width80_fall_frames', ...
            'STA_height','STA_AHP','STA_AHP_delay_ms', ...
            'STA_upstroke','STA_downstroke'});

    if store_cell_STA_trace
        cell_STA_features.STA_trace = record_trace;
    end
end

%% ========================================================================
% 9. Population STA waveforms and population/cell-STA feature summaries
% =========================================================================
empty_bin = struct( ...
    'bin_number',NaN, ...
    'start_frame',NaN, ...
    'end_frame',NaN, ...
    'n_cells',0, ...
    'n_spikes_total',0, ...
    'mean_trace',single([]), ...
    'std_trace',single([]), ...
    'sem_trace',single([]), ...
    'population_features',struct());

empty_group = struct( ...
    'group_id',"", ...
    'n_cells_analyzed',0, ...
    'bin',repmat(empty_bin,nBins,1));

group_STA = repmat(empty_group,nGroups,1);

nShapeFeatures = numel(shape_feature_names);
nSTAFeatureRows = nGroups * nBins * nShapeFeatures;

sta_group_col = strings(nSTAFeatureRows,1);
sta_bin_col = zeros(nSTAFeatureRows,1);
sta_start_col = zeros(nSTAFeatureRows,1);
sta_end_col = zeros(nSTAFeatureRows,1);
sta_feature_col = strings(nSTAFeatureRows,1);
sta_population_value_col = nan(nSTAFeatureRows,1);
sta_mean_col = nan(nSTAFeatureRows,1);
sta_std_col = nan(nSTAFeatureRows,1);
sta_sem_col = nan(nSTAFeatureRows,1);
sta_median_col = nan(nSTAFeatureRows,1);
sta_q25_col = nan(nSTAFeatureRows,1);
sta_q75_col = nan(nSTAFeatureRows,1);
sta_n_col = zeros(nSTAFeatureRows,1);
sta_n_cells_col = zeros(nSTAFeatureRows,1);
sta_n_spikes_col = zeros(nSTAFeatureRows,1);

sta_row = 0;

for iGroup = 1:nGroups
    g = groups(iGroup);
    group_STA(iGroup).group_id = g;
    group_STA(iGroup).n_cells_analyzed = ...
        group_info.n_cells_analyzed(iGroup);

    for iBin = 1:nBins
        record_idx = find(record_group_id == g & record_bin == iBin);
        nCellsSTA = numel(record_idx);
        nSpikesTotal = sum(record_n_spikes(record_idx));

        group_STA(iGroup).bin(iBin).bin_number = iBin;
        group_STA(iGroup).bin(iBin).start_frame = bin_start(iBin);
        group_STA(iGroup).bin(iBin).end_frame = bin_end(iBin);
        group_STA(iGroup).bin(iBin).n_cells = nCellsSTA;
        group_STA(iGroup).bin(iBin).n_spikes_total = nSpikesTotal;

        population_shape_values = nan(1,nShapeFeatures);

        if nCellsSTA > 0
            trace_matrix = vertcat(record_trace{record_idx});

            switch population_STA_method
                case "mean"
                    population_trace = mean(trace_matrix,1,'omitnan');
                case "median"
                    population_trace = median(trace_matrix,1,'omitnan');
                otherwise
                    error('Unknown population_config.STA.population_STA_method: %s', ...
                        population_STA_method);
            end

            trace_std = std(double(trace_matrix),0,1,'omitnan');
            trace_sem = trace_std ./ sqrt(nCellsSTA);

            group_STA(iGroup).bin(iBin).mean_trace = single(population_trace);
            group_STA(iGroup).bin(iBin).std_trace = single(trace_std);
            group_STA(iGroup).bin(iBin).sem_trace = single(trace_sem);

            if nCellsSTA >= min_cells_for_population_STA
                PF = characterize_STA_waveform( ...
                    population_trace,fps,peak_idx,onset_cfg);

                group_STA(iGroup).bin(iBin).population_features = PF;
                population_shape_values = [ ...
                    PF.width20_ms, ...
                    PF.width20_rise_ms, ...
                    PF.width20_fall_ms, ...
                    PF.width20_rise_frames, ...
                    PF.width20_fall_frames, ...
                    PF.width50_ms, ...
                    PF.width50_rise_ms, ...
                    PF.width50_fall_ms, ...
                    PF.width50_rise_frames, ...
                    PF.width50_fall_frames, ...
                    PF.width80_ms, ...
                    PF.width80_rise_ms, ...
                    PF.width80_fall_ms, ...
                    PF.width80_rise_frames, ...
                    PF.width80_fall_frames, ...
                    PF.height, ...
                    PF.AHP, ...
                    PF.AHP_delay_ms, ...
                    PF.upstroke, ...
                    PF.downstroke];
            end
        else
            group_STA(iGroup).bin(iBin).mean_trace = single([]);
            group_STA(iGroup).bin(iBin).std_trace = single([]);
            group_STA(iGroup).bin(iBin).sem_trace = single([]);
            group_STA(iGroup).bin(iBin).population_features = struct();
        end

        for iShape = 1:nShapeFeatures
            sta_row = sta_row + 1;

            if nCellsSTA > 0
                x = record_shape(record_idx,iShape);
            else
                x = [];
            end
            S = local_summary_stats(x);

            sta_group_col(sta_row) = g;
            sta_bin_col(sta_row) = iBin;
            sta_start_col(sta_row) = bin_start(iBin);
            sta_end_col(sta_row) = bin_end(iBin);
            sta_feature_col(sta_row) = shape_feature_names(iShape);
            sta_population_value_col(sta_row) = population_shape_values(iShape);
            sta_mean_col(sta_row) = S.mean;
            sta_std_col(sta_row) = S.std;
            sta_sem_col(sta_row) = S.sem;
            sta_median_col(sta_row) = S.median;
            sta_q25_col(sta_row) = S.q25;
            sta_q75_col(sta_row) = S.q75;
            sta_n_col(sta_row) = S.n;
            sta_n_cells_col(sta_row) = nCellsSTA;
            sta_n_spikes_col(sta_row) = nSpikesTotal;
        end
    end
end

group_STA_features = table( ...
    sta_group_col,sta_bin_col,sta_start_col,sta_end_col,sta_feature_col, ...
    sta_population_value_col, ...
    sta_mean_col,sta_std_col,sta_sem_col,sta_median_col,sta_q25_col,sta_q75_col, ...
    sta_n_col,sta_n_cells_col,sta_n_spikes_col, ...
    'VariableNames',{ ...
        'group_id','bin_number','start_frame','end_frame','feature_name', ...
        'population_STA_value', ...
        'mean','std','sem','median','q25','q75','n', ...
        'n_cells_with_STA','n_spikes_total'});

%% ========================================================================
% 10. Population-level QC summary
% =========================================================================
assigned_mask = input_group_label ~= "" & ~ismissing(input_group_label);

population_qc = struct;
population_qc.mode = mode;
population_qc.use_qc_pass = use_qc_pass;
population_qc.n_feature_cells = nFeatureCells;
population_qc.n_spike_data_cells = nSpikeCells;
population_qc.n_cells_assigned = nnz(assigned_mask);
population_qc.n_cells_unassigned = nFeatureCells - nnz(assigned_mask);
population_qc.n_cells_selected = nnz(selected);

if use_qc_pass
    population_qc.n_cells_excluded_qc = ...
        nnz(assigned_mask & ~qc_pass_for_feature);
else
    population_qc.n_cells_excluded_qc = 0;
end

population_qc.n_cells_selected_with_spike_data = ...
    nnz(selected & has_spike_data);
population_qc.n_cells_selected_without_spike_data = ...
    nnz(selected & ~has_spike_data);
population_qc.n_groups = nGroups;
population_qc.n_STA_bins = nBins;
population_qc.n_cell_bin_STA_records = nRecords;
population_qc.n_available_cell_bin_STA = nnz(record_available);
population_qc.n_unavailable_cell_bin_STA = ...
    nRecords - nnz(record_available);
population_qc.n_total_spikes_in_STA_records = sum(record_n_spikes);
population_qc.created_datetime = datetime("now");

fprintf('------------------------------------------------------------\n');
fprintf('Layer-2 population analysis complete\n');
fprintf('Cell-bin STAs:  %d\n',nRecords);
fprintf('Elapsed time:   %.2f min\n',toc(t_start)/60);
fprintf('============================================================\n');

end


%% =========================================================================
% LOCAL HELPERS
% =========================================================================
function keys = local_id_keys(x)
% Convert a table column of IDs into stable string join keys.
if iscell(x)
    keys = strings(numel(x),1);
    for i = 1:numel(x)
        keys(i) = local_scalar_id_key(x{i});
    end
else
    keys = string(x(:));
end
keys = strip(keys);
end


function key = local_scalar_id_key(x)
if isstring(x)
    if isempty(x)
        key = "";
    else
        key = string(x(1));
    end
elseif ischar(x)
    key = string(x);
elseif isnumeric(x) || islogical(x) || iscategorical(x)
    if isempty(x)
        key = "";
    else
        key = string(x(1));
    end
else
    try
        key = string(x);
        if numel(key) > 1
            key = key(1);
        end
    catch
        error('Unsupported cell_id type: %s',class(x));
    end
end
key = strip(key);
end


function S = local_summary_stats(x)
x = double(x(:));
x = x(isfinite(x));

S = struct('mean',NaN,'std',NaN,'sem',NaN, ...
    'median',NaN,'q25',NaN,'q75',NaN,'n',0);

S.n = numel(x);
if S.n == 0
    return
end

S.mean = mean(x);
S.median = median(x);
S.q25 = local_quantile(x,0.25);
S.q75 = local_quantile(x,0.75);

if S.n >= 2
    S.std = std(x,0);
    S.sem = S.std ./ sqrt(S.n);
else
    S.std = NaN;
    S.sem = NaN;
end
end


function q = local_quantile(x,p)
% Toolbox-independent linear-interpolated quantile.
x = sort(double(x(:)));
n = numel(x);

if n == 0
    q = NaN;
    return
elseif n == 1
    q = x(1);
    return
end

p = min(max(double(p),0),1);
pos = 1 + (n-1)*p;
lo = floor(pos);
hi = ceil(pos);

if lo == hi
    q = x(lo);
else
    frac = pos-lo;
    q = x(lo) + frac*(x(hi)-x(lo));
end
end

