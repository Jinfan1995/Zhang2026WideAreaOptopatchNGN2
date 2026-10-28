function population_config = make_population_analysis_config()
% MAKE_POPULATION_ANALYSIS_CONFIG
% Configuration for Layer-2 population analysis.
%
% This function is intentionally independent of protocol so it can be
% created before protocol loading. build_population_features receives the
% protocol separately and uses protocol.fps / protocol.nFrames directly.

%% ========================================================================
% 1. Population summary statistics
% =========================================================================
population_config.summary.include_quantiles = true;
population_config.summary.exclude_feature_names = [ ...
    "cell_id", ...
    "cell_index", ...
    "analysis_row"];

%% ========================================================================
% 2. Optional cell-level QC filtering
% =========================================================================
% If true, build_population_features requires cells_analysis with:
%   cells_analysis.cell_id
%   cells_analysis.qc_pass
%
% Silent/non-spiking cells should remain qc_pass when preprocessing is
% valid; their AP-count features remain valid biological zeros.
% Recommended default for Layer 2:
% only cells that passed the finalized Layer-1 QC contribute to population
% summaries and population STA analysis.
population_config.use_qc_pass = true;

% The pass variable expected in cells_analysis.
population_config.qc.pass_variable = "qc_pass";

% Missing cell_id/QC records are treated as failed and excluded.
population_config.qc.exclude_missing_qc = true;

%% ========================================================================
% 3. Time-binned STA analysis
% =========================================================================
% Use the same 100-frame definition as the Layer-1 AP-count trajectory.
population_config.STA.bin_frames = 100;

% Build one robust cell-level STA per cell/bin using the median across APs.
population_config.STA.cell_STA_method = "median";

% Minimum APs required to create a cell-bin STA. Keep at 1 so sparse bins
% are retained; n_cells and n_spikes are stored for downstream filtering.
population_config.STA.min_spikes_per_cell_bin = 1;

% Population STA is the mean across CELL STAs, not across pooled APs.
population_config.STA.population_STA_method = "mean";

% Minimum number of cell STAs required to characterize the population mean
% waveform. Keep at 1 by default; downstream plotting/statistics can impose
% a stricter n threshold.
population_config.STA.min_cells_for_population_STA = 1;

% When true, save the cell-bin STA trace in cell_STA_features. This can make
% the table very large and is therefore false by default. The traces needed
% for the population output are always retained internally and the group
% mean/std/sem waveforms are saved in group_STA.mat.
population_config.STA.store_cell_STA_trace_in_table = false;

%% ========================================================================
% 4. Progress reporting
% =========================================================================
population_config.progress.enable = true;
population_config.progress.percent_interval = 10;

end
