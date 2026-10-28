%% RUN_POPULATION_ANALYSIS_EXAMPLE
% Layer-2 population analysis using outputs from analysis_manuscript_V2.
%
% Required files in analysis_root:
%   analysis_v2.mat   -> features, cells_analysis, protocol, analysis_config
%   spike_data.mat    -> spike_data
%
% Outputs are written to:
%   <analysis_root>/PopulationAnalysis/
%popu
% Choose ONE of the two population modes below.

clearvars -except analysis_root

%% ========================================================================
% 1. Point to one completed V2 analysis folder
% =========================================================================
% Example:
% analysis_root = "D:\Data\Experiment_001\Analysis\TraceAnalysis";

if ~exist('analysis_root','var') || strlength(string(analysis_root)) == 0
    error('Set analysis_root before running this script.');
end

analysis_root = string(analysis_root);

analysis_file = fullfile(analysis_root,"analysis_v2.mat");
spike_file = fullfile(analysis_root,"spike_data.mat");

if ~isfile(analysis_file)
    error('Cannot find %s',analysis_file);
end
if ~isfile(spike_file)
    error('Cannot find %s',spike_file);
end

A = load(analysis_file, ...
    'features','cells_analysis','protocol','analysis_config');
S = load(spike_file,'spike_data');

features = A.features;
cells_analysis = A.cells_analysis;
protocol = A.protocol;
analysis_config = A.analysis_config;
spike_data = S.spike_data;

%% ========================================================================
% 2. Layer-2 configuration
% =========================================================================
population_config = make_population_analysis_config();

% Layer-2 QC filter.
% Only cells with cells_analysis.qc_pass == true contribute to Layer 2.
% Silent/non-spiking cells are still retained if qc_pass preserves them.
population_config.use_qc_pass = true;
population_config.qc.pass_variable = "qc_pass";

% Reuse the SAME onset settings and peak location used in Layer 0.
population_config.onset = analysis_config.onset;
population_config.STA.peak_idx = ...
    analysis_config.features.spike_extract_pre_frames + 1;

% Fixed time-bin definition.
population_config.STA.bin_frames = 100;

%% ========================================================================
% 3. Choose population assignment mode
% =========================================================================
mode = "all_cells";   % "all_cells" or "cell_id"

switch mode
    %% --------------------------------------------------------------------
    % MODE A: whole analysis result is one arrayed perturbation
    % ---------------------------------------------------------------------
    case "all_cells"

        population_definition.mode = "all_cells";

        % EDIT THIS LABEL for the current arrayed perturbation.
        population_definition.group_id = "KCNQ2_KD";

    %% --------------------------------------------------------------------
    % MODE B: sequencing/decoding assigns each cell_id to a gene/group
    % ---------------------------------------------------------------------
    case "cell_id"

        population_definition.mode = "cell_id";

        % Example mapping table columns:
        %   cell_id    gene
        %   cell_001   KCNQ2
        %   cell_002   SCN2A
        %   cell_003   control
        %
        % Replace with your sequencing-decoding result.
        mapping_file = fullfile(analysis_root,"cell_to_gene_mapping.csv");

        if ~isfile(mapping_file)
            error('Cannot find mapping file: %s',mapping_file);
        end

        population_definition.mapping = readtable(mapping_file);

        % Any mapping column can be used, e.g. "gene", "sgRNA",
        % "perturbation", or a combined perturbation label.
        population_definition.group_variable = "gene";

    otherwise
        error('Unknown mode: %s',mode);
end

%% ========================================================================
% 4. Build Layer-2 population phenotypes
% =========================================================================
[group_features, ...
 group_STA, ...
 group_STA_features, ...
 group_info, ...
 cell_STA_features, ...
 population_qc] = ...
    build_population_features( ...
        features, ...
        spike_data, ...
        protocol, ...
        population_definition, ...
        population_config, ...
        cells_analysis);

%% ========================================================================
% 5. Save Layer-2 outputs
% =========================================================================
population_root = fullfile(analysis_root,"PopulationAnalysis");
if ~isfolder(population_root)
    mkdir(population_root);
end

parquetwrite( ...
    fullfile(population_root,"group_features.parquet"), ...
    group_features);

parquetwrite( ...
    fullfile(population_root,"group_STA_features.parquet"), ...
    group_STA_features);

parquetwrite( ...
    fullfile(population_root,"group_info.parquet"), ...
    group_info);


% Save population QC
population_qc_table = struct2table( ...
    population_qc, ...
    'AsArray', true);

parquetwrite( ...
    fullfile(population_root,"population_qc.parquet"), ...
    population_qc_table);








if ~isempty(cell_STA_features)
    parquetwrite( ...
        fullfile(population_root,"cell_STA_features.parquet"), ...
        cell_STA_features);
end

save( ...
    fullfile(population_root,"group_STA.mat"), ...
    'group_STA', ...
    'population_definition', ...
    'population_config', ...
    '-v7.3');

save( ...
    fullfile(population_root,"population_analysis.mat"), ...
    'group_features', ...
    'group_STA_features', ...
    'group_info', ...
    'population_qc', ...
    'population_definition', ...
    'population_config', ...
    '-v7.3');

fprintf('\nPopulation QC summary:\n');
fprintf('  Total cells:        %d\n', population_qc.n_feature_cells);
fprintf('  Assigned cells:     %d\n', population_qc.n_cells_assigned);
fprintf('  Unassigned cells:   %d\n', population_qc.n_cells_unassigned);

if population_qc.use_qc_pass
    fprintf('  QC-excluded cells:  %d\n', population_qc.n_cells_excluded_qc);
    fprintf('  QC-pass assigned:   %d\n', ...
        population_qc.n_cells_assigned - population_qc.n_cells_excluded_qc);
end

fprintf('  Included in Layer2: %d\n', population_qc.n_cells_selected);
fprintf('  With spike data:    %d\n', population_qc.n_cells_selected_with_spike_data);
fprintf('  Without spike data: %d\n', population_qc.n_cells_selected_without_spike_data);
fprintf('\nPopulation outputs saved to:\n%s\n',population_root);
