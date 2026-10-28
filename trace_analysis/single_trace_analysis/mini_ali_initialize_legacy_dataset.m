function data = mini_ali_initialize_legacy_dataset(raw_traces, protocol, varargin)
% MINI_ALI_INITIALIZE_LEGACY_DATASET
% Convert a legacy raw-trace matrix into the data structure expected by
% the Optopatch Trace Analysis V2 pipeline.
%
% INPUT
% raw_traces
%       Ncell x Nframe raw fluorescence trace matrix
% protocol
%       experiment protocol structure
% OPTIONAL
%
% 'DatasetID'
%       identifier used to generate synthetic cell IDs
%       default:
%           "LEGACY"
% 'StrictFrameCount'
%       require trace length == protocol.nFrames
%       default:
%           true
%
% OUTPUT
% data.cells
%       minimal cell metadata table
% data.traces
%       raw trace matrix
% data.fov_index
%       empty table for compatibility
% data.loaded_fovs
%       empty table for compatibility
%
% IMPORTANT
% Unknown legacy metadata are intentionally stored as:
%       string metadata -> ""
%       numeric metadata -> NaN
% No biological/spatial information is invented.


%% =========================================================
% Parse inputs
% ==========================================================

p = inputParser;
addParameter( p, 'DatasetID', "LEGACY");
addParameter( p, 'StrictFrameCount', true);
parse(p,varargin{:});
opt = p.Results;
dataset_id = string(opt.DatasetID);

%% =========================================================
% Validate raw traces
% ==========================================================
if isempty(raw_traces)
    error('raw_traces is empty.');
end
if ~ismatrix(raw_traces)
    error('raw_traces must be a 2-D matrix: Ncell x Nframe.');
end
nCells = size(raw_traces,1);
nFrames = size(raw_traces,2);
fprintf('Legacy raw trace matrix: %d cells x %d frames\n', nCells, nFrames);

%% =========================================================
% Validate frame count
% ==========================================================
if opt.StrictFrameCount
    if nFrames ~= protocol.nFrames
        error( ...
            ['Legacy trace length does not match protocol.\n' ...
             'Legacy traces: %d frames\n' ...
             'Protocol: %d frames'], ...
            nFrames, ...
            protocol.nFrames);
    end

end


%% =========================================================
% Convert traces to storage-efficient format
% ==========================================================
raw_traces = single(raw_traces);

%% =========================================================
% Generate minimum cell table
% ==========================================================
cells = table;

%% ---------------------------------------------------------
% Permanent synthetic cell identity
% ---------------------------------------------------------
prefix = repmat(dataset_id, nCells, 1);
cells.cell_id = compose("%s_C%06d", prefix, (1:nCells)');

%% ---------------------------------------------------------
% Row mapping
% This is the critical field:
% cells.analysis_row(i)
% corresponds to:
% raw_traces(i,:)
% ---------------------------------------------------------
cells.analysis_row = uint32((1:nCells)');

%% Keep original legacy trace row explicitly
cells.legacy_trace_row = uint32((1:nCells)');

%% ---------------------------------------------------------
% Source information
% ---------------------------------------------------------
cells.source_type = repmat("legacy_raw_trace_matrix", nCells, 1);

%% =========================================================
% Unknown spatial / experimental metadata
% ==========================================================
% These are intentionally left empty.
% Do NOT invent plate/well/FOV information.
cells.plate = repmat("",nCells,1);
cells.well = repmat("",nCells,1);
cells.fov_id = repmat("",nCells,1);
cells.fov_name = repmat("",nCells,1);

%% Unknown coordinates
cells.x_px = nan(nCells,1,'single');
cells.y_px = nan(nCells,1,'single');
%% Unknown original local cell number
cells.local_cell_id = nan(nCells,1,'single');

%% =========================================================
% Final alignment validation
% ==========================================================
assert(height(cells) == size(raw_traces,1), ...
    'Legacy cell table does not match trace matrix.');
assert( ...
    isequal( ...
        double(cells.analysis_row), ...
        (1:nCells)'), ...
    'analysis_row is incorrectly initialized.');

%% =========================================================
% Output using same interface as mini_ali_load_dataset
% ==========================================================
data = struct;
data.cells = cells;
data.traces = raw_traces;

%% These are unavailable for legacy matrix input
data.fov_index = table;
data.loaded_fovs = table;
%% Input information
data.source_type = "legacy_raw_trace_matrix";
data.dataset_id = dataset_id;
fprintf( ...
    'Legacy dataset initialized successfully.\n');
fprintf( ...
    '%d synthetic cell IDs generated.\n', ...
    nCells);
end