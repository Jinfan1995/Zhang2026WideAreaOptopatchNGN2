function protocol = make_protocol_double_ramp_800Hz(saveFolder)
% MAKE_PROTOCOL_DOUBLE_RAMP_800HZ
%
% Generate the stimulus/analysis protocol structure for the
% double-ramp Optopatch experiment acquired at 800 Hz.
%
% Usage:
%   protocol = make_protocol_double_ramp_800Hz();
%
% Or save protocol to disk:
%   protocol = make_protocol_double_ramp_800Hz('D:\Optopatch\protocols');
%
% Output:
%   protocol : MATLAB structure containing acquisition, stimulus,
%              baseline, binning, and analysis information.
%
% -------------------------------------------------------------

%% ============================================================
% 1. Protocol identity
% =============================================================

protocol.name        = "double_ramp_800Hz";
protocol.version     = "1.1";
protocol.description = "Optopatch double-ramp stimulation protocol";

protocol.created_date = datetime('today');


%% ============================================================
% 2. Acquisition parameters
% =============================================================

protocol.fps = 800;

% Number of frames used by the standard electrophysiology analysis.
% Keep protocol.nFrames as the analysis length because downstream functions
% use this field when defining bins and validating epochs.
protocol.nFrames = 8399;

protocol.duration_sec = protocol.nFrames / protocol.fps;


%% ============================================================
% 2B. Waveform / raw-acquisition layout
% =============================================================
%
% These fields describe how the acquired raw waveform is partitioned.
% They are protocol facts rather than analysis-algorithm settings.
%
% For the current waveform:
%   frames 1:8399   = standard electrophysiology analysis
%   frames 8400:8408 = VMAP absolute-voltage reference segment

protocol.waveform.raw_nFrames = 8408;

protocol.waveform.analysis_frame_range = [1 8399];

protocol.waveform.analysis_nFrames = ...
    diff(protocol.waveform.analysis_frame_range) + 1;

protocol.waveform.raw_duration_sec = ...
    protocol.waveform.raw_nFrames / protocol.fps;


%% ============================================================
% 2C. Absolute-voltage / VMAP measurement
% =============================================================
%
% VMAP_ratio is calculated from RAW traces before exponential correction or
% normalization:
%
%   VMAP_ratio =
%       median(raw numerator window) / median(raw denominator window)
%
% Store frame windows as [start end] ranges so a future waveform only needs
% this protocol section changed.

protocol.VMAP.enabled = true;

protocol.VMAP.feature_name = "VMAP_ratio";

protocol.VMAP.source = "raw_traces";

protocol.VMAP.method = "median_ratio";

protocol.VMAP.numerator_frame_range = [8200 8399];

protocol.VMAP.denominator_frame_range = [8400 8408];


%% ============================================================
% 3. General analysis parameters
% =============================================================

% Median filter used in trace preprocessing
protocol.analysis.median_filter_size = 80;

% Number of frames used for temporal bins
protocol.analysis.bin_frames = 200;

% Equivalent time in seconds
protocol.analysis.bin_sec = protocol.analysis.bin_frames / protocol.fps;


%% ============================================================
% 4. Analysis bin locations
% =============================================================

% This replaces: pins = 1:200:8399;

protocol.analysis.bin_starts = ...
    1 : protocol.analysis.bin_frames : protocol.nFrames;

% More useful if bins are defined as:
%
% frameStart = bin_edges(i)
% frameEnd   = bin_edges(i+1)-1
%
% protocol.analysis.bin_edges = ...
%     [protocol.analysis.bin_starts, protocol.nFrames + 1];
% 

%% ============================================================
% 5. Baseline frames
% =============================================================

% Keep individual baseline periods separately.
% This makes the protocol much easier to inspect later.

protocol.baseline.period(1).name   = "baseline_1";
protocol.baseline.period(1).frames = 50:799;

protocol.baseline.period(2).name   = "baseline_2";
protocol.baseline.period(2).frames = 1550:1599;

protocol.baseline.period(3).name   = "baseline_3";
protocol.baseline.period(3).frames = 2350:2399;

protocol.baseline.period(4).name   = "baseline_4";
protocol.baseline.period(4).frames = 3150:3199;

protocol.baseline.period(5).name   = "baseline_5";
protocol.baseline.period(5).frames = 3950:3999;

protocol.baseline.period(6).name   = "baseline_6";
protocol.baseline.period(6).frames = 5950:5999;

protocol.baseline.period(7).name   = "baseline_7";
protocol.baseline.period(7).frames = 8200:8379;


%% Combine all baseline periods
protocol.baseline.frames = ...
    unique([protocol.baseline.period.frames]);

% Also store baseline time in seconds
protocol.baseline.time_sec = ...
    (protocol.baseline.frames - 1) / protocol.fps;


%% ============================================================
% 6. Stimulus summary
% =============================================================
% Stimulus calibration
protocol.stimulus.device = "blue_LED";
protocol.stimulus.control = "analog_voltage";
protocol.stimulus.full_scale_voltage_V = 5; % voltage corresponding to 100%
protocol.stimulus.full_scale_percent = 100;
protocol.stimulus.unit = "%";

protocol.epoch(1).name = "step_1";
protocol.epoch(1).type = "step";
protocol.epoch(1).start_frame = 800;
protocol.epoch(1).end_frame = 1199;
protocol.epoch(1).stimulus_percent = 25;

protocol.epoch(2).name = "step_2";
protocol.epoch(2).type = "step";
protocol.epoch(2).start_frame = 1600;
protocol.epoch(2).end_frame = 1999;
protocol.epoch(2).stimulus_percent = 55;

protocol.epoch(3).name = "step_3";
protocol.epoch(3).type = "step";
protocol.epoch(3).start_frame = 2400;
protocol.epoch(3).end_frame = 2799;
protocol.epoch(3).stimulus_percent = 75;

protocol.epoch(4).name = "step_4";
protocol.epoch(4).type = "step";
protocol.epoch(4).start_frame = 3200;
protocol.epoch(4).end_frame = 3599;
protocol.epoch(4).stimulus_percent = 100;

protocol.epoch(5).name = "ramp_1";
protocol.epoch(5).type = "ramp";
protocol.epoch(5).start_frame = 4000;
protocol.epoch(5).end_frame = 5599;
protocol.epoch(5).stimulus_percent = 100/3;


protocol.epoch(6).name = "ramp_2";
protocol.epoch(6).type = "ramp";
protocol.epoch(6).start_frame = 6000;
protocol.epoch(6).end_frame = 7599;
protocol.epoch(6).stimulus_percent = 100;
% ============================================================
% Blue-off epochs following each stimulus
% =============================================================
protocol.epoch(7).name = "blue_off_1";
protocol.epoch(7).type = "blue_off";
protocol.epoch(7).start_frame = 1200;
protocol.epoch(7).end_frame = 1599;
protocol.epoch(7).stimulus_percent = 0;

protocol.epoch(8).name = "blue_off_2";
protocol.epoch(8).type = "blue_off";
protocol.epoch(8).start_frame = 2000;
protocol.epoch(8).end_frame = 2399;
protocol.epoch(8).stimulus_percent = 0;

protocol.epoch(9).name = "blue_off_3";
protocol.epoch(9).type = "blue_off";
protocol.epoch(9).start_frame = 2800;
protocol.epoch(9).end_frame = 3199;
protocol.epoch(9).stimulus_percent = 0;

protocol.epoch(10).name = "blue_off_4";
protocol.epoch(10).type = "blue_off";
protocol.epoch(10).start_frame = 3600;
protocol.epoch(10).end_frame = 3999;
protocol.epoch(10).stimulus_percent = 0;

protocol.epoch(11).name = "blue_off_5";
protocol.epoch(11).type = "blue_off";
protocol.epoch(11).start_frame = 5600;
protocol.epoch(11).end_frame = 5999;
protocol.epoch(11).stimulus_percent = 0;

protocol.epoch(12).name = "blue_off_6";
protocol.epoch(12).type = "blue_off";
protocol.epoch(12).start_frame = 7600;
protocol.epoch(12).end_frame = 8399;
protocol.epoch(12).stimulus_percent = 0;

for i = 1:numel(protocol.epoch)
    protocol.epoch(i).stimulus_voltage_V = ...
        protocol.epoch(i).stimulus_percent / 100 * ...
        protocol.stimulus.full_scale_voltage_V;
end

%% ============================================================
% 7. Protocol validation
% =============================================================

% Make sure baseline frames are valid
assert(all(protocol.baseline.frames >= 1), ...
    'Protocol contains frame numbers < 1.');

assert(all(protocol.baseline.frames <= protocol.nFrames), ...
    'Protocol contains frames larger than nFrames.');

% Make sure bins are valid
assert(protocol.analysis.bin_frames > 0, ...
    'bin_frames must be > 0.');

% Make sure FPS is valid
assert(protocol.fps > 0, ...
    'fps must be > 0.');

% Make sure the analysis waveform agrees with protocol.nFrames
assert(protocol.waveform.analysis_frame_range(1) >= 1, ...
    'Waveform analysis start frame must be >= 1.');

assert(protocol.waveform.analysis_frame_range(2) == protocol.nFrames, ...
    ['protocol.waveform.analysis_frame_range must end at ' ...
     'protocol.nFrames.']);

assert(protocol.waveform.raw_nFrames >= protocol.nFrames, ...
    'Raw waveform cannot be shorter than the analysis waveform.');

% Validate VMAP windows when enabled
if protocol.VMAP.enabled
    assert(numel(protocol.VMAP.numerator_frame_range) == 2, ...
        'VMAP numerator_frame_range must be [start end].');

    assert(numel(protocol.VMAP.denominator_frame_range) == 2, ...
        'VMAP denominator_frame_range must be [start end].');

    assert(protocol.VMAP.numerator_frame_range(1) >= 1 && ...
           protocol.VMAP.numerator_frame_range(2) >= ...
           protocol.VMAP.numerator_frame_range(1), ...
        'Invalid VMAP numerator frame range.');

    assert(protocol.VMAP.denominator_frame_range(1) >= 1 && ...
           protocol.VMAP.denominator_frame_range(2) >= ...
           protocol.VMAP.denominator_frame_range(1), ...
        'Invalid VMAP denominator frame range.');

    assert(protocol.VMAP.denominator_frame_range(2) <= ...
           protocol.waveform.raw_nFrames, ...
        ['VMAP denominator window exceeds the configured raw waveform ' ...
         'length.']);
end


%% ============================================================
% 8. Save protocol
% =============================================================

if nargin >= 1 && ~isempty(saveFolder)

    if ~exist(saveFolder,'dir')
        mkdir(saveFolder);
    end

    fileName = sprintf( ...
        'protocol_%s_v%s.mat', ...
        protocol.name, ...
        protocol.version);

    savePath = fullfile(saveFolder, fileName);

    save(savePath, 'protocol');

    fprintf('Protocol saved:\n%s\n', savePath);

end

end