function [kernel, sta_info] = traces_STA_kernel(traces_nm, protocol, analysis_config, valid_preprocessing)
% TRACES_STA_KERNEL
% Generate one peak-aligned spike-triggered-average (STA) kernel per cell.
%
% V2.1 onset strategy:
%   1. Detect high-SNR spikes.
%   2. Extract each AP aligned to its detected AP peak.
%   3. Average the peak-aligned APs. No individual-spike onset alignment is
%      used for kernel construction, avoiding derivative-selection bias.
%   4. Detect AP onset ONCE on the averaged STA for QC/diagnostics using
%      detect_AP_onset().
%   5. Build the matched-filter kernel around the known peak position.
%
% No cells are deleted. If a kernel cannot be generated, a neutral delta
% kernel is retained so downstream indexing remains unchanged.

%% ========================================================================
% 1. Input information
% =========================================================================
[nCells, nFrames] = size(traces_nm);
valid_preprocessing = logical(valid_preprocessing(:));

if numel(valid_preprocessing) ~= nCells
    error('traces_STA_kernel:InvalidPreprocessingMask', ...
        'valid_preprocessing must contain one value per cell.');
end

%% ========================================================================
% 2. Load protocol / analysis parameters
% =========================================================================
baseline_frames = unique(protocol.baseline.frames(:));
baseline_frames = baseline_frames( ...
    baseline_frames >= 1 & baseline_frames <= nFrames);

if isempty(baseline_frames)
    error('traces_STA_kernel:NoBaselineFrames', ...
        'No valid baseline frames were supplied.');
end

cfg = analysis_config.sta;
onset_cfg = analysis_config.onset;

median_filter_size = protocol.analysis.median_filter_size;
spike_SNR_th = cfg.spike_SNR_threshold;
min_peak_distance = cfg.min_peak_distance_frames;
min_spike_interval = cfg.min_spike_interval_frames;
edge_start = cfg.edge_start_exclusion_frames;
edge_end = cfg.edge_end_exclusion_frames;
extract_pre = cfg.extract_pre_frames;
extract_post = cfg.extract_post_frames;
kernel_half_width = cfg.kernel_half_width_frames;

% All extracted spikes are peak-aligned by construction.
expected_peak_idx = extract_pre + 1;

if extract_pre < kernel_half_width || extract_post < kernel_half_width
    error('traces_STA_kernel:KernelWindowTooLarge', ...
        ['STA extraction window must contain at least kernel_half_width ' ...
         'frames on both sides of the AP peak.']);
end

%% ========================================================================
% 3. Kernel dimensions
% =========================================================================
kernel_length = 2 * kernel_half_width + 1;
kernel_center = kernel_half_width + 1;

%% ========================================================================
% 4. Initialize outputs
% =========================================================================
% Default = neutral delta kernel.
kernel = zeros(nCells, kernel_length, 'like', traces_nm);
kernel(:,kernel_center) = 1;

sta_info.initial_spike_count = zeros(nCells,1);
sta_info.valid_spike_count = zeros(nCells,1);
sta_info.kernel_available = false(nCells,1);
sta_info.no_spike = false(nCells,1);

% New onset diagnostics. These describe onset on the averaged, peak-aligned
% STA; they are NOT used to align the individual spikes.
sta_info.onset_available = false(nCells,1);
sta_info.onset_index = nan(nCells,1);
sta_info.onset_to_peak_frames = nan(nCells,1);
sta_info.onset_used_fallback = false(nCells,1);
sta_info.STA_observed_peak_index = nan(nCells,1);

%% ========================================================================
% 5. High-pass filter all traces
% =========================================================================
med_traces = movmedian(traces_nm, median_filter_size, 2, 'omitnan');
highpass_traces = traces_nm - med_traces;

%% ========================================================================
% 6. Process each cell
% =========================================================================
for iCell = 1:nCells

    % Skip preprocessing-invalid cells without changing indexing.
    if ~valid_preprocessing(iCell)
        continue
    end

    local_trace = double(highpass_traces(iCell,:));

    %% --------------------------------------------------------------------
    % 6A. Estimate baseline noise
    % ---------------------------------------------------------------------
    baseline_trace = local_trace(baseline_frames);
    noise_std = std(baseline_trace, 0, 'omitnan');

    if ~isfinite(noise_std) || noise_std <= 0
        sta_info.no_spike(iCell) = true;
        continue
    end

    %% --------------------------------------------------------------------
    % 6B. Initial high-SNR spike detection
    % ---------------------------------------------------------------------
    detection_trace = local_trace;
    detection_trace(~isfinite(detection_trace)) = 0;

    threshold = spike_SNR_th * noise_std;
    detection_trace(detection_trace < threshold) = 0;

    [~, local_spikes] = findpeaks( ...
        detection_trace, ...
        'MinPeakDistance', min_peak_distance);

    local_spikes = local_spikes(:)';
    sta_info.initial_spike_count(iCell) = numel(local_spikes);

    if isempty(local_spikes)
        sta_info.no_spike(iCell) = true;
        continue
    end

    %% --------------------------------------------------------------------
    % 6C. Remove BOTH members of any too-close spike pair
    % ---------------------------------------------------------------------
    if numel(local_spikes) > 1
        close_pairs = find(diff(local_spikes) <= min_spike_interval);
        remove_spike = false(size(local_spikes));
        remove_spike(close_pairs) = true;
        remove_spike(close_pairs + 1) = true;
        local_spikes(remove_spike) = [];
    end

    if isempty(local_spikes)
        sta_info.no_spike(iCell) = true;
        continue
    end

    %% ====================================================================
    % 7. Extract PEAK-ALIGNED spike traces
    % =====================================================================
    % Critical change from the old method:
    % do NOT locate d2V/dt max for every noisy spike and do NOT align each
    % spike to that value. The detected AP peak is the alignment anchor.
    local_peak_traces = [];

    for iSpike = 1:numel(local_spikes)
        spike_frame = local_spikes(iSpike);

        % Edge exclusion.
        if spike_frame <= edge_start || spike_frame >= nFrames - edge_end
            continue
        end

        extract_start = spike_frame - extract_pre;
        extract_end = spike_frame + extract_post;

        if extract_start < 1 || extract_end > nFrames
            continue
        end

        single_peak_trace = local_trace(extract_start:extract_end);

        if any(~isfinite(single_peak_trace))
            continue
        end

        local_peak_traces(end+1,:) = single_peak_trace; %#ok<AGROW>
    end

    %% ====================================================================
    % 8. Record valid spike count
    % =====================================================================
    nValidSpikes = size(local_peak_traces,1);
    sta_info.valid_spike_count(iCell) = nValidSpikes;

    if nValidSpikes == 0
        sta_info.no_spike(iCell) = true;
        continue
    end

    %% ====================================================================
    % 9. Average PEAK-ALIGNED spikes
    % =====================================================================
    STA_peak_aligned = mean(local_peak_traces, 1, 'omitnan');

    % Diagnostic: where is the maximum of the averaged STA?
    [~, observed_peak_idx] = max(STA_peak_aligned);
    sta_info.STA_observed_peak_index(iCell) = observed_peak_idx;

    %% ====================================================================
    % 10. Detect onset ONCE on the averaged STA
    % =====================================================================
    [STA_onset_idx, onset_info] = detect_AP_onset( ...
        STA_peak_aligned, ...
        expected_peak_idx, ...
        onset_cfg);

    if ~onset_info.fail && isfinite(STA_onset_idx)
        sta_info.onset_available(iCell) = true;
        sta_info.onset_index(iCell) = STA_onset_idx;
        sta_info.onset_to_peak_frames(iCell) = ...
            expected_peak_idx - STA_onset_idx;
        sta_info.onset_used_fallback(iCell) = onset_info.used_fallback;
    end

    %% ====================================================================
    % 11. Build peak-centered matched-filter kernel
    % =====================================================================
    source_start = expected_peak_idx - kernel_half_width;
    source_end = expected_peak_idx + kernel_half_width;

    new_kernel = STA_peak_aligned(source_start:source_end);

    if any(~isfinite(new_kernel))
        continue
    end

    kernel(iCell,:) = cast(new_kernel, 'like', traces_nm);
    sta_info.kernel_available(iCell) = true;
end

end
