function [new_idx_list, unit_stats] = process_final_clusters(initial_idx_list, spikes, Times, params, maxTime,  minSpikeAmp, maxSpikeAmp)
% PROCESS_FINAL_CLUSTERS2 - Process final clusters after merging and perform burst analysis
%
% INPUTS:
%   initial_idx_list - Cell array of indices for each cluster
%   spikes - Spike waveforms matrix (timepoints × spikes)
%   Times - Spike timing vector
%   params - Parameter structure with fields:
%       .std_cutoff - Threshold factor for outlier removal
%       .refractoryT - Refractory period threshold
%       .cutoff_frequency - Minimum firing rate for active units
%       .maximum_inter_spike_interval - Max ISI for burst detection
%       .minimum_spike_per_burst - Min spikes per burst
%   maxTime - Maximum recording time
%   idxs_outlier - Indices of outlier WFs in initial_idx_list
%   minSpikeAmp - Minimum spike amplitude threshold
%   maxSpikeAmp - Maximum spike amplitude threshold
%
% OUTPUTS:
%   new_idx_list - Cell array of indices for active units
%   unit_stats - Structure with statistics for each unit

    
    % Initialize outputs
    unit_count = 0;
    new_idx_list = {};
    unit_stats = struct();
    unit_stats.active_units = [];
    unit_stats.inactive_units = [];
    unit_stats.active_units_fr = [];
    unit_stats.inactive_units_fr = [];

    % Process each cluster
    for k_list = 1:length(initial_idx_list)
        if isempty(initial_idx_list{k_list})
            continue;  % Skip empty clusters
        end
        
        % Get indices for this unit
        index_of_idx = sort(initial_idx_list{k_list});
        
        % Remove outliers based on RMSE (iterative cleaning)
        A = spikes(:, index_of_idx);
        max_iterations = size(A, 1);  % Maximum iterations based on number of time points
        iteration = 0;
        while iteration < max_iterations
            iteration = iteration + 1;
            A_RMSE = calculate_WF_RMSE(A);
            [~, TFoutlier] = rmoutliers(A_RMSE, 'mean', 'ThresholdFactor', params.std_cutoff);
            outlier_locs = find(TFoutlier == 1);
            if isempty(outlier_locs)
                break;  % No more outliers to remove
            end
            A(:, outlier_locs) = [];
            index_of_idx(outlier_locs) = [];
            if isempty(index_of_idx)
                break;  % No spikes left
            end
        end
        
        % Remove spikes violating refractory period
        if length(index_of_idx) > 1
            t_temp = Times(1, index_of_idx);
            ISI_s = diff(t_temp);
            
            % Find all violations
            violation_idx = find(ISI_s < params.refractoryT);
            num_refractory_violations = length(violation_idx);
            
            % Remove the second spike of each violating pair
            if ~isempty(violation_idx)
                index_of_idx(violation_idx + 1) = [];
            end
        else
            num_refractory_violations = 0;
        end
        
        % Remove invalid time points (<=0)
        if ~isempty(index_of_idx)
            t_temp = Times(1, index_of_idx);
            valid_idx = t_temp > 0;
            index_of_idx = index_of_idx(valid_idx);
        end
        
        % Update the indices in the initial_idx_list
        initial_idx_list{k_list} = index_of_idx;
        
        % Count number of spikes for this unit
        number_of_spikes = numel(index_of_idx);
        
        % Skip if no spikes remain
        if number_of_spikes == 0
            continue;
        end
        
        % Calculate firing rate and check if unit is active
        firing_rate = number_of_spikes / maxTime;
        
        if firing_rate > params.cutoff_frequency
            % This is an active unit
            unit_count = unit_count + 1;
            
            % Initialize flags
            unit_stats.active_units(unit_count).flagLargeWFs = 0;
            unit_stats.active_units(unit_count).flagOutlierWFs = 0;
            
            % Add to active units list
            new_idx_list{unit_count} = index_of_idx;
            
            % Get clean spike times for this unit
            t = Times(1, index_of_idx);
            
            % ===== GENERAL SPIKE STATISTICS =====
            % Calculate basic firing statistics
            total_spike_count = number_of_spikes;
            mean_firing_rate = firing_rate;
            mean_ISI_overall = maxTime / number_of_spikes;  % Average ISI across entire recording
            
            % Calculate overall ISI variability (CV of all ISIs)
            CV_ISI_overall = NaN;  % Default value
            if length(t) > 1
                all_ISIs = diff(t);
                mean_all_ISI = mean(all_ISIs);
                std_all_ISI = std(all_ISIs);
                if mean_all_ISI > 0
                    CV_ISI_overall = std_all_ISI / mean_all_ISI;
                end
            end
            
            % ===== AMPLITUDE STATISTICS =====
            spikes_amplitude = spikes(:, index_of_idx);
            
            
            % Convert to mV for statistics
            spikes_amplitude = spikes_amplitude * 1000;
            spikes_amplitude_mean = mean(spikes_amplitude, 2);
            
            % Peak amplitude statistics
            [amplitude_peak_mean, amplitude_peak_loc] = max(abs(spikes_amplitude_mean));
            amplitude_peak_std = std(spikes_amplitude(amplitude_peak_loc, :));
            
            % Peak-to-trough amplitude
            [waveform_peak, waveform_peak_loc] = max(spikes_amplitude_mean);
            [waveform_trough, waveform_trough_loc] = min(spikes_amplitude_mean);
            amplitude_peak_to_trough_mean = waveform_peak - waveform_trough;
            
            % Calculate std 
            if length(index_of_idx) > 1
                amplitude_peak_to_trough_std = std(spikes_amplitude(waveform_peak_loc, :) - ...
                    spikes_amplitude(waveform_trough_loc, :));
            else
                amplitude_peak_to_trough_std = 0;
            end
            
            % Calculate shortest ISI
            min_ISI_ms = NaN;
            if length(t) > 1
                ISI_ms = diff(t) * 1000;
                min_ISI_ms = min(ISI_ms);
            end
            
            % ===== BURST ANALYSIS =====
            % Get burst information using enhanced get_burst (with 14 rows)
            burst_info = get_burst(t, params.maximum_inter_spike_interval, params.minimum_spike_per_burst);
            
            % Initialize statistics vector with 28 elements
            stats_vector = zeros(28, 1);
            
            % GENERAL SPIKE STATISTICS (indices 1-4)
            stats_vector(1) = total_spike_count;          % Total number of spikes
            stats_vector(2) = mean_firing_rate;           % Mean firing rate (Hz)
            stats_vector(3) = mean_ISI_overall;           % Mean ISI across entire recording (s)
            stats_vector(4) = CV_ISI_overall;             % CV of all ISIs (NEW POSITION)
            
            if ~isempty(burst_info) && size(burst_info, 2) > 0
                % ===== BURST-RELATED STATISTICS =====
                num_bursts = size(burst_info, 2);
                
                % BURST OCCURRENCE (indices 5-6)
                stats_vector(5) = num_bursts;                           % Total number of bursts
                stats_vector(6) = num_bursts / maxTime;                 % Burst frequency (bursts/s)
                
                % BURST DURATION (indices 7-11)
                stats_vector(7) = mean(burst_info(3, :));               % Mean burst duration (s)
                stats_vector(8) = std(burst_info(3, :));                % STD burst duration (s)
                stats_vector(9) = median(burst_info(3, :));             % Median burst duration (s)
                stats_vector(10) = sum(burst_info(2, :)) / number_of_spikes * 100;  % Percent spikes in bursts
                
                % Normalized duration IQR
                if stats_vector(9) > 0
                    stats_vector(11) = iqr(burst_info(3, :)) / stats_vector(9);  % Normalized IQR of duration
                else
                    stats_vector(11) = NaN;
                end
                
                % SPIKES PER BURST (indices 12-14)
                stats_vector(12) = mean(burst_info(2, :));              % Mean spikes per burst
                stats_vector(13) = std(burst_info(2, :));               % STD spikes per burst
                stats_vector(14) = median(burst_info(2, :));            % Median spikes per burst
                
                % WITHIN-BURST ISI (indices 15-18)
                stats_vector(15) = mean(burst_info(4, :));              % Mean of mean ISIs within bursts (s)
                stats_vector(16) = std(burst_info(4, :));               % STD of mean ISIs within bursts (s)
                stats_vector(17) = median(burst_info(4, :));            % Median of mean ISIs within bursts (s)
                stats_vector(18) = mean(burst_info(7, :));              % Mean CV of ISI within bursts
                
                % INTER-BURST INTERVALS (indices 19-22)
                if num_bursts > 1
                    % Row 5 contains inter-burst intervals, last one is NaN
                    valid_IBIs = burst_info(5, ~isnan(burst_info(5, :)));
                    
                    if ~isempty(valid_IBIs)
                        stats_vector(19) = mean(valid_IBIs);            % Mean inter-burst interval (s)
                        stats_vector(20) = std(valid_IBIs);             % STD inter-burst interval (s)
                        stats_vector(21) = median(valid_IBIs);          % Median inter-burst interval (s)
                        
                        if stats_vector(19) > 0
                            stats_vector(22) = stats_vector(20) / stats_vector(19);  % CV inter-burst interval
                        else
                            stats_vector(22) = NaN;
                        end
                    else
                        stats_vector(19:22) = NaN;
                    end
                else
                    % Only one burst - no inter-burst intervals
                    stats_vector(19:22) = NaN;
                end
                
            else
                % No bursts detected
                stats_vector(5:18) = 0;     % Set burst stats to 0
                stats_vector(19:22) = NaN;   % Inter-burst intervals to NaN
            end
            
            % AMPLITUDE STATISTICS (indices 23-26)
            stats_vector(23) = amplitude_peak_mean;                     % Mean peak amplitude (mV)
            stats_vector(24) = amplitude_peak_std;                      % STD peak amplitude (mV)
            stats_vector(25) = amplitude_peak_to_trough_mean;           % Mean peak-to-trough amplitude (mV)
            stats_vector(26) = amplitude_peak_to_trough_std;            % STD peak-to-trough amplitude (mV)
            
            % QUALITY METRICS (indices 27-28)
            stats_vector(27) = min_ISI_ms;                              % Minimum ISI (ms)
            stats_vector(28) = num_refractory_violations;               % Number of refractory violations
            
            % Store statistics
            unit_stats.active_units(unit_count).stats = stats_vector;
            unit_stats.active_units(unit_count).index = k_list;
            unit_stats.active_units_fr(unit_count) = mean_firing_rate;
            
        else
            % This is an inactive unit (firing rate below cutoff)
            inactive_idx = length(unit_stats.inactive_units) + 1;
            unit_stats.inactive_units(inactive_idx).index = k_list;
            unit_stats.inactive_units(inactive_idx).indices = index_of_idx;
            unit_stats.inactive_units_fr(inactive_idx) = firing_rate;
        end
    end
    
    % Final validation
    if unit_count == 0
        warning('No active units found after processing');
    else
        fprintf('Processed %d active units\n', unit_count);
    end
end