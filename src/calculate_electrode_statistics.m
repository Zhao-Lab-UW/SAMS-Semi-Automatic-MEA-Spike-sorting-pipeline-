function electrode_stats = calculate_electrode_statistics(electrode_results, electrode_unit_stats, current_spikes, current_times, maxTime)
% CALCULATE_ELECTRODE_STATISTICS - Calculate comprehensive electrode-level statistics
%
% INPUTS:
%   electrode_results - Structure with electrode results from process_electrode
%   electrode_unit_stats - Structure with unit statistics
%   current_spikes - Spike waveforms matrix
%   current_times - Spike timing matrix
%   maxTime - Total recording duration
%
% OUTPUTS:
%   electrode_stats - Structure with comprehensive electrode statistics

    electrode_stats = struct();
    
    % Basic counts
    total_spikes = size(current_spikes, 2);
    new_idx_list = electrode_results.new_idx_list;
    num_units = length(new_idx_list);
    
    if num_units == 0
        % No units detected - return default values
        electrode_stats = create_empty_electrode_stats();
        electrode_stats.total_spikes = total_spikes;
        return;
    end
    
    % % Calculate sorted and unsorted spikes
    % sorted_spikes = [new_idx_list{:}];
    % unsorted_spikes = setdiff(1:total_spikes, sorted_spikes);
    % num_excluded_spikes = length(unsorted_spikes);
    % num_sorted_spikes = length(sorted_spikes);
    % Calculate sorted and unsorted spikes
    if isempty(new_idx_list)
        sorted_spikes = [];
        num_sorted_spikes = 0;
    else
        % Safe concatenation handling mixed orientations
        sorted_spikes = [];
        for unit_idx = 1:length(new_idx_list)
            if ~isempty(new_idx_list{unit_idx})
                unit_spikes = new_idx_list{unit_idx};
                % Ensure it's a row vector and concatenate
                sorted_spikes = [sorted_spikes, unit_spikes(:)'];
            end
        end
        num_sorted_spikes = length(sorted_spikes);
    end

    % Calculate unsorted spikes
    if total_spikes > 0 && ~isempty(sorted_spikes)
        unsorted_spikes = setdiff(1:total_spikes, sorted_spikes);
        num_excluded_spikes = length(unsorted_spikes);
    else
        unsorted_spikes = 1:total_spikes;
        num_excluded_spikes = total_spikes;
    end
    % Basic electrode metrics
    electrode_stats.total_spikes = total_spikes;
    electrode_stats.num_excluded_spikes = num_excluded_spikes;
    electrode_stats.unsorted_firing_rate = total_spikes / maxTime;
    electrode_stats.num_units_detected = num_units;
    
    % Calculate unit-specific statistics
    unit_firing_rates = [];
    unit_spike_counts = [];
    unit_burst_counts = [];
    unit_burst_durations = [];
    unit_spikes_per_burst = [];
    unit_inter_burst_intervals = [];
    all_unit_spike_times = [];
    
    for unit_idx = 1:num_units
        unit_spikes = new_idx_list{unit_idx};
        unit_spike_times = current_times(1, unit_spikes);
        all_unit_spike_times = [all_unit_spike_times, unit_spike_times];
        
        % Unit firing rate calculation (using actual spike time span)
        if length(unit_spike_times) > 1
            time_span = maxTime;
            if time_span > 0
                firing_rate = length(unit_spike_times) / time_span;
            else
                firing_rate = length(unit_spike_times) / 0.001; % Default 1ms if all spikes at same time
            end
        else
            firing_rate = NaN; % Cannot calculate meaningful firing rate for single spike
        end
        
        unit_firing_rates(end+1) = firing_rate;
        unit_spike_counts(end+1) = length(unit_spikes);
        % Extract burst statistics if available
        if isfield(electrode_unit_stats, 'active_units') && unit_idx <= length(electrode_unit_stats.active_units)
            unit_stats = electrode_unit_stats.active_units(unit_idx).stats;
            
            % Extract burst-related statistics (assuming they are in the stats array)
            if length(unit_stats) >= 5   
                num_bursts = unit_stats(5);  % Number of bursts
                unit_burst_counts(end+1) = num_bursts;
                
                if length(unit_stats) >= 7 && num_bursts > 0  
                    avg_burst_duration = unit_stats(7);  % Burst duration mean
                    unit_burst_durations(end+1) = avg_burst_duration;
                else
                    unit_burst_durations(end+1) = 0;
                end
                
                if length(unit_stats) >= 12 && num_bursts > 0  
                    avg_spikes_per_burst = unit_stats(12);  % Spikes per burst mean
                    unit_spikes_per_burst(end+1) = avg_spikes_per_burst;
                else
                    unit_spikes_per_burst(end+1) = 0;
                end
                
                if length(unit_stats) >= 19 && num_bursts > 1  
                    avg_inter_burst_interval = unit_stats(19);  % Inter-burst interval mean
                    unit_inter_burst_intervals(end+1) = avg_inter_burst_interval;
                else
                    unit_inter_burst_intervals(end+1) = 0;
                end
            else
                unit_burst_counts(end+1) = 0;
                unit_burst_durations(end+1) = 0;
                unit_spikes_per_burst(end+1) = 0;
                unit_inter_burst_intervals(end+1) = 0;
            end
        else
            unit_burst_counts(end+1) = 0;
            unit_burst_durations(end+1) = 0;
            unit_spikes_per_burst(end+1) = 0;
            unit_inter_burst_intervals(end+1) = 0;
        end
        
    end
    
    % Remove NaN values for statistics calculation
    valid_firing_rates = unit_firing_rates(~isnan(unit_firing_rates));
    
    % Average spike frequency per electrode
    electrode_stats.avg_spike_frequency = sum(unit_spike_counts) / maxTime;
    
    % Combined firing rate for all units
    if ~isempty(all_unit_spike_times)
        electrode_time_span = max(all_unit_spike_times) - min(all_unit_spike_times);
        if electrode_time_span > 0
            electrode_stats.electrode_firing_rate = length(all_unit_spike_times) / electrode_time_span;
        else
            electrode_stats.electrode_firing_rate = length(all_unit_spike_times) / 0.001;
        end
        
        % Electrode activity metrics
        electrode_stats.electrode_activity_span = electrode_time_span;
        electrode_stats.electrode_activity_percentage = (electrode_time_span / maxTime) * 100;
    else
        electrode_stats.electrode_firing_rate = 0;
        electrode_stats.electrode_activity_span = 0;
        electrode_stats.electrode_activity_percentage = 0;
    end
    
    % Unit firing rate distributions
    if ~isempty(valid_firing_rates)
        electrode_stats.unit_firing_rates_mean = mean(valid_firing_rates);
        electrode_stats.unit_firing_rates_std = std(valid_firing_rates);
        electrode_stats.unit_firing_rates_median = median(valid_firing_rates);
        electrode_stats.unit_firing_rates_min = min(valid_firing_rates);
        electrode_stats.unit_firing_rates_max = max(valid_firing_rates);
    else
        electrode_stats.unit_firing_rates_mean = NaN;
        electrode_stats.unit_firing_rates_std = NaN;
        electrode_stats.unit_firing_rates_median = NaN;
        electrode_stats.unit_firing_rates_min = NaN;
        electrode_stats.unit_firing_rates_max = NaN;
    end
    
    % Unit spike count distributions
    electrode_stats.unit_spike_counts_mean = mean(unit_spike_counts);
    electrode_stats.unit_spike_counts_std = std(unit_spike_counts);
    electrode_stats.unit_spike_counts_median = median(unit_spike_counts);
    electrode_stats.unit_spike_counts_min = min(unit_spike_counts);
    electrode_stats.unit_spike_counts_max = max(unit_spike_counts);
    
    % Burst characteristics
    bursting_units = sum(unit_burst_counts > 0);
    electrode_stats.num_bursting_units = bursting_units;
    electrode_stats.percentage_bursting_units = (bursting_units / num_units) * 100;
    
    if bursting_units > 0
        electrode_stats.avg_bursts_per_unit = mean(unit_burst_counts(unit_burst_counts > 0));
        electrode_stats.avg_burst_duration_electrode = mean(unit_burst_durations(unit_burst_durations > 0));
        electrode_stats.avg_spikes_per_burst_electrode = mean(unit_spikes_per_burst(unit_spikes_per_burst > 0));
        
        valid_inter_burst = unit_inter_burst_intervals(unit_inter_burst_intervals > 0);
        if ~isempty(valid_inter_burst)
            electrode_stats.avg_inter_burst_interval_electrode = mean(valid_inter_burst);
        else
            electrode_stats.avg_inter_burst_interval_electrode = 0;
        end
        
        % Electrode burst frequency (total bursts per time)
        total_bursts = sum(unit_burst_counts);
        electrode_stats.electrode_burst_frequency = total_bursts / maxTime;
    else
        electrode_stats.avg_bursts_per_unit = 0;
        electrode_stats.avg_burst_duration_electrode = 0;
        electrode_stats.avg_spikes_per_burst_electrode = 0;
        electrode_stats.avg_inter_burst_interval_electrode = 0;
        electrode_stats.electrode_burst_frequency = 0;
    end
    
    % Recording duration
    electrode_stats.total_recording_duration = maxTime;
end

function empty_stats = create_empty_electrode_stats()
    % Create structure with default/empty values for electrodes with no units
    empty_stats = struct();
    empty_stats.total_spikes = 0;
    empty_stats.num_excluded_spikes = 0;
    empty_stats.unsorted_firing_rate = 0;
    empty_stats.num_units_detected = 0;
    empty_stats.avg_spike_frequency = 0;
    empty_stats.electrode_firing_rate = 0;
    empty_stats.unit_firing_rates_mean = NaN;
    empty_stats.unit_firing_rates_std = NaN;
    empty_stats.unit_firing_rates_median = NaN;
    empty_stats.unit_firing_rates_min = NaN;
    empty_stats.unit_firing_rates_max = NaN;
    empty_stats.unit_spike_counts_mean = NaN;
    empty_stats.unit_spike_counts_std = NaN;
    empty_stats.unit_spike_counts_median = NaN;
    empty_stats.unit_spike_counts_min = NaN;
    empty_stats.unit_spike_counts_max = NaN;
    empty_stats.num_bursting_units = 0;
    empty_stats.percentage_bursting_units = 0;
    empty_stats.avg_bursts_per_unit = 0;
    empty_stats.avg_burst_duration_electrode = 0;
    empty_stats.avg_spikes_per_burst_electrode = 0;
    empty_stats.avg_inter_burst_interval_electrode = 0;
    empty_stats.electrode_burst_frequency = 0;
    empty_stats.total_recording_duration = 0;
    empty_stats.electrode_activity_span = 0;
    empty_stats.electrode_activity_percentage = 0;
end