function [T, T_electrode, T_parameters] = initialize_tables(params)
    % Create table for unit measurements - UPDATED ORDER
    Measurements = {
        % GENERAL SPIKE STATISTICS (1-4)
        'Number of Spikes';                                    % 1: total_spike_count
        'Mean Firing Rate (Hz)';                               % 2: mean_firing_rate
        'Mean ISI Overall (sec)';                              % 3: mean_ISI_overall (renamed from "ISI in general")
        'Overall ISI Coefficient of Variation';                % 4: CV_ISI_overall (NEW POSITION)
        
        % BURST OCCURRENCE (5-6)
        'Number of Bursts';                                    % 5: num_bursts
        'Burst Frequency (Hz)';                                % 6: burst_frequency
        
        % BURST DURATION (7-11)
        'Burst Duration - Mean (sec)';                         % 7: burst_duration_mean (renamed from Avg)
        'Burst Duration - Std (sec)';                          % 8: burst_duration_std
        'Burst Duration - Median (sec)';                       % 9: burst_duration_median
        'Percent Spikes in Bursts (%)';                        % 10: percent_spikes_in_bursts (renamed from "Burst Percentage")
        'Normalized Burst Duration IQR';                       % 11: burst_duration_IQR_normalized
        
        % SPIKES PER BURST (12-14)
        'Spikes per Burst - Mean';                             % 12: spikes_per_burst_mean (renamed from "Number of Spikes per Burst -Avg")
        'Spikes per Burst - Std';                              % 13: spikes_per_burst_std
        'Spikes per Burst - Median';                           % 14: spikes_per_burst_median
        
        % WITHIN-BURST ISI (15-18)
        'Within-Burst Mean ISI - Mean (sec)';                  % 15: within_burst_ISI_mean (renamed for clarity)
        'Within-Burst Mean ISI - Std (sec)';                   % 16: within_burst_ISI_std
        'Within-Burst Mean ISI - Median (sec)';                % 17: within_burst_ISI_median
        'Within-Burst ISI Coefficient of Variation';           % 18: within_burst_ISI_CV (renamed from "ISI Coefficient of Variation")
        
        % INTER-BURST INTERVALS (19-22)
        'Inter-Burst Interval - Mean (sec)';                   % 19: inter_burst_interval_mean (renamed from Avg)
        'Inter-Burst Interval - Std (sec)';                    % 20: inter_burst_interval_std
        'Inter-Burst Interval - Median (sec)';                 % 21: inter_burst_interval_median
        'Inter-Burst Interval Coefficient of Variation';       % 22: inter_burst_interval_CV
        
        % AMPLITUDE STATISTICS (23-26)
        'Peak Amplitude - Mean (mV)';                          % 23: amplitude_peak_mean (renamed from "Amplitude -Avg")
        'Peak Amplitude - Std (mV)';                           % 24: amplitude_peak_std
        'Peak-to-Trough Amplitude - Mean (mV)';                % 25: amplitude_peak_to_trough_mean (renamed for clarity)
        'Peak-to-Trough Amplitude - Std (mV)';                 % 26: amplitude_peak_to_trough_std
        
        % QUALITY METRICS (27-28)
        'Minimum ISI (ms)';                                    % 27: min_ISI_ms (renamed from "min ISI")
        'Number of Refractory Violations'                      % 28: num_refractory_violations (renamed from "number of removed short ISIs")
    };
    T = table(Measurements);

    % Create expanded table for electrode measurements
    Measurements_electrode = {
        % BASIC ELECTRODE STATISTICS
        'Total Number of Spikes';                              % Total spikes detected on electrode
        'Number of Excluded Spikes';                           % Spikes removed during processing
        'Unsorted Firing Rate (Hz)';                           % Raw firing rate before sorting
        'Number of Units Detected';                            % Total units identified
        
        % UNIT-LEVEL STATISTICS (aggregated across electrode)
        'Average Spike Frequency per Electrode (Hz)';          % Mean across all units
        'Electrode Firing Rate (all units combined, Hz)';      % Combined firing rate
        'Unit Firing Rates - Mean (Hz)';                       % Mean of unit firing rates
        'Unit Firing Rates - Std (Hz)';                        % Std of unit firing rates
        'Unit Firing Rates - Median (Hz)';                     % Median of unit firing rates
        'Unit Firing Rates - Min (Hz)';                        % Minimum unit firing rate
        'Unit Firing Rates - Max (Hz)';                        % Maximum unit firing rate
        'Unit Spike Counts - Mean';                            % Mean spike count per unit
        'Unit Spike Counts - Std';                             % Std of spike counts
        'Unit Spike Counts - Median';                          % Median spike count
        'Unit Spike Counts - Min';                             % Minimum spike count
        'Unit Spike Counts - Max';                             % Maximum spike count
        
        % BURST STATISTICS (electrode-level)
        'Number of Bursting Units';                            % Units showing burst activity
        'Percentage of Bursting Units (%)';                    % Percent of units that burst
        'Average Bursts per Unit';                             % Mean burst count across units
        'Average Burst Duration per Electrode (sec)';          % Mean burst duration
        'Average Spikes per Burst per Electrode';              % Mean spikes per burst
        'Inter-Burst Interval - Electrode Average (sec)';      % Mean IBI across electrode
        'Electrode Burst Frequency (Hz)';                      % Overall burst frequency
        
        % TEMPORAL STATISTICS
        'Total Recording Duration for Electrode (sec)';        % Total recording time
        'Electrode Activity Span (first to last spike, sec)';  % Time from first to last spike
        'Electrode Activity Percentage (%)'                    % Percent of time with activity
    };
    T_electrode = table(Measurements_electrode);
    
    % Create parameters table
    parameter_list = {'threshold_to_merge', 'refractoryT', ...
        'min_spikes_E', 'max_ISI_E', 'min_spikes_N', 'max_ISI_N',...
        'network_participation_threshold',...
        'cutoff_frequency', 'flag_threshold', 'std_cutoff',...
        'overlap_threshold'};
    T_parameters = table(parameter_list');
    
    % Add parameter values
    parameter_values = [params.threshold_to_merge; params.refractoryT; ...
        params.min_spikes_E; params.max_ISI_E; params.min_spikes_N; params.max_ISI_N; ...
        params.network_participation_threshold; ...
        params.cutoff_frequency; params.flag_threshold; params.std_cutoff;...
        params.overlap_threshold];
    Tleft_parameters = table(parameter_values);
    Tleft_parameters.Properties.VariableNames(1) = {'parameter Value'};
    T_parameters = [T_parameters, Tleft_parameters];
end