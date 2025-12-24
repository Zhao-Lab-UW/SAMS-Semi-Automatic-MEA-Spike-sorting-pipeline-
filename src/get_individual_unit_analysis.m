function [T_unit, T_electrode] = get_individual_unit_analysis(sorting_results, allData, maxTime, ...
    cutoffFrequency, maxISI, minSpikesPerBurst, refractoryT, outputFolder)
% GET_INDIVIDUAL_UNIT_ANALYSIS - Extract and analyze individual unit statistics
% This should exactly matches calculations from process_final_clusters.m and calculate_electrode_statistics.m
%
% INPUTS:
%   sorting_results - Cell array containing sorted unit indices (manually reclassified)
%   allData - Cell array containing spike data from all electrodes
%   maxTime - Maximum recording duration in seconds
%   cutoffFrequency - Cutoff frequency for filtering (Hz) - units with firing rates 
%                     below this value will be excluded from active unit analysis
%   maxISI - Maximum inter-spike interval for burst detection (ms)
%   minSpikesPerBurst - Minimum spikes per burst
%   outputFolder - Folder to save results
%
% OUTPUTS:
%   T_unit - Table with individual unit statistics (27 measurements per unit)
%   T_electrode - Table with electrode-level statistics (26 measurements per electrode)

    % Define the measurements to track for individual units (28 measurements matching initialize_tables.m)
    unitMeasurements = {
        % GENERAL SPIKE STATISTICS (1-4)
        'Number of Spikes';                                    % 1
        'Mean Firing Rate (Hz)';                               % 2
        'Mean ISI Overall (sec)';                              % 3
        'Overall ISI Coefficient of Variation';                % 4
        
        % BURST OCCURRENCE (5-6)
        'Number of Bursts';                                    % 5
        'Burst Frequency (Hz)';                                % 6
        
        % BURST DURATION (7-11)
        'Burst Duration - Mean (sec)';                         % 7
        'Burst Duration - Std (sec)';                          % 8
        'Burst Duration - Median (sec)';                       % 9
        'Percent Spikes in Bursts (%)';                        % 10
        'Normalized Burst Duration IQR';                       % 11
        
        % SPIKES PER BURST (12-14)
        'Spikes per Burst - Mean';                             % 12
        'Spikes per Burst - Std';                              % 13
        'Spikes per Burst - Median';                           % 14
        
        % WITHIN-BURST ISI (15-18)
        'Within-Burst Mean ISI - Mean (sec)';                  % 15
        'Within-Burst Mean ISI - Std (sec)';                   % 16
        'Within-Burst Mean ISI - Median (sec)';                % 17
        'Within-Burst ISI Coefficient of Variation';           % 18
        
        % INTER-BURST INTERVALS (19-22)
        'Inter-Burst Interval - Mean (sec)';                   % 19
        'Inter-Burst Interval - Std (sec)';                    % 20
        'Inter-Burst Interval - Median (sec)';                 % 21
        'Inter-Burst Interval Coefficient of Variation';       % 22
        
        % AMPLITUDE STATISTICS (23-26)
        'Peak Amplitude - Mean (mV)';                          % 23
        'Peak Amplitude - Std (mV)';                           % 24
        'Peak-to-Trough Amplitude - Mean (mV)';                % 25
        'Peak-to-Trough Amplitude - Std (mV)';                 % 26
        
        % QUALITY METRICS (27-28)
        'Minimum ISI (ms)';                                    % 27
        'Number of Refractory Violations'                      % 28
    };
    
    % Define electrode-level measurements (26 measurements matching initialize_tables.m)
    electrodeMeasurements = {
        'Total Number of Spikes';
        'Number of Excluded Spikes';
        'Unsorted Firing Rate (Hz)';  
        'Number of Units Detected';
        'Average Spike Frequency per Electrode (Hz)';
        'Electrode Firing Rate (all units combined, Hz)';
        'Unit Firing Rates - Mean (Hz)';
        'Unit Firing Rates - Std (Hz)';
        'Unit Firing Rates - Median (Hz)';
        'Unit Firing Rates - Min (Hz)';
        'Unit Firing Rates - Max (Hz)';
        'Unit Spike Counts - Mean';
        'Unit Spike Counts - Std';
        'Unit Spike Counts - Median';
        'Unit Spike Counts - Min';
        'Unit Spike Counts - Max';
        'Number of Bursting Units';
        'Percentage of Bursting Units (%)';
        'Average Bursts per Unit';
        'Average Burst Duration per Electrode (sec)';
        'Average Spikes per Burst per Electrode';
        'Inter-Burst Interval - Electrode Average (sec)';
        'Electrode Burst Frequency (Hz)';
        'Total Recording Duration for Electrode (sec)';
        'Electrode Activity Span (first to last spike, sec)';
        'Electrode Activity Percentage (%)'
    };
    
    % Initialize output matrices
    unitData = zeros(length(unitMeasurements), 0);
    electrodeData = zeros(length(electrodeMeasurements), 0);
    
    % Arrays to store column names
    unitNames = {};
    electrodeNames = {};
    
    % Get dimensions of sorting_results array
    [nwr, nwc, nec, ner, ~] = size(sorting_results);
    
    % Convert maxISI from ms to seconds for internal calculations
    maxISI_seconds = maxISI / 1000;


    % Loop through all electrodes
    for i = 1:nwr
        for j = 1:nwc
            for m = 1:nec
                for n = 1:ner
                    % Skip if no sorting results for this electrode
                    if isempty(sorting_results{i,j,m,n,1})
                        continue;
                    end

                    % Get electrode name
                    electrodeName = [char(i+'A'-1), num2str(j), '_', num2str(m), num2str(n)];

                    % Get spike data for this electrode
                    Spikes = allData{i,j,m,n}(:);
                    if isempty(Spikes)
                        continue;
                    end

                    % Get spike times and waveforms
                    [Times, spikes] = Spikes.GetTimeVoltageVector;

                    % Get sorted units and unsorted waveforms
                    sorted_units = sorting_results{i,j,m,n,1};
                    unsorted_wfs = sorting_results{i,j,m,n,2};

                    % Count total units detected (including inactive ones)
                    total_units_detected = length(sorted_units);

                    % Initialize electrode-level statistics collectors
                    unit_firing_rates = [];  % For electrode stats calculation
                    unit_spike_counts = [];
                    unit_burst_counts = [];
                    unit_burst_durations = [];
                    unit_spikes_per_burst = [];
                    unit_inter_burst_intervals = [];
                    all_electrode_spike_times = [];

                    % Count active units for this electrode
                    active_unit_count = 0;

                    % Process each unit for this electrode
                    for unit_idx = 1:length(sorted_units)
                        if isempty(sorted_units{unit_idx})
                            continue;
                        end

                        % Get spike indices for this unit
                        unit_spike_indices = sorted_units{unit_idx};

                        % Get spike times for this unit
                        unit_times_raw = Times(1,unit_spike_indices);
                        
                        % Remove spikes violating refractory period (matching process_final_clusters2.m)
                        removed_short_isis = 0;
                        if length(unit_spike_indices) > 1
                            % Sort indices and times together
                            [unit_times_sorted, sort_idx] = sort(unit_times_raw);
                            unit_spike_indices_sorted = unit_spike_indices(sort_idx);
                            
                            % Calculate ISIs
                            isis_seconds = diff(unit_times_sorted);
      
                            % Find violations
                            violation_idx = find(isis_seconds < refractoryT);
                            removed_short_isis = length(violation_idx);
                            
                            % Remove the second spike of each violating pair
                            if ~isempty(violation_idx)
                                unit_spike_indices_sorted(violation_idx + 1) = [];
                                unit_times_sorted(violation_idx + 1) = [];
                            end
                            
                            % Update cleaned indices and times
                            unit_spike_indices = unit_spike_indices_sorted;
                            unit_times = unit_times_sorted;
                        else
                            unit_times = unit_times_raw;
                        end

                        % Number of spikes (after removing violations)
                        num_spikes = length(unit_times);

                        % Calculate mean firing rate (matching process_final_clusters2.m)
                        mean_fr = num_spikes / maxTime;
                        
                        % Check if unit meets cutoff frequency
                        if mean_fr < cutoffFrequency
                            fprintf('Unit %s firing rate (%.2f Hz) below cutoff (%.2f Hz) - excluding from analysis\n', ...
                                [electrodeName, '_', num2str(unit_idx)], mean_fr, cutoffFrequency);
                            continue;
                        end

                        % This unit is active
                        active_unit_count = active_unit_count + 1;
                        unitName = [electrodeName, '_', num2str(active_unit_count)];

                        % Get waveforms for this unit
                        unit_waveforms = spikes(:, unit_spike_indices);

                        % Add to electrode-level spike time collection (for activity span calc)
                        all_electrode_spike_times = [all_electrode_spike_times, unit_times];

                        % Calculate ISI in general (matching process_final_clusters2.m line 43)
                        % This is mean ISI across entire recording
                        isi_in_general = maxTime / num_spikes;

                        % Calculate ISIs in ms (after refractory violations removed)
                        if length(unit_times) > 1
                            isis_seconds = diff(unit_times);
                            isis_ms = isis_seconds * 1000;
                            shortest_isi_ms = min(isis_ms);
                            
                            % Calculate overall ISI CV (matching process_final_clusters2.m)
                            mean_all_ISI = mean(isis_seconds);
                            std_all_ISI = std(isis_seconds);
                            if mean_all_ISI > 0
                                overall_isi_cv = std_all_ISI / mean_all_ISI;
                            else
                                overall_isi_cv = NaN;
                            end
                        else
                            shortest_isi_ms = NaN;
                            overall_isi_cv = NaN;
                        end

                        % Calculate waveform amplitude statistics (matching process_final_clusters.m)
                        if ~isempty(unit_waveforms)
                            % Convert to mV first (matching process_final_clusters.m line 127-128)
                            spikes_amplitude = unit_waveforms * 1000;
                            spikes_amplitude_mean = mean(spikes_amplitude, 2);

                            % Find maximum absolute value location (peak amplitude)
                            [amplitude_avg, amplitude_avg_loc] = max(abs(spikes_amplitude_mean));
                            amplitude_std = std(spikes_amplitude(amplitude_avg_loc, :));

                            % Calculate peak-to-trough amplitude
                            [spikes_amplitude_mean_peak, spikes_amplitude_mean_peak_loc] = max(spikes_amplitude_mean);
                            [spikes_amplitude_mean_trough, spikes_amplitude_mean_trough_loc] = min(spikes_amplitude_mean);
                            amplitude_avg_peak_to_trough = spikes_amplitude_mean_peak - spikes_amplitude_mean_trough;
                            if length(unit_spike_indices) > 1
                                amplitude_std_peak_to_trough = std(spikes_amplitude(spikes_amplitude_mean_peak_loc, :) - ...
                                    spikes_amplitude(spikes_amplitude_mean_trough_loc, :));
                            else
                                amplitude_std_peak_to_trough = 0;
                            end
                        else
                            amplitude_avg = NaN;
                            amplitude_std = NaN;
                            amplitude_avg_peak_to_trough = NaN;
                            amplitude_std_peak_to_trough = NaN;
                        end

                        % Detect bursts using the get_burst function
                        if length(unit_times) >= minSpikesPerBurst
                            burstInfo = get_burst(unit_times, maxISI_seconds, minSpikesPerBurst);
                        else
                            burstInfo = [];
                        end

                        % Calculate burst statistics (matching process_final_clusters2.m)
                        if ~isempty(burstInfo) && size(burstInfo, 2) > 0
                            num_bursts = size(burstInfo, 2);

                            % Burst frequency (matching process_final_clusters2.m)
                            burst_freq = num_bursts / maxTime;

                            % Burst duration stats
                            burst_durations = burstInfo(3, :);
                            burst_duration_avg = mean(burst_durations);
                            burst_duration_std = std(burst_durations);
                            burst_duration_median = median(burst_durations);
                            burst_duration_iqr = iqr(burst_durations);
                            if burst_duration_median > 0
                                norm_burst_duration_iqr = burst_duration_iqr / burst_duration_median;
                            else
                                norm_burst_duration_iqr = NaN;
                            end

                            % Spikes per burst stats
                            spikes_per_burst = burstInfo(2, :);
                            spikes_per_burst_avg = mean(spikes_per_burst);
                            spikes_per_burst_std = std(spikes_per_burst);
                            spikes_per_burst_median = median(spikes_per_burst);

                            % Calculate burst percentage
                            total_spikes_in_bursts = sum(spikes_per_burst);
                            burst_percentage = (total_spikes_in_bursts / num_spikes) * 100;

                            % ISI within burst stats (mean of mean ISIs within bursts)
                            isi_within_burst = burstInfo(4, :);
                            isi_within_burst_avg = mean(isi_within_burst);
                            isi_within_burst_std = std(isi_within_burst);
                            isi_within_burst_median = median(isi_within_burst);

                            % Within-burst ISI CV (from row 7 of burstInfo if available)
                            if size(burstInfo, 1) >= 7
                                within_burst_isi_cv = mean(burstInfo(7, :));  % Mean CV of ISI within bursts
                            else
                                % Calculate if not available
                                if isi_within_burst_avg > 0
                                    within_burst_isi_cv = isi_within_burst_std / isi_within_burst_avg;
                                else
                                    within_burst_isi_cv = NaN;
                                end
                            end

                            % Inter-burst interval stats
                            if num_bursts > 1
                                % Get valid IBIs from row 5 (excluding NaN)
                                valid_IBIs = burstInfo(5, ~isnan(burstInfo(5, :)));
                                
                                if ~isempty(valid_IBIs)
                                    ibi_avg = mean(valid_IBIs);
                                    ibi_std = std(valid_IBIs);
                                    ibi_median = median(valid_IBIs);
                                    if ibi_avg > 0
                                        ibi_cv = ibi_std / ibi_avg;
                                    else
                                        ibi_cv = NaN;
                                    end
                                else
                                    ibi_avg = NaN;
                                    ibi_std = NaN;
                                    ibi_median = NaN;
                                    ibi_cv = NaN;
                                end
                            else
                                ibi_avg = NaN;
                                ibi_std = NaN;
                                ibi_median = NaN;
                                ibi_cv = NaN;
                            end

                            % Store for electrode-level stats
                            unit_burst_counts(end+1) = num_bursts;
                            unit_burst_durations(end+1) = burst_duration_avg;
                            unit_spikes_per_burst(end+1) = spikes_per_burst_avg;
                            if ~isnan(ibi_avg) && ibi_avg > 0
                                unit_inter_burst_intervals(end+1) = ibi_avg;
                            end
                        else
                            % No bursts detected - use zeros for burst stats, NaN for IBI (matching process_final_clusters.m)
                            num_bursts = 0;
                            burst_freq = 0;
                            burst_duration_avg = 0;
                            burst_duration_std = 0;
                            burst_duration_median = 0;
                            norm_burst_duration_iqr = 0;
                            spikes_per_burst_avg = 0;
                            spikes_per_burst_std = 0;
                            spikes_per_burst_median = 0;
                            isi_within_burst_avg = 0;
                            isi_within_burst_std = 0;
                            isi_within_burst_median = 0;
                            within_burst_isi_cv = 0;  % Within-burst ISI CV is 0 when no bursts
                            ibi_avg = NaN;  % Inter-burst intervals are NaN when no bursts
                            ibi_std = NaN;
                            ibi_median = NaN;
                            ibi_cv = NaN;
                            burst_percentage = 0;

                            % Store zeros for electrode-level stats
                            unit_burst_counts(end+1) = 0;
                        end

                        % Store unit firing rate for electrode stats (using time span like calculate_electrode_statistics.m)
                        if length(unit_times) > 1
                            time_span = maxTime;
                            if time_span > 0
                                unit_fr_for_electrode = length(unit_times) / time_span;
                            else
                                unit_fr_for_electrode = length(unit_times) / 0.001;
                            end
                            unit_firing_rates(end+1) = unit_fr_for_electrode;
                        elseif length(unit_times) == 1
                            % Single spike - don't include in firing rate stats
                            % This matches calculate_electrode_statistics.m which skips NaN values
                        end

                        unit_spike_counts(end+1) = num_spikes;

                        % Compile unit statistics in exact order matching initialize_tables.m (28 measurements)
                        unitNames{end+1} = unitName;
                        unitData(:, end+1) = [
                            % GENERAL SPIKE STATISTICS (1-4)
                            num_spikes;                      % 1: Number of Spikes
                            mean_fr;                         % 2: Mean Firing Rate (Hz)
                            isi_in_general;                  % 3: Mean ISI Overall (sec)
                            overall_isi_cv;                  % 4: Overall ISI Coefficient of Variation
                            
                            % BURST OCCURRENCE (5-6)
                            num_bursts;                      % 5: Number of Bursts
                            burst_freq;                      % 6: Burst Frequency (Hz)
                            
                            % BURST DURATION (7-11)
                            burst_duration_avg;              % 7: Burst Duration - Mean (sec)
                            burst_duration_std;              % 8: Burst Duration - Std (sec)
                            burst_duration_median;           % 9: Burst Duration - Median (sec)
                            burst_percentage;                % 10: Percent Spikes in Bursts (%)
                            norm_burst_duration_iqr;         % 11: Normalized Burst Duration IQR
                            
                            % SPIKES PER BURST (12-14)
                            spikes_per_burst_avg;            % 12: Spikes per Burst - Mean
                            spikes_per_burst_std;            % 13: Spikes per Burst - Std
                            spikes_per_burst_median;         % 14: Spikes per Burst - Median
                            
                            % WITHIN-BURST ISI (15-18)
                            isi_within_burst_avg;            % 15: Within-Burst Mean ISI - Mean (sec)
                            isi_within_burst_std;            % 16: Within-Burst Mean ISI - Std (sec)
                            isi_within_burst_median;         % 17: Within-Burst Mean ISI - Median (sec)
                            within_burst_isi_cv;             % 18: Within-Burst ISI Coefficient of Variation
                            
                            % INTER-BURST INTERVALS (19-22)
                            ibi_avg;                         % 19: Inter-Burst Interval - Mean (sec)
                            ibi_std;                         % 20: Inter-Burst Interval - Std (sec)
                            ibi_median;                      % 21: Inter-Burst Interval - Median (sec)
                            ibi_cv;                          % 22: Inter-Burst Interval Coefficient of Variation
                            
                            % AMPLITUDE STATISTICS (23-26)
                            amplitude_avg;                   % 23: Peak Amplitude - Mean (mV)
                            amplitude_std;                   % 24: Peak Amplitude - Std (mV)
                            amplitude_avg_peak_to_trough;    % 25: Peak-to-Trough Amplitude - Mean (mV)
                            amplitude_std_peak_to_trough;    % 26: Peak-to-Trough Amplitude - Std (mV)
                            
                            % QUALITY METRICS (27-28)
                            shortest_isi_ms;                 % 27: Minimum ISI (ms)
                            removed_short_isis               % 28: Number of Refractory Violations
                        ];
                    end

                    % Create electrode statistics even if no active units (matching calculate_electrode_statistics.m)
                    if total_units_detected > 0  % At least one unit was detected (active or inactive)
                        % Calculate electrode-level statistics (matching calculate_electrode_statistics.m)
                        total_spikes = size(spikes, 2);
                        excluded_spikes = length(unsorted_wfs);

                        % Unsorted firing rate - uses total_spikes (matching calculate_electrode_statistics.m)
                        unsorted_fr = total_spikes / maxTime;  

                        % Average spike frequency (sum of unit spike counts / maxTime)
                        if ~isempty(unit_spike_counts)
                            avg_spike_frequency = sum(unit_spike_counts) / maxTime;
                        else
                            avg_spike_frequency = 0;
                        end

                        % Electrode firing rate based on actual activity span
                        if ~isempty(all_electrode_spike_times)
                            electrode_time_span = max(all_electrode_spike_times) - min(all_electrode_spike_times);
                            if electrode_time_span > 0
                                electrode_firing_rate = length(all_electrode_spike_times) / electrode_time_span;
                            else
                                electrode_firing_rate = length(all_electrode_spike_times) / 0.001;
                            end
                            electrode_activity_span = electrode_time_span;
                            electrode_activity_percentage = (electrode_time_span / maxTime) * 100;
                        else
                            electrode_firing_rate = 0;
                            electrode_activity_span = 0;
                            electrode_activity_percentage = 0;
                        end

                        % Unit firing rate statistics (matching calculate_electrode_statistics.m)
                        if ~isempty(unit_firing_rates)
                            unit_fr_mean = mean(unit_firing_rates);
                            unit_fr_std = std(unit_firing_rates);
                            unit_fr_median = median(unit_firing_rates);
                            unit_fr_min = min(unit_firing_rates);
                            unit_fr_max = max(unit_firing_rates);
                        else
                            unit_fr_mean = NaN;
                            unit_fr_std = NaN;
                            unit_fr_median = NaN;
                            unit_fr_min = NaN;
                            unit_fr_max = NaN;
                        end

                        % Unit spike count statistics
                        if ~isempty(unit_spike_counts)
                            unit_sc_mean = mean(unit_spike_counts);
                            unit_sc_std = std(unit_spike_counts);
                            unit_sc_median = median(unit_spike_counts);
                            unit_sc_min = min(unit_spike_counts);
                            unit_sc_max = max(unit_spike_counts);
                        else
                            unit_sc_mean = NaN;
                            unit_sc_std = NaN;
                            unit_sc_median = NaN;
                            unit_sc_min = NaN;
                            unit_sc_max = NaN;
                        end

                        % Burst statistics
                        if active_unit_count > 0
                            bursting_units = sum(unit_burst_counts > 0);
                            percentage_bursting_units = (bursting_units / active_unit_count) * 100;

                            if bursting_units > 0
                                avg_bursts_per_unit = mean(unit_burst_counts(unit_burst_counts > 0));
                                avg_burst_duration_electrode = mean(unit_burst_durations(unit_burst_durations > 0));
                                avg_spikes_per_burst_electrode = mean(unit_spikes_per_burst(unit_spikes_per_burst > 0));

                                if ~isempty(unit_inter_burst_intervals)
                                    avg_inter_burst_interval_electrode = mean(unit_inter_burst_intervals);
                                else
                                    avg_inter_burst_interval_electrode = 0;
                                end

                                total_bursts = sum(unit_burst_counts);
                                electrode_burst_frequency = total_bursts / maxTime;
                            else
                                avg_bursts_per_unit = 0;
                                avg_burst_duration_electrode = 0;
                                avg_spikes_per_burst_electrode = 0;
                                avg_inter_burst_interval_electrode = 0;
                                electrode_burst_frequency = 0;
                            end
                        else
                            bursting_units = 0;
                            percentage_bursting_units = 0;
                            avg_bursts_per_unit = 0;
                            avg_burst_duration_electrode = 0;
                            avg_spikes_per_burst_electrode = 0;
                            avg_inter_burst_interval_electrode = 0;
                            electrode_burst_frequency = 0;
                        end

                        % Add electrode statistics to data matrix (26 measurements)
                        electrodeNames{end+1} = electrodeName;
                        electrodeData(:, end+1) = [
                            total_spikes;                         % 1
                            excluded_spikes;                      % 2
                            unsorted_fr;                          % 3
                            total_units_detected;                 % 4 - All detected units (active + inactive)
                            avg_spike_frequency;                  % 5
                            electrode_firing_rate;                % 6
                            unit_fr_mean;                         % 7
                            unit_fr_std;                          % 8
                            unit_fr_median;                       % 9
                            unit_fr_min;                          % 10
                            unit_fr_max;                          % 11
                            unit_sc_mean;                         % 12
                            unit_sc_std;                          % 13
                            unit_sc_median;                       % 14
                            unit_sc_min;                          % 15
                            unit_sc_max;                          % 16
                            bursting_units;                       % 17
                            percentage_bursting_units;            % 18
                            avg_bursts_per_unit;                  % 19
                            avg_burst_duration_electrode;         % 20
                            avg_spikes_per_burst_electrode;       % 21
                            avg_inter_burst_interval_electrode;   % 22
                            electrode_burst_frequency;            % 23
                            maxTime;                              % 24
                            electrode_activity_span;              % 25
                            electrode_activity_percentage         % 26
                        ];
                    end
                end
            end
        end
    end
    
    T_unitMeasurement = table(unitMeasurements);
    % Convert data matrices to tables
    if isempty(unitData)
        % Create empty table with correct variable names if no units pass the filter
        T_unit = array2table(zeros(length(unitMeasurements), 0), 'RowNames', unitMeasurements);
    else
        % Create table with unit data
        T_unit = array2table(unitData, 'RowNames', unitMeasurements, 'VariableNames', unitNames);
    end
    T_unit = [T_unitMeasurement T_unit];

    T_electrodeMeasurement = table(electrodeMeasurements);
    if isempty(electrodeData)
        % Create empty table with correct variable names if no electrodes
        T_electrode = array2table(zeros(length(electrodeMeasurements), 0), 'RowNames', electrodeMeasurements);
    else
        % Create table with electrode data
        T_electrode = array2table(electrodeData, 'RowNames', electrodeMeasurements, 'VariableNames', electrodeNames);
    end
    T_electrode = [T_electrodeMeasurement T_electrode];
end