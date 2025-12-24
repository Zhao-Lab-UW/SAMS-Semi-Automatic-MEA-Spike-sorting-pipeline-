function [T, T_electrode, ...
                        raster_raw, sorting_results, num_electrode, total_num_cell, ...
                        num_total_electrodes_with_spikes, num_inactive_units, num_total_detected_units, ...
                        fr_inactive_units, fr_active_units, possibleMUAList, flagTooMuchList] = ...
                        update_results(flg_too_much_threshold, electrode_results, electrode_unit_stats, T, T_electrode, ...
                        raster_raw, ...
                        sorting_results, num_electrode, total_num_cell, num_well, i, j, m, n, ...
                        num_total_electrodes_with_spikes, num_inactive_units, num_total_detected_units, ...
                        fr_inactive_units, fr_active_units, possibleMUAList, flagTooMuchList, pptx, ...
                        maxTime, fileName, current_spikes, current_times, current_score)
% UPDATE_RESULTS - Update results tables and collect statistics
%
% This function updates the results tables and statistical data structures
% based on the electrode processing results
%
% INPUTS:
%   electrode_results - Structure with electrode results
%   electrode_stats - Structure with electrode statistics
%   T - Table of unit measurements
%   T_electrode - Table of electrode measurements
%   raster_raw - Cell array for storing spike times
%   sorting_results - Cell array for storing sorting results
%   num_electrode, total_num_cell - Counters
%   num_well, i, j, m, n - Indices
%   num_total_electrodes_with_spikes, num_inactive_units, num_total_detected_units - Counters
%   fr_inactive_units, fr_active_units - Firing rate arrays
%   possibleMUAList, flagTooMuchList - Lists for problematic electrodes
%   pptx - PowerPoint object
%   maxTime - Maximum recording time
%   maxSpikeAmp, minSpikeAmp - Amplitude range for plotting
%   fileName - Current file name
%   current_spikes - Spike waveforms for current electrode
%   current_times - Spike times for current electrode
%   current_score - PCA scores for current electrode

% OUTPUTS:
%   Updated versions of input variables




    % Current electrode name
    currElectrodeName = [char(i+'A'-1), num2str(j), '_', num2str(m), num2str(n)];
    % Recalculate PCA on aligned/trimmed spikes (waveforms have been trimmed during alignment)
    [~, score, ~] = pca(current_spikes');
    
    % Add to MUA or Flag list if needed
    if electrode_results.possibleMUA == 1
        possibleMUAList{end+1} = currElectrodeName;
    end
    


    unsorted_WFs = [];
    % Process active units
    new_idx_list = electrode_results.new_idx_list;
    if ~isempty(new_idx_list)
        % ... (rest of existing code remains the same) ...
        
        % Update PowerPoint if available
        if ~isempty(pptx)
            % Add a slide for this electrode
            pptx.addSlide('Master', 1, 'Layout', '5 unit');
            % Plot each unit
            sorted_WFs = [];
            T_network = [];
            for unit_ = 1:length(new_idx_list)
                if unit_ ==6

                    pptx.addSlide('Master', 1, 'Layout', '5 unit');
                end
                
                % Track sorted waveforms
                sorted_WFs = [sorted_WFs; new_idx_list{unit_}];
                
                % Create figure for this unit
                f = figure('Visible', 'off'); 
                hold on;
                plot(1:size(current_spikes, 1), zeros(1, size(current_spikes, 1)), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 5);
                
                % Plot waveforms with different colors
                colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'];
                plot(current_spikes(:, new_idx_list{unit_}), colors(rem(unit_, 5) + 1));
                
                % Set title and axis
                title(['unit ', num2str(unit_)]);
                % axis([0 40 minSpikeAmp maxSpikeAmp]);
                hold off;

                % Add to PowerPoint
                pptx.addPicture(f, 'Position', ['Picture Placeholder ', num2str(6 + 2 * (mod(unit_-1, 5) + 1))], 'Scale', 'max');
                close(f);
                
                % Collect spike times for network analysis
                T_network = [T_network current_times(1, new_idx_list{unit_})];
                
                % Update counter
                total_num_cell = total_num_cell + 1;
            end
            
            % Store raster data for network analysis
            raster_raw{num_electrode, 1, num_well} = sort(T_network); 
            raster_raw{num_electrode, 2, num_well} = [i, j, m, n, unit_];
            num_electrode = num_electrode + 1;
            
            % Find unsorted waveforms
            unsorted_WFs = setdiff(1:size(current_spikes, 2), sorted_WFs);


            if length(unsorted_WFs)/size(current_spikes, 2) >= flg_too_much_threshold
                flagTooMuchList{end+1} = currElectrodeName;
            end

            % Store sorting results
            sorting_results{i, j, m, n, 1} = new_idx_list;
            sorting_results{i, j, m, n, 2} = unsorted_WFs;
            
            % Add table to PowerPoint
            pptx_tableData = { ...
                'All Valid WFs: ', size(current_spikes, 2); ...
                'Unsorted WFs: ', length(unsorted_WFs); };
            pptx.addTable(pptx_tableData, 'Position', 'Table Placeholder 11', ...
                'Vert', 'middle', 'Horiz', 'center', 'EdgeColor', 'w');
            
            % Plot PCA figure (using recalculated score from aligned/trimmed spikes)
            F_PCA = figure('visible', 'off');
            hold on;
            plot(score(unsorted_WFs, 1), score(unsorted_WFs, 2), '.', 'Color', [0.5, 0.5, 0.5]);
            for unit_ = 1:length(new_idx_list)
                plot(score(new_idx_list{unit_}, 1), score(new_idx_list{unit_}, 2), '.', 'Color', colors(rem(unit_, 5) + 1));
            end
            title('PCA'); 
            xlabel('PC1'); 
            ylabel('PC2');
            hold off;
            % Add to PowerPoint
            pptx.addPicture(F_PCA, 'Position', 'Picture Placeholder 18', 'Scale', 'max');
            close(F_PCA);
            
            % Add title
            pptx.addTextbox([fileName, ' Electrode ', currElectrodeName], 'Position', 'Title 1');
        end
        
        % Calculate comprehensive electrode statistics
        electrode_comprehensive_stats = calculate_electrode_statistics(electrode_results, electrode_unit_stats, current_spikes, current_times, maxTime);
        
        % Update electrode statistics table with comprehensive data
        electrode_stats_array = [
            electrode_comprehensive_stats.total_spikes;
            electrode_comprehensive_stats.num_excluded_spikes;
            electrode_comprehensive_stats.unsorted_firing_rate;
            electrode_comprehensive_stats.num_units_detected;
            electrode_comprehensive_stats.avg_spike_frequency;
            electrode_comprehensive_stats.electrode_firing_rate;
            electrode_comprehensive_stats.unit_firing_rates_mean;
            electrode_comprehensive_stats.unit_firing_rates_std;
            electrode_comprehensive_stats.unit_firing_rates_median;
            electrode_comprehensive_stats.unit_firing_rates_min;
            electrode_comprehensive_stats.unit_firing_rates_max;
            electrode_comprehensive_stats.unit_spike_counts_mean;
            electrode_comprehensive_stats.unit_spike_counts_std;
            electrode_comprehensive_stats.unit_spike_counts_median;
            electrode_comprehensive_stats.unit_spike_counts_min;
            electrode_comprehensive_stats.unit_spike_counts_max;
            electrode_comprehensive_stats.num_bursting_units;
            electrode_comprehensive_stats.percentage_bursting_units;
            electrode_comprehensive_stats.avg_bursts_per_unit;
            electrode_comprehensive_stats.avg_burst_duration_electrode;
            electrode_comprehensive_stats.avg_spikes_per_burst_electrode;
            electrode_comprehensive_stats.avg_inter_burst_interval_electrode;
            electrode_comprehensive_stats.electrode_burst_frequency;
            electrode_comprehensive_stats.total_recording_duration;
            electrode_comprehensive_stats.electrode_activity_span;
            electrode_comprehensive_stats.electrode_activity_percentage
        ];
        
        Tleft_electrode = table(electrode_stats_array);
        Tleft_electrode.Properties.VariableNames(1) = {currElectrodeName};
        T_electrode = [T_electrode, Tleft_electrode];

        % Update unit statistics table
        for unit_ = 1:length(electrode_unit_stats.active_units)
            stats = electrode_unit_stats.active_units(unit_).stats;
            Tleft = table(stats);
            var_name_unit = [currElectrodeName, '_', num2str(unit_)];
            Tleft.Properties.VariableNames(1) = {var_name_unit};

            T = [T, Tleft];
            
        end
        
        % Update counters
        num_total_electrodes_with_spikes = num_total_electrodes_with_spikes + 1;
        num_total_detected_units = num_total_detected_units + length(new_idx_list) + length(electrode_unit_stats.inactive_units);
        num_inactive_units = num_inactive_units + length(electrode_unit_stats.inactive_units);
        
        % Update firing rate lists
        fr_active_units = [fr_active_units, electrode_unit_stats.active_units_fr];
        fr_inactive_units = [fr_inactive_units, electrode_unit_stats.inactive_units_fr];
    end
end