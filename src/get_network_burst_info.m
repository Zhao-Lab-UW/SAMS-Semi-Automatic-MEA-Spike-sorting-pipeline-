function get_network_burst_info(raster_raw, maxTime, samplingRate, networkParticipationThreshold, ...
    minSpikesElectrode, maxISIElectrode, minSpikesNetwork, maxISINetwork, outputFolder, sorting_results)
% GET_NETWORK_BURST_INFO - Analyze network bursting activity across electrodes
%
% This function detects and analyzes network bursts by identifying synchronized 
% bursting activity across multiple electrodes in the recording.
%
% INPUTS:
%   raster_raw - Cell array containing spike times for each electrode (FIXED: parameter name)
%   maxTime - Total recording duration in seconds
%   samplingRate - Recording sampling frequency in Hz
%   networkParticipationThreshold - Fraction of electrodes required to participate in network burst (0-1)
%   minSpikesElectrode - Minimum number of spikes per burst for electrode bursts
%   maxISIElectrode - Maximum inter-spike interval for electrode bursts (ms)
%   minSpikesNetwork - Minimum number of spikes per network burst
%   maxISINetwork - Maximum inter-spike interval for network bursts (ms)
%   outputFolder - Path to save results
%   sorting_results - Cell array containing spike sorting results
%
% OUTPUTS:
%   Excel file with network burst metrics
%   MAT file with burst information


    % Define network burst measurements
    networkMeasurements = {'number of network bursts',...
        'network burst frequency (Hz)',...
        'Network burst duration -avg (s)',...
        'Network burst duration -std (s)',...
        'Number of spikes per network burst - avg',...
        'Number of spikes per network burst - std',...
        'Number of electrodes participation in burst -avg',...
        'Number of electrodes participation in burst -std',...
        'Number of spikes per network burst per channel -avg',...
        'Number of spikes per network burst per channel -std',...
        'Number of spikes per unit per network burst -avg',...
        'Number of spikes per unit per network burst -std',...
        'Network burst percentage',...
        'Network IBI coefficient of variation',...
        'Network ISI coefficient of variation',...
        'Network normalized duration IQR',...
        'Synchrony index',...
        'Number of Bursting Electrodes',...
        'Number of Active Electrodes'};

    % Initialize results table
    resultsTable = table(networkMeasurements');
    
    % Initialize storage for burst information
    burst_info_all = {};

    % Check if raster_raw is empty
    if isempty(raster_raw) || size(raster_raw, 3) == 0
        fprintf('No data available for network burst analysis.\n');
        % Save empty results
        writetable(resultsTable, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'network burst');
        return;
    end
    
    % Process each well
    for wellIndex = 1:size(raster_raw, 3)
        % Skip if no data for this well
        if isempty(raster_raw{1, 2, wellIndex})
            continue;
        end
        
        % Initialize electrode counter and data structures
        electrodeCount = 1;
        spikeTrain = zeros(size(raster_raw, 1), ceil(maxTime * samplingRate)); % electrodes × time
        electrodeIndices = []; % Stores electrode indices
        activeElectrodeCount = 0;
        
        % Process each electrode 
        validElectrodeIndices = [];
        for electrodeNum = 1:size(raster_raw, 1)
            % Skip empty electrodes
            if ~isempty(raster_raw{electrodeNum, 1, wellIndex})
                % Count active electrodes
                activeElectrodeCount = activeElectrodeCount + 1;
                validElectrodeIndices = [validElectrodeIndices, electrodeNum];
                
                % Get spike times and convert to timepoints
                spikeTimes = round(raster_raw{electrodeNum, 1, wellIndex} * samplingRate);
                spikeTimes(spikeTimes == 0) = 1;  % Ensure valid indices
                spikeTimes(spikeTimes > size(spikeTrain, 2)) = []; % Remove out-of-bounds indices
                
                % Create binary spike train
                spikeTrain(electrodeCount, spikeTimes) = 1;
                
                % Track electrode indices - simplified version
                electrodeIndices(electrodeCount, :) = electrodeNum;
                
                electrodeCount = electrodeCount + 1;
            end
        end
        
        % Trim unused rows from spikeTrain
        spikeTrain = spikeTrain(1:electrodeCount-1, :);
        electrodeIndices = electrodeIndices(1:electrodeCount-1, :);
        
        % Only analyze if there are multiple electrodes
        if size(spikeTrain, 1) > 1
            % Detect network bursts
            networkBurstInfo = get_network_spike_participation(samplingRate, ...
                spikeTrain, electrodeIndices, networkParticipationThreshold, ...
                maxISINetwork, minSpikesNetwork);
            
            burst_info_all{wellIndex, 1} = networkBurstInfo;
            
            % Calculate network metrics
            if isempty(networkBurstInfo)
                fprintf('No network bursts detected for well %d\n', wellIndex);
                wellMetrics = zeros(19, 1);
            else
                % Count electrodes that participate in bursts
                electrodesWithBursts = 0;
                for e = 1:size(spikeTrain, 1)
                    for b = 1:size(networkBurstInfo, 2)
                        burstStart = networkBurstInfo(1, b);
                        burstEnd = burstStart + (networkBurstInfo(3, b) * samplingRate);
                        burstEnd = min(burstEnd, size(spikeTrain, 2));
                        if any(spikeTrain(e, burstStart:burstEnd))
                            electrodesWithBursts = electrodesWithBursts + 1;
                            break; % This electrode has at least one burst
                        end
                    end
                end
                burstingElectrodeCount = electrodesWithBursts;
                
                % Basic network burst counts and rates
                numNetworkBursts = size(networkBurstInfo, 2);
                networkBurstFrequency = numNetworkBursts / maxTime;
                
                % Duration statistics
                networkBurstDurationAvg = mean(networkBurstInfo(3, :));
                networkBurstDurationStd = std(networkBurstInfo(3, :));
                
                % Spike count statistics
                spikesPerNetworkBurstAvg = mean(networkBurstInfo(5, :));
                spikesPerNetworkBurstStd = std(networkBurstInfo(5, :));
                
                % Electrode participation statistics
                electrodesInBurstAvg = mean(networkBurstInfo(2, :));
                electrodesInBurstStd = std(networkBurstInfo(2, :));
                
                % Per-channel spike statistics with division by zero check
                if any(networkBurstInfo(2, :) > 0)
                    validIndices = networkBurstInfo(2, :) > 0;
                    spikesPerChannelAvg = mean(networkBurstInfo(5, validIndices) ./ networkBurstInfo(2, validIndices));
                    spikesPerChannelStd = std(networkBurstInfo(5, validIndices) ./ networkBurstInfo(2, validIndices));
                else
                    spikesPerChannelAvg = 0;
                    spikesPerChannelStd = 0;
                end
                
                % Per-unit spike statistics with division by zero check
                if any(networkBurstInfo(6, :) > 0)
                    validIndices = networkBurstInfo(6, :) > 0;
                    spikesPerUnitAvg = mean(networkBurstInfo(5, validIndices) ./ networkBurstInfo(6, validIndices));
                    spikesPerUnitStd = std(networkBurstInfo(5, validIndices) ./ networkBurstInfo(6, validIndices));
                else
                    spikesPerUnitAvg = 0;
                    spikesPerUnitStd = 0;
                end
                
                % Overall activity metrics (with division by zero check)
                totalSpikes = sum(spikeTrain, 'all');
                if totalSpikes > 0
                    networkBurstPercentage = sum(networkBurstInfo(5, :)) / totalSpikes * 100;
                else
                    networkBurstPercentage = 0;
                end
                
                % Interval variability metrics with division by zero checks
                if size(networkBurstInfo, 2) > 1
                    burstIntervals = diff(networkBurstInfo(1, :));
                    if ~isempty(burstIntervals) && mean(burstIntervals) > 0
                        networkIBICoV = std(burstIntervals) / mean(burstIntervals);
                    else
                        networkIBICoV = 0;
                    end
                else
                    networkIBICoV = 0;
                end
                
                if mean(networkBurstInfo(4, :)) > 0
                    networkISICoV = std(networkBurstInfo(4, :)) / mean(networkBurstInfo(4, :));
                else
                    networkISICoV = 0;
                end
                
                if median(networkBurstInfo(3, :)) > 0
                    networkNormalizedDurationIQR = iqr(networkBurstInfo(3, :)) / median(networkBurstInfo(3, :));
                else
                    networkNormalizedDurationIQR = 0;
                end
                
                % Calculate synchrony index with function check
                synchronyIndex = calculate_multivariate_synchrony(spikeTrain');

                
                % Compile all metrics
                wellMetrics = [
                    numNetworkBursts;
                    networkBurstFrequency;
                    networkBurstDurationAvg;
                    networkBurstDurationStd;
                    spikesPerNetworkBurstAvg;
                    spikesPerNetworkBurstStd;
                    electrodesInBurstAvg;
                    electrodesInBurstStd;
                    spikesPerChannelAvg;
                    spikesPerChannelStd;
                    spikesPerUnitAvg;
                    spikesPerUnitStd;
                    networkBurstPercentage;
                    networkIBICoV;
                    networkISICoV;
                    networkNormalizedDurationIQR;
                    synchronyIndex;
                    burstingElectrodeCount;
                    activeElectrodeCount
                ];
            end
        else
            % Not enough electrodes for network analysis
            wellMetrics = zeros(19, 1);
        end
        
        % Add well metrics to results table
        if ~isempty(raster_raw{1, 2, wellIndex})
            wellRow = raster_raw{1, 2, wellIndex}(1);
            wellCol = raster_raw{1, 2, wellIndex}(2);
            wellName = [char(wellRow + 'A' - 1), num2str(wellCol)];
            
            wellResultsTable = table(wellMetrics);
            wellResultsTable.Properties.VariableNames(1) = {wellName};
            resultsTable = [resultsTable, wellResultsTable];
        end
        
        % Store well identifier
        burst_info_all{wellIndex, 2} = raster_raw{1, 2, wellIndex};
    end
    
    % Save results to Excel
    writetable(resultsTable, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'network burst');
    
    % Save detailed burst information to MAT file
    save([outputFolder, '\burst_info_all.mat'], 'burst_info_all', 'raster_raw', 'maxTime', 'sorting_results', '-v7.3');
    
    fprintf('Network burst analysis completed: %d wells processed\n', size(raster_raw, 3));
end
