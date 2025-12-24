function [networkBurstInfo] = get_network_spike_participation(samplingRate, spikeMatrix, electrodeIndices, networkThreshold, maxNetworkISI, minNetworkSpikes)
% GET_NETWORK_SPIKE_PARTICIPATION_FIXED - Detect network bursts with correct burst boundary detection
%
% This function identifies synchronized bursting across multiple electrodes
% based on network participation threshold and burst criteria.
%
% INPUTS:
%   samplingRate - Recording sampling rate in Hz
%   spikeMatrix - Binary spike matrix (electrodes × timepoints)
%   electrodeIndices - Indices of electrodes in recording
%   networkThreshold - Minimum fraction of electrodes that must participate
%   maxNetworkISI - Maximum inter-spike interval for network burst (ms)
%   minNetworkSpikes - Minimum spikes required for a network burst
%
% OUTPUTS:
%   networkBurstInfo - Matrix containing network burst information:
%       Row 1: Starting timepoint of each burst
%       Row 2: Number of electrodes participating in each burst
%       Row 3: Duration of each burst (seconds)
%       Row 4: Mean ISI within each burst (seconds)
%       Row 5: Number of spikes per network burst
%       Row 6: Number of firing cells in each burst

    % Convert maxNetworkISI from ms to timepoints
    maxNetworkISI_timepoints = maxNetworkISI * samplingRate / 1000;
% Ensure electrodeIndices is a row vector
    electrodeIndices = electrodeIndices(:);
    % Find timepoints with any firing
    totalFiringByTimepoint = sum(spikeMatrix);
    firingTimepoints = find(totalFiringByTimepoint > 0);

    % Find explicit burst groups based on gap size
    burstGroups = {};

    if ~isempty(firingTimepoints)
        currentGroup = 1;  % Start with first timepoint index

        for idx = 2:length(firingTimepoints)
            gap = firingTimepoints(idx) - firingTimepoints(idx-1);

            if gap < maxNetworkISI_timepoints
                % Small gap - add to current burst group
                currentGroup = [currentGroup, idx];
            else
                % Large gap - save current burst and start new one
                if length(currentGroup) >= 2  % At least 2 spikes for a burst
                    burstGroups{end+1} = currentGroup;
                end
                currentGroup = idx;  % Start new potential burst
            end
        end

        % Add the last group
        if length(currentGroup) >= 2
            burstGroups{end+1} = currentGroup;
        elseif length(currentGroup) == 1 && length(firingTimepoints) == 1
            % Special case: only one spike total - not a burst
            burstGroups = {};
        end
    end

    % ===== PROCESS EACH BURST GROUP =====
    % Initialize network burst statistics
    burstCount = 0;
    burstStartTimepoints = [];
    electrodesPerBurst = [];
    burstDurations = [];
    meanISIWithinBurst = [];
    spikesPerBurst = [];
    cellsPerBurst = [];
    % Process each identified burst group
    for groupIdx = 1:length(burstGroups)
        burstIndices = burstGroups{groupIdx};
        burstTimepoints = firingTimepoints(burstIndices);

        % Collect all participating electrodes and cells in this burst
        firingElectrodeIndices = [];
        firingCellIndices = [];

        for idx = burstIndices
            currentTimepoint = firingTimepoints(idx);
            currentFiringCells = find(spikeMatrix(:, currentTimepoint) == 1);

            if ~isempty(currentFiringCells)
                firingElectrodeIndices = [firingElectrodeIndices; electrodeIndices(currentFiringCells)];
                firingCellIndices = [firingCellIndices; currentFiringCells];
            end
        end

        % Get unique electrode and cell indices
        uniqueElectrodeIndices = unique(firingElectrodeIndices);
        uniqueCellIndices = unique(firingCellIndices);

        % Calculate total spikes in this burst
        totalSpikes = sum(spikeMatrix(:, burstTimepoints), 'all');

        % Check if this qualifies as a network burst
        totalElectrodes = length(unique(electrodeIndices));

        if (length(uniqueElectrodeIndices) > networkThreshold * totalElectrodes) && ...
           (totalSpikes >= minNetworkSpikes)

            % Valid network burst detected
            burstCount = burstCount + 1;

            % Store burst information
            burstStartTimepoints(burstCount) = burstTimepoints(1);
            electrodesPerBurst(burstCount) = length(uniqueElectrodeIndices);

            % Duration in timepoints (will convert to seconds later)
            if length(burstTimepoints) > 1
                burstDurations(burstCount) = burstTimepoints(end) - burstTimepoints(1);
            else
                burstDurations(burstCount) = 0;
            end

            % Calculate mean ISI within burst
            if length(burstTimepoints) > 1
                meanISIWithinBurst(burstCount) = mean(diff(burstTimepoints));
            else
                meanISIWithinBurst(burstCount) = 0;
            end

            % Store spike and cell counts
            spikesPerBurst(burstCount) = totalSpikes;
            cellsPerBurst(burstCount) = length(uniqueCellIndices);
        end
    end

    % Convert durations and ISIs from timepoints to seconds
    if burstCount > 0
        burstDurations = burstDurations / samplingRate;
        meanISIWithinBurst = meanISIWithinBurst / samplingRate;
        % Combine network burst information
        networkBurstInfo = [burstStartTimepoints; electrodesPerBurst; burstDurations; ...
                             meanISIWithinBurst; spikesPerBurst; cellsPerBurst];
    else
        % No bursts found, return empty matrix with correct dimensions
        networkBurstInfo = zeros(6, 0);
    end

    % Replace any NaN values with zeros
    networkBurstInfo(isnan(networkBurstInfo)) = 0;
end
