function [networkBurstInfo, isBurstingElectrode] = get_network_spike_participation(samplingRate, spikeSamples, electrodeIndices, networkThreshold, maxNetworkISI, minNetworkSpikes)
% GET_NETWORK_SPIKE_PARTICIPATION - Detect network bursts with correct burst boundary detection
%
% This function identifies synchronized bursting across multiple electrodes
% based on network participation threshold and burst criteria.
%
% Spikes are passed as sample indices per electrode instead of an electrodes x samples
% matrix, so memory scales with the number of spikes, not the recording length.
%
% INPUTS:
%   samplingRate - Recording sampling rate in Hz
%   spikeSamples - Cell array, one entry per electrode: sorted unique sample indices of
%                  its spikes (may be empty)
%   electrodeIndices - Indices of electrodes in recording, one unique index per
%                      spikeSamples entry
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
%   isBurstingElectrode - Logical column, one entry per electrode: true if the electrode
%       fires in at least one network burst

    % Convert maxNetworkISI from ms to timepoints
    maxNetworkISI_timepoints = maxNetworkISI * samplingRate / 1000;
    % Ensure electrodeIndices is a column vector
    electrodeIndices = electrodeIndices(:);
    numElectrodes = numel(spikeSamples);
    isBurstingElectrode = false(numElectrodes, 1);

    % List every spike with its electrode, then find timepoints with any firing
    spikeCounts = cellfun(@numel, spikeSamples(:));
    allSpikeSamples = zeros(sum(spikeCounts), 1);
    spikeElectrode = reshape(repelem((1:numElectrodes)', spikeCounts), [], 1);
    spikeOffset = 0;
    for electrodeNum = 1:numElectrodes
        allSpikeSamples(spikeOffset + (1:spikeCounts(electrodeNum))) = spikeSamples{electrodeNum};
        spikeOffset = spikeOffset + spikeCounts(electrodeNum);
    end
    [firingTimepoints, ~, spikeTimepoint] = unique(allSpikeSamples);

    % Find explicit burst groups based on gap size: runs of firing timepoints whose
    % gaps are all shorter than the max ISI, keeping runs of at least 2 timepoints
    if isempty(firingTimepoints)
        groupStarts = zeros(0, 1);
        groupEnds = zeros(0, 1);
    else
        isGap = ~(diff(firingTimepoints) < maxNetworkISI_timepoints);
        groupStarts = [1; find(isGap) + 1];
        groupEnds = [find(isGap); numel(firingTimepoints)];
        isBurstGroup = groupEnds - groupStarts >= 1;
        groupStarts = groupStarts(isBurstGroup);
        groupEnds = groupEnds(isBurstGroup);
    end
    numGroups = numel(groupStarts);

    % Assign every spike to its burst group (0 = outside any group)
    timepointGroup = zeros(numel(firingTimepoints), 1);
    for groupIdx = 1:numGroups
        timepointGroup(groupStarts(groupIdx):groupEnds(groupIdx)) = groupIdx;
    end
    spikeGroup = timepointGroup(spikeTimepoint);
    inGroup = spikeGroup > 0;
    % Keep these as columns: with a single spike, indexing would give 0x0
    groupOfSpike = reshape(spikeGroup(inGroup), [], 1);
    electrodeOfSpike = reshape(spikeElectrode(inGroup), [], 1);

    % Spikes and participating electrodes per group. Each raster row is one electrode,
    % so the number of firing cells equals the number of participating electrodes
    spikesPerGroup = accumarray(groupOfSpike, 1, [numGroups, 1]);
    groupElectrodePairs = unique([groupOfSpike, electrodeOfSpike], 'rows');
    electrodesPerGroup = accumarray(groupElectrodePairs(:, 1), 1, [numGroups, 1]);

    % Check which groups qualify as network bursts
    totalElectrodes = length(unique(electrodeIndices));
    isNetworkBurst = (electrodesPerGroup > networkThreshold * totalElectrodes) & ...
        (spikesPerGroup >= minNetworkSpikes);
    burstGroupIndices = find(isNetworkBurst)';
    burstCount = numel(burstGroupIndices);

    % ===== PROCESS EACH NETWORK BURST =====
    burstStartTimepoints = zeros(1, burstCount);
    burstDurations = zeros(1, burstCount);
    meanISIWithinBurst = zeros(1, burstCount);
    for burstIdx = 1:burstCount
        groupIdx = burstGroupIndices(burstIdx);
        burstTimepoints = firingTimepoints(groupStarts(groupIdx):groupEnds(groupIdx));

        % Duration in timepoints (will convert to seconds later)
        burstStartTimepoints(burstIdx) = burstTimepoints(1);
        burstDurations(burstIdx) = burstTimepoints(end) - burstTimepoints(1);

        % Calculate mean ISI within burst
        meanISIWithinBurst(burstIdx) = mean(diff(burstTimepoints));
    end
    electrodesPerBurst = electrodesPerGroup(burstGroupIndices)';
    spikesPerBurst = spikesPerGroup(burstGroupIndices)';
    cellsPerBurst = electrodesPerBurst;

    % Electrodes that fire in at least one network burst
    isBurstPair = isNetworkBurst(groupElectrodePairs(:, 1));
    isBurstingElectrode(groupElectrodePairs(isBurstPair, 2)) = true;

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
