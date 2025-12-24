function [alignedWaveforms, alignedClusterIndices, timeLags] = align_and_trim_spikes(logFile,waveforms, clusterIndices, mergeGroups)
% ALIGN_AND_TRIM_SPIKES - Align spikes by peak and trim to same length
%
% Within each merge group:
%   1. Find the cluster where the peak occurs earliest (reference)
%   2. Align other clusters by trimming from the front to match peak position
% For all clusters:
%   - Merged clusters: trim from front (for alignment) and end (for length)
%   - Non-merged clusters: trim only from end
%   - All clusters end up with the same final length
%
% INPUTS:
%   waveforms - Matrix of spike waveforms (timepoints × spikes)
%   clusterIndices - Cell array containing spike indices for each cluster
%   mergeGroups - Cell array of cluster indices to be merged (from template_comparison)
%
% OUTPUTS:
%   alignedWaveforms - Peak-aligned and trimmed waveforms
%   alignedClusterIndices - Updated cluster indices after alignment
%   timeLags - Structure containing lag information for each merge group

% Initialize outputs
fprintf(logFile,'\nAlignment started!\n');
numClusters = length(clusterIndices);
timeLags = struct('mergeGroup', {}, 'clusterIndices', {}, 'lags', {}, 'referenceIdx', {});

% Check if waveforms are too short for alignment
MIN_WAVEFORM_LENGTH = 20;  % Minimum waveform length to attempt alignment
originalLength = size(waveforms, 1);

if originalLength < MIN_WAVEFORM_LENGTH
    fprintf(logFile,'WARNING: Waveforms too short (%d samples < %d minimum). Skipping alignment.\n', ...
        originalLength, MIN_WAVEFORM_LENGTH);
    % Return original waveforms without alignment
    alignedWaveforms = waveforms;
    alignedClusterIndices = clusterIndices;
    return;
end

% Get all clusters involved in merging
allMergeClusters = [];
if ~isempty(mergeGroups)
    for g = 1:length(mergeGroups)
        allMergeClusters = [allMergeClusters, mergeGroups{g}];
    end
    allMergeClusters = unique(allMergeClusters);
end

% Track which clusters are being merged
isMergedCluster = false(1, numClusters);
isMergedCluster(allMergeClusters) = true;

% Store trim info for each cluster
clusterTrimFront = zeros(1, numClusters);  % Trim from front for each cluster
clusterTrimEnd = zeros(1, numClusters);  % Trim from end for each cluster
maxAbsLag = 0;  % Maximum front trim needed across all merge groups
groupMaxAbsLags = [];  % Store max front trim for each group

% Maximum allowed lag (samples with lag > this will be excluded)
MAX_ALLOWED_LAG = 10;

% Track clusters/waveforms to exclude due to extreme lag
excludedClusters = [];

% Process each merge group to find time lags
for g = 1:length(mergeGroups)
    groupClusters = mergeGroups{g};
    numInGroup = length(groupClusters);

    % Calculate mean templates for this group
    groupTemplates = zeros(size(waveforms, 1), numInGroup);
    peakTimes = zeros(1, numInGroup);  % Track peak time for each cluster

    for i = 1:numInGroup
        clusterIdx = groupClusters(i);
        groupTemplates(:, i) = mean(waveforms(:, clusterIndices{clusterIdx}), 2);

        % Find peak time (index of minimum value for negative spikes)
        [~, peakTimes(i)] = min(groupTemplates(:, i));
    end

    % Use cluster with earliest peak as reference
    [~, refTemplateIdx] = min(peakTimes);
    referenceIdx = groupClusters(refTemplateIdx);
    referenceTemplate = groupTemplates(:, refTemplateIdx);

    % Find the earliest peak time in this group (for alignment)
    minPeakTime = min(peakTimes);

    % Find time lag for each cluster relative to earliest peak
    groupLags = zeros(1, numInGroup);

    fprintf(logFile,'Merge Group %d:\n', g);

    % Check lags and exclude clusters with excessive lag
    for i = 1:numInGroup
        clusterIdx = groupClusters(i);

        % Lag is the difference from earliest peak
        lag = peakTimes(i) - peakTimes(refTemplateIdx);
        groupLags(i) = lag;

        % Check if lag exceeds threshold
        if abs(lag) > MAX_ALLOWED_LAG
            fprintf(logFile,'  Cluster %d: lag = %d samples - EXCLUDED (lag > %d)\n', ...
                clusterIdx, lag, MAX_ALLOWED_LAG);
            excludedClusters = [excludedClusters, clusterIdx];
            % Don't set trim values for excluded clusters
        else
            fprintf(logFile,'  Cluster %d: lag = %d samples\n', clusterIdx, lag);
            % Trim from front to align with earliest peak
            clusterTrimFront(clusterIdx) = peakTimes(i) - minPeakTime;
        end
    end

    % Recalculate maxAbsLag only considering non-excluded clusters
    validPeakTimes = [];
    for i = 1:numInGroup
        clusterIdx = groupClusters(i);
        if ~ismember(clusterIdx, excludedClusters)
            validPeakTimes = [validPeakTimes, peakTimes(i)];
        end
    end

    if ~isempty(validPeakTimes)
        groupMaxTrimFront = max(validPeakTimes) - min(validPeakTimes);
        maxAbsLag = max(maxAbsLag, groupMaxTrimFront);
        groupMaxAbsLags(g) = groupMaxTrimFront;
    else
        groupMaxAbsLags(g) = 0;
    end

    % Store lag information
    timeLags(g).mergeGroup = g;
    timeLags(g).clusterIndices = groupClusters;
    timeLags(g).lags = groupLags;
    timeLags(g).referenceIdx = referenceIdx;
    timeLags(g).excludedClusters = excludedClusters;
end

% Empty excluded clusters (waveforms with excessive lag)
for clusterIdx = excludedClusters
    fprintf(logFile,'Emptying cluster %d due to excessive lag\n', clusterIdx);
    alignedClusterIndices{clusterIdx} = [];
end

% Ensure all clusters have the same final length
% Strategy: Total trim = maxAbsLag, distributed between front and end
%           trimFront (already set based on lag) + trimEnd = maxAbsLag
%           This ensures: (originalLength - trimFront - trimEnd) = finalLength for all
for clusterIdx = 1:numClusters
    if ~isempty(clusterIndices{clusterIdx}) && ~ismember(clusterIdx, excludedClusters)
        % Set end trim so that total trim equals maxAbsLag
        clusterTrimEnd(clusterIdx) = maxAbsLag - clusterTrimFront(clusterIdx);

        fprintf(logFile,'Cluster %d: trim front = %d, trim end = %d (total = %d)\n', ...
            clusterIdx, clusterTrimFront(clusterIdx), clusterTrimEnd(clusterIdx), ...
            clusterTrimFront(clusterIdx) + clusterTrimEnd(clusterIdx));
    end
end

% Cap the maximum trim to prevent excessive waveform reduction
MAX_TRIM_PER_ALIGNMENT = 10;  % Maximum samples to trim from each end
if maxAbsLag > MAX_TRIM_PER_ALIGNMENT
    fprintf(logFile,'WARNING: Calculated trim (%d samples) exceeds maximum (%d samples). Capping trim amount.\n', ...
        maxAbsLag, MAX_TRIM_PER_ALIGNMENT);
    maxAbsLag = MAX_TRIM_PER_ALIGNMENT;

    % Recalculate trim amounts with capped maxAbsLag
    for clusterIdx = 1:numClusters
        if ~isempty(clusterIndices{clusterIdx}) && ~ismember(clusterIdx, excludedClusters)
            % Cap individual cluster trims
            clusterTrimFront(clusterIdx) = min(clusterTrimFront(clusterIdx), maxAbsLag);
            clusterTrimEnd(clusterIdx) = maxAbsLag - clusterTrimFront(clusterIdx);
        end
    end
end

% Determine final waveform length after trimming
originalLength = size(waveforms, 1);
% All clusters trim total of maxAbsLag (distributed between front and end)
% Final length = original - maxAbsLag
finalLength = originalLength - maxAbsLag;

% Validate that we have enough waveform length
if finalLength <= 0
    error('Alignment error: Cannot trim %d samples from %d-sample waveforms (would result in %d samples). Peak lag of %d samples is too large for current waveform length.', ...
        maxAbsLag, originalLength, finalLength, maxAbsLag);
end

fprintf(logFile,'\nAlignment Summary:\n');
fprintf(logFile,'Original waveform length: %d samples\n', originalLength);
fprintf(logFile,'Maximum lag (peak span): %d samples\n', maxAbsLag);
fprintf(logFile,'Total trim per cluster: %d samples\n', maxAbsLag);
fprintf(logFile,'Final waveform length: %d samples\n', finalLength);
fprintf(logFile,'All clusters trim: %d from end\n', maxAbsLag);
if ~isempty(excludedClusters)
    fprintf(logFile,'Excluded clusters (lag > %d): %d clusters\n', MAX_ALLOWED_LAG, length(excludedClusters));
end

% Create aligned and trimmed waveforms
alignedWaveforms = zeros(finalLength, size(waveforms, 2));
alignedClusterIndices = clusterIndices;  % Indices stay the same, just waveforms change

% Track which spikes are in clusters
allClusteredSpikes = [];
for clusterIdx = 1:numClusters
    if ~isempty(clusterIndices{clusterIdx})
        allClusteredSpikes = [allClusteredSpikes; clusterIndices{clusterIdx}(:)];
    end
end

% Process clustered spikes
for clusterIdx = 1:numClusters
    if isempty(clusterIndices{clusterIdx}) || ismember(clusterIdx, excludedClusters)
        continue;
    end

    spikeIndices = clusterIndices{clusterIdx};
    trimFront = clusterTrimFront(clusterIdx);
    trimEnd = clusterTrimEnd(clusterIdx);

    % Ensure spikeIndices is a column vector for proper iteration
    if size(spikeIndices, 1) < size(spikeIndices, 2)
        spikeIndices = spikeIndices';
    end

    for spikeIdx = spikeIndices'
        originalWaveform = waveforms(:, spikeIdx);

        % Trim trimFront from front, trimEnd from end
        startIdx = trimFront + 1;
        endIdx = originalLength - trimEnd;

        % Ensure the trimmed waveform is a column vector
        trimmedWaveform = originalWaveform(startIdx:endIdx);
        if size(trimmedWaveform, 2) > size(trimmedWaveform, 1)
            trimmedWaveform = trimmedWaveform';
        end

        alignedWaveforms(:, spikeIdx) = trimmedWaveform;
    end
end

% Process unclustered spikes (trim only from end)
totalSpikes = size(waveforms, 2);
unclusteredSpikes = setdiff(1:totalSpikes, allClusteredSpikes');

if ~isempty(unclusteredSpikes)
    fprintf(logFile,'Processing %d unclustered spikes (trim end only)...\n', length(unclusteredSpikes));
    for spikeIdx = unclusteredSpikes
        originalWaveform = waveforms(:, spikeIdx);

        % Trim only from end (no front trim for unclustered spikes)
        startIdx = 1;
        endIdx = originalLength - maxAbsLag;

        % Ensure the trimmed waveform is a column vector
        trimmedWaveform = originalWaveform(startIdx:endIdx);
        if size(trimmedWaveform, 2) > size(trimmedWaveform, 1)
            trimmedWaveform = trimmedWaveform';
        end

        alignedWaveforms(:, spikeIdx) = trimmedWaveform;
    end
end

fprintf(logFile,'\nAlignment complete!\n');
fprintf(logFile,'Clusters merged: %d\n', length(allMergeClusters));
fprintf(logFile,'Clusters not merged: %d\n', sum(~isMergedCluster));

end
