function [meanDistance, mergeGroups, mergeDistances, nonMergeDistances] = template_comparison(waveforms, clusterIndices, mergeThreshold,overlap_threshold)
% TEMPLATE_COMPARISON - Compare spike templates to identify similar units for merging
%
% This function calculates distances between mean templates of spike clusters
% and identifies which clusters should be merged based on similarity.
%
% INPUTS:
%   waveforms - Matrix of spike waveforms (timepoints × spikes)
%   clusterIndices - Cell array containing spike indices for each cluster
%   mergeThreshold - Distance threshold below which clusters are merged
%
% OUTPUTS:
%   meanDistance - Mean distance between all template pairs
%   mergeGroups - Cell array of cluster indices to be merged
%   mergeDistances - Distances between templates that should be merged
%   nonMergeDistances - Distances between templates that should not be merged

% Calculate mean template for each cluster
templateCount = 0;
meanTemplates = [];
validClusterIndices = [];

clusterIndices = clusterIndices(~cellfun('isempty',clusterIndices));
for clusterIdx = 1:length(clusterIndices)
    templateCount = templateCount + 1;
    validClusterIndices(templateCount) = clusterIdx;
    meanTemplates(:, templateCount) = mean(waveforms(:, clusterIndices{clusterIdx}), 2);
end

% Store original templates before normalization for peak/trough comparison
originalTemplates = meanTemplates;

% Normalize templates for better comparison
templateMin = min(meanTemplates(:));
templateMax = max(meanTemplates(:));
% save('results.mat')
normalizedTemplates = (meanTemplates - templateMin) / (templateMax - templateMin);

% Calculate distances between templates
pairwiseDistances = [];
mergeCandidate = [];
mergeCount = 0;
nonMergeCount = 0;
mergeDistances = [];
nonMergeDistances = [];

% Compare each template pair
for template1_idx = 1:templateCount-1
    for template2_idx = template1_idx+1:templateCount
        % Calculate DTW distance using normalized templates
        currentDistance = dtw(normalizedTemplates(:, template1_idx), normalizedTemplates(:, template2_idx));
        pairwiseDistances(end+1) = currentDistance;

        % Check if templates have similar shape features using normalized templates
        template1_normalized = normalizedTemplates(:, template1_idx);
        template2_normalized = normalizedTemplates(:, template2_idx);

        % Find peaks and troughs from original templates
        peak1 = max(template1_normalized);
        trough1 = min(template1_normalized);
        peak2 = max(template2_normalized);
        trough2 = min(template2_normalized);

        % Calculate height from original templates
        height1 = peak1 - trough1;
        height2 = peak2 - trough2;

        % Use minimum height for strictest comparison (must satisfy both templates' tolerances)
        minHeight = min(height1, height2);

        % Define threshold for peak/trough matching
        shapeSimilarityThreshold = 1-overlap_threshold; %using params.overlap_threshold=0.7,


        peakMatch = 0; troughMatch = 0;
        if peak1<(peak2+height1*shapeSimilarityThreshold ) && peak1>(peak2-height1*shapeSimilarityThreshold )
            peakMatch = 1;
        end
        if trough1<(trough2+height1*shapeSimilarityThreshold ) && trough1>(trough2-height1*shapeSimilarityThreshold )
            troughMatch = 1;
        end

        % % Check if peaks match (using minimum height for strictest tolerance)
        % peakMatch = abs(peak1 - peak2) < (minHeight * shapeSimilarityThreshold);
        %
        % % Check if troughs match (using minimum height for strictest tolerance)
        % troughMatch = abs(trough1 - trough2) < (minHeight * shapeSimilarityThreshold);

        % Determine if templates should be merged
        shouldMerge = (currentDistance < mergeThreshold) && peakMatch && troughMatch;
        if shouldMerge
            mergeCount = mergeCount + 1;
            mergeCandidate(mergeCount, :) = [template1_idx, template2_idx];
            mergeDistances(mergeCount) = currentDistance;
        else
            nonMergeCount = nonMergeCount + 1;
            nonMergeDistances(nonMergeCount) = currentDistance;
        end
    end
end

% Calculate mean distance
if ~isempty(pairwiseDistances)
    meanDistance = mean(pairwiseDistances);
else
    meanDistance = 999; % Default high value if no comparisons
end

% Group merge candidates into connected components
mergeGroups = {};

for m = 1:size(mergeCandidate,1)
    if m==1
        mergeGroups{1} = mergeCandidate(m,:);
    else
        merged = 0;
        n=1;
        while n < size(mergeGroups,2)+1
            % search if first already exists
            if ~isempty(intersect(mergeCandidate(m,:),mergeGroups{n}))
                mergeGroups{n} = unique([mergeGroups{n},mergeCandidate(m,:)]);
                merged = 1;
                n = n+1;
                break;
            end
            n = n+1;
        end
        if merged==0
            mergeGroups{size(mergeGroups,2)+1} = mergeCandidate(m,:);
        end
    end
end

end