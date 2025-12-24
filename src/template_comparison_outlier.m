function [merge_outlier_list,initial_idx_list_unsorted_WFs,spikes_aligned] = template_comparison_outlier(logFile,spikes,idx_list,outlier_WFs,myfunc,klist,score,params)
% TEMPLATE_COMPARISON_OUTLIER - Compare outlier waveforms to existing templates
% This function loops through outlier waveforms and compares them to existing
% unit templates to determine if they should be merged or clustered separately.
%
% Inputs:
%   spikes       - Waveform matrix (N x M), N: time frames, M: total number of waveforms
%   idx_list     - Cell array of indices for existing units
%   outlier_WFs  - Vector of indices for outlier waveforms
%   myfunc       - Clustering function handle
%   klist        - Range of k values for clustering evaluation
%   score        - PCA score matrix
%   params       - Parameter struct containing threshold_to_merge and std_cutoff
%
% Outputs:
%   merge_outlier_list            - Cell array mapping outliers to existing units
%   initial_idx_list_unsorted_WFs - Cell array of new clusters from unmerged outliers

std_cutoff = params.std_cutoff;
% Initialize variables
n = length(idx_list);
mean_spikes = [];

% Calculate mean templates for existing units
valid_unit_count = 0;
for k_cluster = 1:n
    if ~isempty(idx_list{k_cluster})
        valid_unit_count = valid_unit_count + 1;
        mean_spikes(:, valid_unit_count) = mean(spikes(:, idx_list{k_cluster}), 2);
    end
end

% Normalize mean_spikes for comparison
min_val = min(mean_spikes(:));
max_val = max(mean_spikes(:));
if max_val - min_val == 0
    mean_spikes_norm = mean_spikes;
else
    mean_spikes_norm = (mean_spikes - min_val) / (max_val - min_val);
end

% Normalize all spikes
min_val = min(spikes(:));
max_val = max(spikes(:));
if max_val - min_val == 0
    spikes_norm = spikes;
else
    spikes_norm = (spikes - min_val) / (max_val - min_val);
end


% Compare outlier waveforms to existing templates
merge_list = [];
count_merge = 0;
dists_outlier_WFs = zeros(length(outlier_WFs), valid_unit_count);

for i = 1:length(outlier_WFs)  % Loop through outlier waveforms
    outlier_WF_to_compare = spikes_norm(:, outlier_WFs(i));
    dists = zeros(1, valid_unit_count);
    matches = zeros(1, valid_unit_count);

    % Calculate distances and check peak/trough matches
    for j = 1:valid_unit_count  % Loop through existing template waveforms
        % Calculate DTW distance
        dists(j) = dtw(outlier_WF_to_compare, mean_spikes_norm(:, j));

        % Check peak and trough matching
        peakI = max(outlier_WF_to_compare);
        troughI = min(outlier_WF_to_compare);
        heightI = peakI - troughI;

        peakJ = max(mean_spikes_norm(:, j));
        troughJ = min(mean_spikes_norm(:, j));

        thresholdRange = 1-params.overlap_threshold;  % 30% tolerance for peak/trough matching

        % Check if peaks match within threshold
        peakMatch = (peakI < (peakJ + heightI * thresholdRange)) && ...
            (peakI > (peakJ - heightI * thresholdRange));

        % Check if troughs match within threshold
        troughMatch = (troughI < (troughJ + heightI * thresholdRange)) && ...
            (troughI > (troughJ - heightI * thresholdRange));

        % Both peak and trough must match
        if peakMatch && troughMatch
            matches(j) = 1;
        else
            matches(j) = 0;
        end
    end

    % Store distances for this outlier
    dists_outlier_WFs(i, :) = dists;

    % Find the closest template
    [min_dist, idx_min] = min(dists);

    % Check if this outlier should be merged with the closest template
    if matches(idx_min) == 1
        count_merge = count_merge + 1;
        merge_list(count_merge, :) = [i, idx_min];
    end
end

% Organize merge results into cell array format
merge_outlier_list = cell(n, 1);
merged_outliers = [];
spikes_aligned = spikes;  % Initialize with original spikes

if ~isempty(merge_list)
    % Build alignment groups: each main unit with its matched outliers
    merge_groups_for_alignment = {};
    temp_idx_list = idx_list;  % Start with existing units

    for i = 1:n
        merge_indices = find(merge_list(:, 2) == i);
        if ~isempty(merge_indices)
            outlier_indices = outlier_WFs(merge_list(merge_indices, 1));
            merge_outlier_list{i} = outlier_indices;
            merged_outliers = [merged_outliers; outlier_indices];

            % Add outliers as a temporary cluster
            temp_cluster_idx = length(temp_idx_list) + 1;
            temp_idx_list{temp_cluster_idx} = outlier_indices';

            % Create merge group: [main_unit_index, outlier_cluster_index]
            merge_groups_for_alignment{end+1} = [i, temp_cluster_idx];
        end
    end
end

% Identify outliers that weren't merged
unsorted_WFs = setdiff(outlier_WFs, merged_outliers);
initial_idx_list_unsorted_WFs = {};

eva_unsorted_WFs = evalclusters(score(unsorted_WFs,1:2),myfunc,"DaviesBouldin",'klist',klist);
if ~isnan(eva_unsorted_WFs.OptimalK)
    [idx_unsorted_WFs,C] = spectralcluster(score(unsorted_WFs,1:2), eva_unsorted_WFs.OptimalK);
    % [idx_unsorted_WFs,C] = spectralcluster(peak_valleys, eva_unsorted_WFs.OptimalK);
    initial_idx_list_unsorted_WFs = {};
    for k_file=1:max(idx_unsorted_WFs)
        % remove outliers (using aligned spikes)
        index_of_idx = unsorted_WFs(find(idx_unsorted_WFs==k_file));
        A = spikes_aligned(:,unsorted_WFs(idx_unsorted_WFs==k_file));
        for num = 1:size(A,1)
            [B,TFoutlier] = rmoutliers(A(num,:),'mean','ThresholdFactor',std_cutoff);
            % disp('a')
            locs = find(TFoutlier==1);
            A(:,locs) = [];
            index_of_idx(locs) = [];
        end
        initial_idx_list_unsorted_WFs{k_file} = index_of_idx;
    end
end


end
