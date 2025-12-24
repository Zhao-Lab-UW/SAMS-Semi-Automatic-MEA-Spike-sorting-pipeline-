function [electrode_results, electrode_unit_stats,spikes] = process_electrode(logFile, spikes, Times, i, j, m, n, params, maxTime, maxSpikeAmp, minSpikeAmp, score)
% PROCESS_ELECTRODE - Process a single electrode
% change on 10-27-2025: outlier-removal on 4 windows
% INPUTS:
%   spikes - Spike waveforms matrix (timepoints × spikes)
%   Times - Spike timing vector
%   i, j, m, n - Electrode indices
%   params - Parameter structure
%   maxTime - Maximum recording time
%   maxSpikeAmp, minSpikeAmp - Amplitude range
%   score - PCA scores (optional)
%
% OUTPUTS:
%   electrode_results - Structure with electrode results (includes DTW_table)
%   electrode_stats - Structure with electrode statistics

% Initialize results
electrode_results = struct();
electrode_unit_stats = struct();
electrode_results.new_idx_list = {};
electrode_results.HDT_flag = 0;
electrode_results.DTW_flag = 0;
electrode_results.possibleMUA = 0;
electrode_results.flagTooMuch = 0;
electrode_results.flag_outliers = 0;

% Current electrode name
currElectrodeName = [char(i+'A'-1), num2str(j), '_', num2str(m), num2str(n)];
electrode_results.name = currElectrodeName;

% Check if enough spikes
numSpikes = size(spikes, 2);
if numSpikes < 20
    % fprintf(logFile, 'Too few spikes (%d) for processing.  Skip electrode.\n', numSpikes);
    return;
end

% Perform PCA if not provided
if nargin < 11 || isempty(score)
    try
        [~, score, ~] = pca(spikes');
    catch ME
        fprintf('Error in PCA: %s\n', ME.message);
        return;
    end
end

% Check if enough dimensions in score
if size(score, 2) < 5
    % fprintf(logFile, 'Not enough dimensions in PCA score. Skip electrode.\n');
    return;
end

% Initialize spectral clustering parameters
myfunc = @(X, K)(spectralcluster(X, K));
klist = 1:5; % Number of clusters to try

% Evaluate optimal number of clusters
try
    eva = evalclusters(score(:,1:2), myfunc, "DaviesBouldin", 'klist', klist);

    if isnan(eva.OptimalK)
        % fprintf(logFile, 'Could not determine optimal number of clusters. Skip electrode.\n');
        return;
        % K = 2;
    else
        % fprintf('Optimal number of clusters: %d\n', eva.OptimalK);
        K = eva.OptimalK;
    end
catch ME
    fprintf(logFile, 'Error in cluster evaluation: %s\n', ME.message);
    fprintf(logFile, 'Using default K=2.\n');
    K = 2;
end

% Perform spectral clustering
try
    if K > 1
        [idx, ~] = spectralcluster(score(:,1:2), K);
    else
        idx = ones(size(score, 1), 1);
    end
catch ME
    fprintf(logFile, 'Error in spectral clustering: %s\n', ME.message);
    idx = ones(size(score, 1), 1);
end

% clean initial clustering before HDT
[outlier_WFs, initial_idx_list] = remove_outliers_before_HDT(idx,spikes,params);

%% for main unit list
% Fix under-sorting for main unit list
[initial_idx_list, HDT_flag] = fix_undersorting(initial_idx_list, spikes, Times, score, params.refractoryT, params);
electrode_results.HDT_flag = HDT_flag;

% Remove empty cells from list
initial_idx_list = initial_idx_list(~cellfun('isempty', initial_idx_list));

if isempty(initial_idx_list)
    fprintf(logFile, 'Empty clustering idx.\n');
    return;
end

% Perform template comparison and merge similar units

% Initialize variables that may not be set in all code paths
merge_outlier_list = {};
initial_idx_list_unsorted_WFs = {};

[dist, merge, distMerge, distNotMerge] = template_comparison(spikes, initial_idx_list, params.threshold_to_merge, params.overlap_threshold);

% Align and trim waveforms for merge groups BEFORE merging
if ~isempty(merge)
    fprintf(logFile, 'Aligning waveforms for first merge groups...\n');
    [spikes, initial_idx_list, ~] = align_and_trim_spikes(logFile, spikes, initial_idx_list, merge);
end

if ~isempty(outlier_WFs) % try classifying outlier_WFs is existed
    [merge_outlier_list,initial_idx_list_unsorted_WFs,spikes] = template_comparison_outlier(logFile,spikes,initial_idx_list,outlier_WFs,myfunc,klist,score,params);
end

% Merge units based on DTW distance
if ~isempty(merge)
    electrode_results.DTW_flag = 1;
    for ii = 1:size(merge, 2)
        target_idx = merge{ii}(1);  % First cluster is the target
        for jj = 2:size(merge{ii}, 2)
            source_idx = merge{ii}(jj);
            % Combine spike indices
            initial_idx_list{target_idx} = sort([
                initial_idx_list{target_idx};
                initial_idx_list{source_idx}
                ]);

            % Clear the source cluster
            initial_idx_list{source_idx} = [];
        end
    end
end

initial_idx_list = [initial_idx_list,initial_idx_list_unsorted_WFs];

% Remove empty cells from list
initial_idx_list = initial_idx_list(~cellfun('isempty', initial_idx_list));

% second merge based on DTW distance
[dist, merge, distMerge, distNotMerge] = template_comparison(spikes, initial_idx_list, params.threshold_to_merge, params.overlap_threshold);

% Align and trim waveforms for second merge groups BEFORE merging
if ~isempty(merge)
    fprintf(logFile, 'Aligning waveforms for second merge groups \n');
    [spikes, initial_idx_list, ~] = align_and_trim_spikes(logFile, spikes, initial_idx_list, merge);
end

if ~isempty(merge)
    electrode_results.DTW_flag = 1;
    for ii = 1:size(merge, 2)
        target_idx = merge{ii}(1);  % First cluster is the target
        for jj = 2:size(merge{ii}, 2)
            source_idx = merge{ii}(jj);
            % Combine spike indices
            initial_idx_list{target_idx} = sort([
                initial_idx_list{target_idx};
                initial_idx_list{source_idx}
                ]);

            % Clear the source cluster
            initial_idx_list{source_idx} = [];
        end
    end
end
% Remove empty cells from list
initial_idx_list = initial_idx_list(~cellfun('isempty', initial_idx_list));

%% process final cluster list

% Process final clusters
[new_idx_list, electrode_unit_stats] = process_final_clusters(initial_idx_list, spikes, Times, params, maxTime, minSpikeAmp, maxSpikeAmp);
electrode_results.new_idx_list = new_idx_list;

% Check for possible multi-unit activity
if length(new_idx_list) > 1
    [dist, merge, distMerge, distNotMerge] = template_comparison(spikes, initial_idx_list, params.threshold_to_merge, params.overlap_threshold);
    if ~isempty(find(distNotMerge < params.threshold_to_merge, 1))
        electrode_results.possibleMUA = 1;
        electrode_results.possibleMUAName = currElectrodeName;
    end
end

% Check for inconsistency in flags
if (electrode_results.HDT_flag == 1) &&(electrode_results.DTW_flag==1)
    fprintf(logFile, 'Flag consistency in electrode %s: HDT=%d, DTW=%d\n', currElectrodeName, electrode_results.HDT_flag, electrode_results.DTW_flag);
end

% Check for inconsistency in flags
if (electrode_results.HDT_flag == 1) &&(electrode_results.DTW_flag==0)
    fprintf(logFile, 'Flag INconsistency in electrode %s: HDT=%d, DTW=%d\n', currElectrodeName, electrode_results.HDT_flag, electrode_results.DTW_flag);
end
end