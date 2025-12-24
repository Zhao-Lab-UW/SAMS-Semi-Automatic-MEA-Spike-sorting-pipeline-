function [initial_idx_list, HDT_flag] = fix_undersorting(initial_idx_list, spikes, Times, score, refractoryT, params)
    % Fix undersorting using Hartigan's dip test for bimodality
    %
    % INPUTS:
    %   initial_idx_list - Cell array of indices for each cluster after outlier removal
    %   spikes - Spike waveforms matrix (time points × spikes)
    %   Times - Spike timing vector
    %   score - PCA scores for each spike
    %   refractoryT - Refractory period in seconds
    %   params - Parameter structure
    %
    % OUTPUTS:
    %   initial_idx_list - Updated cell array of indices after fixing undersorting
    %   HDT_flag - Flag indicating if Hartigan's dip test found bimodality

    HDT_flag = 0;
    nboot = 500; % Bootstrap sample size for the dip test
    count_fix_undersorting = 0;
    copy_initial_idx_list = initial_idx_list;
    initial_idx_list = cell(1, length(copy_initial_idx_list)*2); % Pre-allocate larger array
    
    for num_unit_initial = 1:length(copy_initial_idx_list)
        try
            if length(copy_initial_idx_list{num_unit_initial}) > 10
                % Get PCA scores for this cluster
                samplePCA1 = score(copy_initial_idx_list{num_unit_initial}, 1)';
                samplePCA2 = score(copy_initial_idx_list{num_unit_initial}, 2)';
                % Perform Hartigan's dip test on the first two PCs
                [~, p_value1, ~, ~] = HartigansDipSignifTest(samplePCA1, nboot);
                [~, p_value2, ~, ~] = HartigansDipSignifTest(samplePCA2, nboot);
                % If either PC shows significant bimodality (p < 0.05)
                if (p_value1 < 0.05) || (p_value2 < 0.05)
                    try
                        % Perform additional clustering on this unit
                        sample2d = score(copy_initial_idx_list{num_unit_initial}, 1:2);

                        % Use evalclusters to determine optimal K
                        klist = 2:5;
                        myfunc = @(X, K) spectralcluster(X, K);
                        eva = evalclusters(sample2d, myfunc, "DaviesBouldin", 'KList', klist);
                        optimal_K = eva.OptimalK;

                        % Perform spectral clustering with optimal K
                        [cluster_idx, ~] = spectralcluster(sample2d, optimal_K);
                    catch ME
                        if contains(ME.message, 'Invalid data type') || contains(ME.message, 'real array')
                            try
                                % Try k-means as fallback
                                optimal_K = 2; % Default to 2 for fallback
                                [cluster_idx, ~] = kmeans(sample2d, optimal_K, 'Replicates', 2);
                                warning('Used k-means fallback');
                            catch
                                % If even k-means fails, skip this iteration
                                warning('Both spectral clustering and k-means failed');
                                continue;
                            end
                        else
                            rethrow(ME);  % Re-throw if it's a different error
                        end
                    end

                    % Split the cluster into K sub-clusters based on optimal K
                    sub_clusters = cell(1, optimal_K);
                    for k = 1:optimal_K
                        sub_clusters{k} = copy_initial_idx_list{num_unit_initial}(cluster_idx == k);
                    end

                    % Validate sub-clusters through pairwise overlap analysis
                    valid_clusters = true(1, optimal_K);
                    merge_groups = cell(1, optimal_K);
                    for k = 1:optimal_K
                        merge_groups{k} = k;
                    end

                    % Check all pairs of sub-clusters for overlap
                    for k1 = 1:optimal_K-1
                        for k2 = k1+1:optimal_K
                            if valid_clusters(k1) && valid_clusters(k2)
                                % Get spike waveforms and times for each sub-cluster
                                spk1 = spikes(:, sub_clusters{k1});
                                t1 = Times(:, sub_clusters{k1});
                                spk2 = spikes(:, sub_clusters{k2});
                                t2 = Times(:, sub_clusters{k2});

                                % Calculate overlap between the two sub-clusters
                                [ovlp1, ovlp2, tdiff] = calculate_overlap(spk1, spk2, t1, t2);
                                ovlp = max(sum(ovlp2 > params.overlap_threshold) / length(ovlp1), sum(ovlp2 > params.overlap_threshold) / length(ovlp2));

                                % Check if the sub-clusters are likely from the same spike
                                sameSpike = (ovlp > params.overlap_threshold) && (tdiff(1) > refractoryT);

                                if sameSpike
                                    % Mark for merging - combine indices
                                    sub_clusters{k1} = [sub_clusters{k1}; sub_clusters{k2}];
                                    valid_clusters(k2) = false;
                                end
                            end
                        end
                    end

                    % Add validated sub-clusters to the output
                    num_valid = sum(valid_clusters);
                    if num_valid > 1
                        % Split detected - add each valid sub-cluster
                        HDT_flag = 1;
                        for k = 1:optimal_K
                            if valid_clusters(k)
                                count_fix_undersorting = count_fix_undersorting + 1;
                                initial_idx_list{count_fix_undersorting} = sub_clusters{k};
                            end
                        end
                    else
                        % All sub-clusters merged back - keep as single unit
                        count_fix_undersorting = count_fix_undersorting + 1;
                        initial_idx_list{count_fix_undersorting} = copy_initial_idx_list{num_unit_initial};
                    end
                else
                    % Keep as a single unit if no bimodality detected
                    count_fix_undersorting = count_fix_undersorting + 1;
                    initial_idx_list{count_fix_undersorting} = copy_initial_idx_list{num_unit_initial};
                end
            end
        catch
            % Keep as is if an error occurs
            fprintf('Error in dip test for unit %d, keeping as is\n', num_unit_initial);
            count_fix_undersorting = count_fix_undersorting + 1;
            initial_idx_list{count_fix_undersorting} = copy_initial_idx_list{num_unit_initial};
        end
    end
    
    % Trim empty cells
    initial_idx_list = initial_idx_list(~cellfun('isempty', initial_idx_list));
end