function [outlier_WFs, initial_idx_list] = remove_outliers_before_HDT(idx,spikes,params)
% This function cleans initial spectral clustering before HDT
% outlier WFs from highly active electrodes would be recollected 
flag_tooManyOutlierRemoved = 0;
outlier_WFs = [];

initial_idx_list = {};
% if total WFs > 10,000: perform 3-std outlier-removal at each time step
if length(idx)>10000
    for k_cluster=1:max(idx)
        index_of_idx = find(idx==k_cluster);
        A = spikes(:,idx==k_cluster);
        num_spikes_before_remove_outliers = numel(index_of_idx);
        for num = 1:size(A,1)
            [B,TFoutlier] = rmoutliers(A(num,:),'mean','ThresholdFactor',params.std_cutoff);
            outlier_locs = find(TFoutlier==1);
            A(:,outlier_locs) = [];
            outlier_WFs = [outlier_WFs; index_of_idx(outlier_locs)];
            index_of_idx(outlier_locs) = [];
        end
        number_of_spikes = numel(index_of_idx);
        initial_idx_list{k_cluster} = index_of_idx;
        if (1-number_of_spikes/num_spikes_before_remove_outliers) >= params.flag_threshold
            flag_tooManyOutlierRemoved = 1;
        end
    end
end

% otherwise, perform 3-std outlier-removal once 
if flag_tooManyOutlierRemoved == 0
    initial_idx_list = {};
    for k_cluster=1:max(idx)
        index_of_idx = find(idx==k_cluster);
        A = spikes(:,idx==k_cluster);
        A_RMSE = calculate_WF_RMSE(A);
        [B,TFoutlier] = rmoutliers(A_RMSE,'mean','ThresholdFactor',params.std_cutoff);
        outlier_locs = find(TFoutlier==1);
        A(:,outlier_locs) = [];
        index_of_idx(outlier_locs) = [];
        initial_idx_list{k_cluster} = index_of_idx;
    end
end

end