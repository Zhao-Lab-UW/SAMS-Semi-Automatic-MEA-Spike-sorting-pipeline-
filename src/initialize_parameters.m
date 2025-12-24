function params = initialize_parameters(threshold_to_merge, refractoryT, ...
    network_participation_threshold, min_spikes_E, max_ISI_E, ...
    min_spikes_N, max_ISI_N, cutoff_frequency, flag_threshold, std_cutoff,...
        DTWfigure_check, SortingResults_check, overlap_threshold)
    % Convert string parameters to numeric
    params.threshold_to_merge = str2double(threshold_to_merge);
    params.refractoryT = str2double(refractoryT) / 1000; % Convert ms to seconds
    
    % Store other parameters
    params.network_participation_threshold = network_participation_threshold / 100; % Convert to fraction
    params.min_spikes_E = min_spikes_E;
    params.max_ISI_E = max_ISI_E;
    params.min_spikes_N = min_spikes_N;
    params.max_ISI_N = max_ISI_N;
    params.maximum_inter_spike_interval = max_ISI_E / 1000; % Convert ms to seconds
    params.minimum_spike_per_burst = min_spikes_E;
    params.cutoff_frequency = cutoff_frequency;
    params.flag_threshold = flag_threshold;
    params.std_cutoff = std_cutoff;
    params.violation_threshold = params.refractoryT; % threshold for refractory violations (in seconds)
    params.DTWfigure_check = DTWfigure_check;
    params.SortingResults_check = SortingResults_check;
    params.overlap_threshold = overlap_threshold;
end