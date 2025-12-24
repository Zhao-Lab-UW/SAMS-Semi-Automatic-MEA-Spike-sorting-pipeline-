function [overlapRatio1, overlapRatio2, timeDifferences] = calculate_overlap(spikeWaveforms1, spikeWaveforms2, spikeTimes1, spikeTimes2)
% CALCULATE_OVERLAP - Calculate overlap between two sets of spike waveforms
%
% This function removes outliers from each set of waveforms, then calculates
% the amount of overlap in amplitude range between the two sets.
%
% INPUTS:
%   spikeWaveforms1 - First set of spike waveforms (timepoints × spikes)
%   spikeWaveforms2 - Second set of spike waveforms (timepoints × spikes)
%   spikeTimes1 - Timing of first set of spikes
%   spikeTimes2 - Timing of second set of spikes
%
% OUTPUTS:
%   overlapRatio1 - Overlap ratio relative to amplitude range of first set
%   overlapRatio2 - Overlap ratio relative to amplitude range of second set
%   timeDifferences - Sorted time differences between spikes

    % % Remove outliers from first set of waveforms
    % cleanWaveforms1 = remove_waveform_outliers(spikeWaveforms1);
    %
    % % Remove outliers from second set of waveforms
    % cleanWaveforms2 = remove_waveform_outliers(spikeWaveforms2);
    for num = 1:size(spikeWaveforms1,1)
        [B,TFoutlier] = rmoutliers(spikeWaveforms1(num,:),'mean');
        locs = find(TFoutlier==1);
        spikeWaveforms1(:,locs) = [];
        spikeTimes1(locs) = [];  % Also remove corresponding spike times
    end

    for num = 1:size(spikeWaveforms2,1)
        [B,TFoutlier] = rmoutliers(spikeWaveforms2(num,:),'mean');
        locs = find(TFoutlier==1);
        spikeWaveforms2(:,locs) = [];
        spikeTimes2(locs) = [];  % Also remove corresponding spike times
    end
    cleanWaveforms1 = spikeWaveforms1;
    cleanWaveforms2 = spikeWaveforms2;
    % Calculate amplitude ranges for each set
    maxAmplitude1 = max(cleanWaveforms1, [], 2);
    minAmplitude1 = min(cleanWaveforms1, [], 2);

    maxAmplitude2 = max(cleanWaveforms2, [], 2);
    minAmplitude2 = min(cleanWaveforms2, [], 2);
    
    % Calculate overlap
    minOfMaxes = min(maxAmplitude1, maxAmplitude2);
    maxOfMins = max(minAmplitude1, minAmplitude2);
    
    % Calculate overlap ratios
    amplitudeRange1 = maxAmplitude1 - minAmplitude1;
    amplitudeRange2 = maxAmplitude2 - minAmplitude2;
    
    overlapRatio1 = (minOfMaxes - maxOfMins) ./ amplitudeRange1;
    overlapRatio2 = (minOfMaxes - maxOfMins) ./ amplitudeRange2;
    
    % Calculate time differences between spikes
    combinedTimes = sort([spikeTimes1, spikeTimes2]);
    timeDifferences = diff(combinedTimes);
end
