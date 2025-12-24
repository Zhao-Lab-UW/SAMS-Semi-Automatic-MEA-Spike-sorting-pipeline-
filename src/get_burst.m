function burstInfo = get_burst(spikeTimes, maxInterSpikeInterval, minSpikesPerBurst)
% GET_BURST - Detect bursts in spike trains based on inter-spike intervals
%
% This function detects bursts in spike times based on maximum inter-spike interval
% and minimum number of spikes per burst criteria.
%
% INPUTS:
%   spikeTimes - Vector of spike times (assumed sorted)
%   maxInterSpikeInterval - Maximum ISI to be considered same burst
%   minSpikesPerBurst - Minimum spikes required for a burst
%
% OUTPUTS:
%   burstInfo - Matrix with burst information:
%       Row 1: Starting spike index of each burst
%       Row 2: Number of spikes in each burst  
%       Row 3: Duration of each burst (seconds)
%       Row 4: Mean ISI within each burst
%       Row 5: Inter-burst interval (time from end of this burst to start of next)
%       Row 6: Std ISI within each burst
%       Row 7: CV of ISI within burst (std/mean)
%       Row 8: Min ISI within each burst
%       Row 9: Max ISI within each burst
%       Row 10: Burst start time (seconds)
%       Row 11: Burst end time (seconds)
%       Row 12: Mean firing rate within burst (Hz)
%       Row 13: Burst surprise index (ratio to baseline firing)
%       Row 14: ISI acceleration (trend of ISIs within burst)

    % Input validation
    if isempty(spikeTimes) || length(spikeTimes) < minSpikesPerBurst
        burstInfo = [];
        return;
    end
    
    % Ensure row vector
    spikeTimes = spikeTimes(:)';
    
    % Pre-allocate arrays (for efficiency)
    maxBursts = ceil(length(spikeTimes) / minSpikesPerBurst);
    burstStartIndices = zeros(1, maxBursts);
    spikesPerBurst = zeros(1, maxBursts);
    burstDurations = zeros(1, maxBursts);
    meanISIWithinBurst = zeros(1, maxBursts);
    stdISIWithinBurst = zeros(1, maxBursts);
    cvISIWithinBurst = zeros(1, maxBursts);
    minISIWithinBurst = zeros(1, maxBursts);
    maxISIWithinBurst = zeros(1, maxBursts);
    burstStartTimes = zeros(1, maxBursts);
    burstEndTimes = zeros(1, maxBursts);
    burstFiringRates = zeros(1, maxBursts);
    isiAcceleration = zeros(1, maxBursts);
    
    burstCount = 0;
    
    % Main burst detection loop
    spikeIdx = 1;
    while spikeIdx <= length(spikeTimes)
        % Start a potential burst at current spike
        burstEndIdx = spikeIdx;
        
        % Extend burst while ISIs are small enough
        while burstEndIdx < length(spikeTimes) && ...
              (spikeTimes(burstEndIdx + 1) - spikeTimes(burstEndIdx)) < maxInterSpikeInterval
            burstEndIdx = burstEndIdx + 1;
        end
        
        % Calculate burst size
        currentBurstSize = burstEndIdx - spikeIdx + 1;
        
        % Check if this qualifies as a burst
        if currentBurstSize >= minSpikesPerBurst
            burstCount = burstCount + 1;
            
            % Basic indices
            burstStartIndices(burstCount) = spikeIdx;
            spikesPerBurst(burstCount) = currentBurstSize;
            
            % Timing
            burstStartTimes(burstCount) = spikeTimes(spikeIdx);
            burstEndTimes(burstCount) = spikeTimes(burstEndIdx);
            burstDurations(burstCount) = burstEndTimes(burstCount) - burstStartTimes(burstCount);
            
            % Firing rate within burst (Hz)
            if burstDurations(burstCount) > 0
                % Use N-1 for rate calculation (number of intervals)
                burstFiringRates(burstCount) = (currentBurstSize - 1) / burstDurations(burstCount);
            else
                % All spikes at same time point
                burstFiringRates(burstCount) = Inf;
            end
            
            % ISI statistics within burst
            if currentBurstSize > 1
                burstISIs = diff(spikeTimes(spikeIdx:burstEndIdx));
                
                % Basic stats
                meanISIWithinBurst(burstCount) = mean(burstISIs);
                stdISIWithinBurst(burstCount) = std(burstISIs);
                minISIWithinBurst(burstCount) = min(burstISIs);
                maxISIWithinBurst(burstCount) = max(burstISIs);
                
                % Coefficient of variation
                if meanISIWithinBurst(burstCount) > 0
                    cvISIWithinBurst(burstCount) = stdISIWithinBurst(burstCount) / ...
                                                   meanISIWithinBurst(burstCount);
                else
                    cvISIWithinBurst(burstCount) = 0;
                end
                
                % ISI acceleration (linear trend in ISIs)
                % Positive = ISIs increasing (burst slowing down)
                % Negative = ISIs decreasing (burst speeding up)
                if length(burstISIs) > 2
                    x = 1:length(burstISIs);
                    p = polyfit(x, burstISIs, 1);
                    isiAcceleration(burstCount) = p(1);  % Slope
                else
                    isiAcceleration(burstCount) = 0;
                end
            else
                % Single spike "burst" (shouldn't happen with minSpikesPerBurst > 1)
                meanISIWithinBurst(burstCount) = 0;
                stdISIWithinBurst(burstCount) = 0;
                cvISIWithinBurst(burstCount) = 0;
                minISIWithinBurst(burstCount) = 0;
                maxISIWithinBurst(burstCount) = 0;
                isiAcceleration(burstCount) = 0;
            end
            
            % Move past this burst
            spikeIdx = burstEndIdx + 1;
        else
            % Not a burst, move to next spike
            spikeIdx = spikeIdx + 1;
        end
    end
    
    % Trim arrays to actual burst count
    if burstCount > 0
        burstStartIndices = burstStartIndices(1:burstCount);
        spikesPerBurst = spikesPerBurst(1:burstCount);
        burstDurations = burstDurations(1:burstCount);
        meanISIWithinBurst = meanISIWithinBurst(1:burstCount);
        stdISIWithinBurst = stdISIWithinBurst(1:burstCount);
        cvISIWithinBurst = cvISIWithinBurst(1:burstCount);
        minISIWithinBurst = minISIWithinBurst(1:burstCount);
        maxISIWithinBurst = maxISIWithinBurst(1:burstCount);
        burstStartTimes = burstStartTimes(1:burstCount);
        burstEndTimes = burstEndTimes(1:burstCount);
        burstFiringRates = burstFiringRates(1:burstCount);
        isiAcceleration = isiAcceleration(1:burstCount);
        
        % Calculate inter-burst intervals
        interBurstIntervals = NaN(1, burstCount);
        if burstCount > 1
            for i = 1:(burstCount - 1)
                interBurstIntervals(i) = burstStartTimes(i + 1) - burstEndTimes(i);
            end
        end
        
        % Calculate burst surprise index (optional - requires baseline)
        % This compares burst firing rate to overall firing rate
        totalDuration = spikeTimes(end) - spikeTimes(1);
        if totalDuration > 0
            overallFiringRate = length(spikeTimes) / totalDuration;
            burstSurpriseIndex = burstFiringRates / overallFiringRate;
        else
            burstSurpriseIndex = ones(1, burstCount);
        end
        
        % Combine all statistics
        burstInfo = [
            burstStartIndices;      % Row 1
            spikesPerBurst;          % Row 2
            burstDurations;          % Row 3
            meanISIWithinBurst;      % Row 4
            interBurstIntervals;     % Row 5
            stdISIWithinBurst;       % Row 6
            cvISIWithinBurst;        % Row 7
            minISIWithinBurst;       % Row 8
            maxISIWithinBurst;       % Row 9
            burstStartTimes;         % Row 10
            burstEndTimes;           % Row 11
            burstFiringRates;        % Row 12
            burstSurpriseIndex;      % Row 13
            isiAcceleration          % Row 14
        ];
    else
        burstInfo = [];
    end
end