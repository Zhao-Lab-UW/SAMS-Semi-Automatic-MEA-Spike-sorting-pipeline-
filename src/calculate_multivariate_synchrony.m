function synchronizationIndex = calculate_multivariate_synchrony(signalMatrix)
% CALCULATE_MULTIVARIATE_SYNCHRONY - Calculate synchronization index across multiple signals
%
% This function calculates the Kuramoto order parameter to measure synchronization
% across multiple neural signals based on their phase relationships.
%
% INPUTS:
%   signalMatrix - Matrix where each column represents a signal time series
%
% OUTPUTS:
%   synchronizationIndex - Synchronization index (0-1), where 1 indicates
%                          perfect synchronization across all signals

    % Get dimensions
    [numSamples, numSignals] = size(signalMatrix);
    
    % Check if there are enough signals to calculate synchrony
    if numSignals < 2
        warning('Need at least 2 signals to calculate synchronization. Returning 0.');
        synchronizationIndex = 0;
        return;
    end
    
    % Calculate the Hilbert transform to get instantaneous phase
    try
        % Apply Hilbert transform to get complex analytic signal
        analyticSignal = hilbert(signalMatrix);
        
        % Extract the instantaneous phase
        instantaneousPhase = angle(analyticSignal);
        
        % Calculate the complex order parameter at each time point
        % e^(i*θ) for each signal, then average across signals
        orderParameter = sum(exp(1i * instantaneousPhase), 2) / numSignals;
        
        % Calculate the synchronization index (mean amplitude of order parameter)
        synchronizationIndex = mean(abs(orderParameter));
    catch ME
        warning('Error calculating synchronization: %s', ME.message);
        synchronizationIndex = 0;
    end
    
    % Ensure output is in valid range
    synchronizationIndex = max(0, min(1, synchronizationIndex));
end