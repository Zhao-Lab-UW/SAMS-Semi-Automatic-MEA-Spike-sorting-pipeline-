function [synchronizationIndex, errorMessage] = calculate_multivariate_synchrony(spikeSamples, numSamples, samplingRate)
% CALCULATE_MULTIVARIATE_SYNCHRONY - Calculate synchronization index across multiple signals
%
% This function calculates the Kuramoto order parameter to measure synchronization
% across multiple neural signals based on their phase relationships.
%
% Each electrode's spikes are counted in 1 ms bins and the instantaneous phase comes from
% the Hilbert transform of those counts, one electrode at a time. Memory scales with the
% recording length in ms rather than in samples. Up to version d4456d3 the Hilbert transform
% ran on the full-rate 0/1 spike train; on the sample recording the two agree within 0.01.
%
% INPUTS:
%   spikeSamples - Cell array, one entry per signal: sorted unique sample indices
%                  (1..numSamples) of its spikes
%   numSamples - Recording length in samples
%   samplingRate - Recording sampling frequency in Hz
%
% OUTPUTS:
%   synchronizationIndex - Synchronization index (0-1), where 1 indicates
%                          perfect synchronization across all signals.
%                          NaN if it could not be calculated
%   errorMessage - Why it could not be calculated ('' otherwise)

    numSignals = numel(spikeSamples);
    errorMessage = '';

    % Check if there are enough signals to calculate synchrony
    if numSignals < 2
        warning('Need at least 2 signals to calculate synchronization. Returning 0.');
        synchronizationIndex = 0;
        return;
    end

    % Calculate the Hilbert transform to get instantaneous phase
    try
        samplesPerBin = samplingRate / 1000;
        numBins = ceil(numSamples / samplesPerBin);

        % The FFT inside hilbert needs about twice the memory and runs several times slower
        % when its length has a large prime factor, so pad with a few empty bins (under
        % 50 ms) up to a length whose prime factors are all at most 1000
        fftLength = numBins;
        while max(factor(fftLength)) > 1000
            fftLength = fftLength + 1;
        end

        orderParameterSum = complex(zeros(numBins, 1));
        for signalNum = 1:numSignals
            % Spike counts per 1 ms bin
            spikeBins = floor((spikeSamples{signalNum}(:) - 1) / samplesPerBin) + 1;
            binCounts = accumarray(spikeBins, 1, [numBins, 1]);

            % Apply Hilbert transform to get complex analytic signal, extract the
            % instantaneous phase and add e^(i*phase) to the order parameter sum
            instantaneousPhase = angle(hilbert(binCounts, fftLength));
            orderParameterSum = orderParameterSum + exp(1i * instantaneousPhase(1:numBins));
        end

        % Calculate the synchronization index (mean amplitude of order parameter)
        synchronizationIndex = mean(abs(orderParameterSum)) / numSignals;
    catch ME
        warning('SAMS:synchronyFailed', ...
            'Synchrony index could not be calculated (%s). Reporting NaN.', ME.message);
        synchronizationIndex = NaN;
        errorMessage = ME.message;
        return;
    end

    % Ensure output is in valid range
    synchronizationIndex = max(0, min(1, synchronizationIndex));
end
