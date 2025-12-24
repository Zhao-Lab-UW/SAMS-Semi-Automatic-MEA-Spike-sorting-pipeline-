function waveformRMSE = calculate_WF_RMSE(waveforms)
% CALCULATE_WF_RMSE - Calculate root mean square error between spike waveforms and their mean
%
% This function calculates the RMSE between each spike waveform and the mean waveform,
% providing a measure of waveform variability within a cluster.
%
% INPUTS:
%   waveforms - Matrix of spike waveforms (timepoints × spikes)
%
% OUTPUTS:
%   waveformRMSE - Vector of RMSE values for each spike waveform

    % Calculate the mean waveform across all spikes
    meanWaveform = mean(waveforms, 2);
    
    % Calculate residuals for each waveform relative to the mean
    waveformResiduals = waveforms - meanWaveform;
    
    % Square the residuals
    squaredResiduals = waveformResiduals.^2;
    
    % Calculate mean squared error for each waveform
    meanSquaredError = sum(squaredResiduals, 1) / size(waveforms, 1);
    
    % Calculate RMSE by taking the square root
    waveformRMSE = sqrt(meanSquaredError);
end