function [dipStatistic, pValue, lowerBound, upperBound] = HartigansDipSignifTest(distribution, bootstrapSamples, varargin)
% HARTIGANSDIPSIGNIFTEST - Calculates Hartigan's dip statistic and its significance
%
% This function tests for multimodality in a distribution using Hartigan's dip test.
% It calculates the dip statistic and its p-value using bootstrap sampling.
%
% INPUTS:
%   distribution - Vector of sample values to test
%   bootstrapSamples - Number of bootstrap samples for significance testing
%   varargin - Optional parameters for plotting:
%     'plot', figureID - Enable plotting with specified figure ID
%
% OUTPUTS:
%   dipStatistic - The dip statistic (measure of multimodality)
%   pValue - Significance of the dip statistic
%   lowerBound - Lower bound of the distribution's modal interval
%   upperBound - Upper bound of the distribution's modal interval

    % Parse optional plotting parameters
    shouldPlot = false;
    figureID = 6;
    paramIdx = 1;
    
    while paramIdx < nargin - 1
        if strncmpi(varargin{paramIdx}, 'plot', 4)
            shouldPlot = true;
            
            % Check if figure ID is provided
            if isnumeric(varargin{paramIdx+1})
                figureID = varargin{paramIdx+1};
                paramIdx = paramIdx + 1;
            else
                figureID = 1;
            end
        end
        paramIdx = paramIdx + 1;
    end
    
    % Calculate dip statistic for the empirical distribution
    [dipStatistic, lowerBound, upperBound, ifault, gcm, lcm, mn, mj] = HartigansDipTest(distribution);
    sampleSize = length(distribution);
    
    % Calculate bootstrap samples for uniform distribution
    bootstrapDipValues = zeros(bootstrapSamples, 1);
    
    for sampleIdx = 1:bootstrapSamples
        % Generate uniform random sample of the same size
        uniformSample = sort(unifrnd(0, 1, 1, sampleSize));
        
        % Calculate dip statistic for uniform sample
        [bootstrapDipValues(sampleIdx)] = HartigansDipTest(uniformSample);
    end
    
    % Sort bootstrap values and calculate p-value
    bootstrapDipValues = sort(bootstrapDipValues);
    pValue = sum(dipStatistic < bootstrapDipValues) / bootstrapSamples;
    
    % Plot bootstrap distribution and test statistic if requested
    if shouldPlot
        % figure(figureID);
        figure();
        
        % Create histogram of bootstrap dip values
        [histCounts, histBins] = hist(bootstrapDipValues);
        bar(histBins, histCounts, 'k');
        hold on;
        
        % Plot the dip statistic line
        plot([dipStatistic, dipStatistic], [0, max(histCounts)*1.1], 'r:');
        
        title('Bootstrap Distribution of Dip Statistic');
        xlabel('Dip Value');
        ylabel('Frequency');
        legend('Bootstrap Distribution', 'Observed Dip');
    end
end