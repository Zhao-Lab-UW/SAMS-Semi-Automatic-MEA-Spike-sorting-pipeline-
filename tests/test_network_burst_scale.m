function tests = test_network_burst_scale
% TEST_NETWORK_BURST_SCALE - Network burst analysis of a long recording stays within a few GB
%
% Same size as the A1_42 report: 4 electrodes x 22,616 s at 12.5 kHz (282.7 million samples).
% The pre-fix code needed a 9 GB dense matrix plus copies for this. Each test takes well under
% a minute and measures this MATLAB process's peak memory, so run the file in a fresh MATLAB
% session (Windows only). A later test's peak can only be overestimated, never hidden.
%
% Run from the repository root:
%   results = runtests('tests/test_network_burst_scale.m')

    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testsFolder = fileparts(mfilename('fullpath'));
    testCase.TestData.originalPath = path;
    addpath(fullfile(testsFolder, '..', 'src'));
end

function teardownOnce(testCase)
    path(testCase.TestData.originalPath);
end

function testShellyRecordingLength(testCase)
    % 22,616,000 one-ms bins, the length in the A1_42 report
    verifyLongRecording(testCase, 22616);
end

function testPrimeBinCount(testCase)
    % 22,616,003 bins is prime, the slowest and most memory-hungry FFT length unless padded
    verifyLongRecording(testCase, 22616.002001);
end

function verifyLongRecording(testCase, maxTime)
    fs = 12500;
    numElectrodes = 4;
    burstTimes = 5:10:maxTime;  % a network burst every 10 s

    % 2 Hz background plus 20 spikes within 100 ms in every burst, per electrode
    randomStream = RandStream('mt19937ar', 'Seed', 1);
    raster_raw = cell(numElectrodes, 2, 1);
    for electrodeNum = 1:numElectrodes
        background = rand(randomStream, 1, round(2 * maxTime)) * maxTime;
        bursts = burstTimes + rand(randomStream, 20, numel(burstTimes)) * 0.1;
        raster_raw{electrodeNum, 1, 1} = sort([background, bursts(:)']);
        raster_raw{electrodeNum, 2, 1} = [1, 1, 4, electrodeNum, 1];
    end

    outputFolder = tempname;
    mkdir(outputFolder);
    testCase.addTeardown(@() rmdir(outputFolder, 's'));

    process = System.Diagnostics.Process.GetCurrentProcess();
    process.Refresh();
    memoryBefore = double(process.WorkingSet64);
    tic;
    evalc(['get_network_burst_info(raster_raw, maxTime, fs, 0.35, 5, 100, 50, 100, ', ...
        'outputFolder, {})']);
    elapsedSeconds = toc;
    process.Refresh();
    peakIncreaseGB = (double(process.PeakWorkingSet64) - memoryBefore) / 1e9;
    fprintf('Long recording (%.6f s): %.0f s, peak memory increase at most %.2f GB\n', ...
        maxTime, elapsedSeconds, peakIncreaseGB);

    sheet = readtable(fullfile(outputFolder, 'spike_sorting.xlsx'), 'Sheet', 'network burst', ...
        'VariableNamingRule', 'preserve');
    metric = @(name) sheet{strcmp(sheet{:, 1}, name), 2};
    testCase.verifyEqual(metric('number of network bursts'), numel(burstTimes));
    testCase.verifyEqual(metric('Number of Bursting Electrodes'), numElectrodes);
    testCase.verifyGreaterThan(metric('Synchrony index'), 0);
    testCase.verifyLessThan(metric('Synchrony index'), 1);
    testCase.verifyLessThan(peakIncreaseGB, 3, 'Network burst analysis should need only a few GB');
    testCase.verifyLessThan(elapsedSeconds, 60);
end
