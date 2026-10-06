function tests = test_network_burst
% TEST_NETWORK_BURST - Tests for get_network_burst_info against the frozen pre-fix copy
%
% The pre-fix functions live in tests/legacy as *_legacy.m. Each comparison test runs the
% current and the legacy function on the same input and compares the 'network burst' sheet
% and the saved burst_info_all exactly, except the synchrony index: its definition changed
% from full-rate to 1 ms bins, so it only has to agree within a tolerance.
%
% Run from the repository root:
%   results = runtests('tests/test_network_burst.m')
% The memory test for a long recording is in test_network_burst_scale.m.

    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testsFolder = fileparts(mfilename('fullpath'));
    testCase.TestData.originalPath = path;
    testCase.TestData.repoRoot = fullfile(testsFolder, '..');
    addpath(fullfile(testsFolder, '..', 'src'));
    addpath(fullfile(testsFolder, 'legacy'));
    addpath(fullfile(testsFolder, 'helpers'));
end

function teardownOnce(testCase)
    path(testCase.TestData.originalPath);
end

function setup(testCase)
    testCase.TestData.outputRoot = tempname;
    mkdir(testCase.TestData.outputRoot);
end

function teardown(testCase)
    if exist(testCase.TestData.outputRoot, 'dir')
        rmdir(testCase.TestData.outputRoot, 's');
    end
end

% ---------------------------------------------------------------------------------------------
% Step 1: spike times before the first sample
% ---------------------------------------------------------------------------------------------

function testNegativeTimeCrashesLegacy(testCase)
    % Reproduces Shelly's A1_42 failure: badsubscript at line 88
    params = defaultParams();
    raster_raw = makeBurstyRaster(4, 10, 1);
    raster_raw{1, 1, 1} = [-0.0005, raster_raw{1, 1, 1}];

    outputFolder = makeOutputFolder(testCase, 'legacy');
    caughtError = [];
    try
        evalc(['get_network_burst_info_legacy(raster_raw, params.maxTime, params.fs, params.threshold, ', ...
            'params.minSpikesE, params.maxISIE, params.minSpikesN, params.maxISIN, outputFolder, {})']);
    catch caughtError
    end

    testCase.assertNotEmpty(caughtError, 'The legacy function should fail on a negative spike time');
    testCase.verifyEqual(caughtError.identifier, 'MATLAB:badsubscript');
    testCase.verifyEqual(caughtError.stack(1).name, 'get_network_burst_info_legacy');
    testCase.verifyEqual(caughtError.stack(1).line, 88);
end

function testTimesAtOrBeforeZeroDropped(testCase)
    % Spikes at -0.5 ms and at t = 0 are dropped, as process_final_clusters.m drops them
    params = defaultParams();
    raster_raw = makeBurstyRaster(4, 10, 1);
    rasterEdited = raster_raw;
    rasterEdited{1, 1, 1} = [-0.0005, 0, raster_raw{1, 1, 1}];

    [newOutput, consoleText] = runCurrent(testCase, rasterEdited, params);
    legacyOutput = runLegacy(testCase, raster_raw, params);

    verifyOutputsEqual(testCase, newOutput, legacyOutput);
    expectedLine = ['electrode A1_11: dropped 2 spike time(s) at or before t = 0 (earliest -0.5000 ms) ', ...
        'and 0 non-finite time(s)'];
    testCase.verifySubstring(consoleText, expectedLine);
    testCase.verifySubstring(fileread(fullfile(makeOutputFolder(testCase, 'current'), 'export_log.txt')), expectedLine);
end

function testNonFiniteTimesDropped(testCase)
    params = defaultParams();
    raster_raw = makeBurstyRaster(4, 10, 2);
    rasterNonFinite = raster_raw;
    rasterNonFinite{2, 1, 1} = [raster_raw{2, 1, 1}, NaN, Inf];
    rasterNonFinite{3, 1, 1} = [-Inf, raster_raw{3, 1, 1}];

    [newOutput, consoleText] = runCurrent(testCase, rasterNonFinite, params);
    legacyOutput = runLegacy(testCase, raster_raw, params);

    verifyOutputsEqual(testCase, newOutput, legacyOutput);
    testCase.verifySubstring(consoleText, 'electrode A1_21: dropped 0 spike time(s) at or before t = 0 and 2 non-finite time(s)');
    testCase.verifySubstring(consoleText, 'electrode A1_31: dropped 0 spike time(s) at or before t = 0 and 1 non-finite time(s)');
end

function testColumnVectorTimes(testCase)
    % raster_raw entries are rows today; a column, or even a matrix, must work the same way
    params = defaultParams();
    raster_raw = makeBurstyRaster(4, 10, 3);
    raster_raw{2, 1, 1} = raster_raw{2, 1, 1}(1:2 * floor(end / 2));
    rasterShapes = raster_raw;
    rasterShapes{1, 1, 1} = [-0.0005; -0.0002; raster_raw{1, 1, 1}(:)];
    rasterShapes{2, 1, 1} = reshape(raster_raw{2, 1, 1}, 2, []);

    [newOutput, consoleText] = runCurrent(testCase, rasterShapes, params);
    legacyOutput = runLegacy(testCase, raster_raw, params);

    verifyOutputsEqual(testCase, newOutput, legacyOutput);
    testCase.verifySubstring(consoleText, 'dropped 2 spike time(s) at or before t = 0 (earliest -0.5000 ms)');
end

function testValidDataUnchanged(testCase)
    params = defaultParams();
    raster_raw = makeBurstyRaster(4, 10, 4);

    [newOutput, consoleText] = runCurrent(testCase, raster_raw, params);
    legacyOutput = runLegacy(testCase, raster_raw, params);

    verifyOutputsEqual(testCase, newOutput, legacyOutput);
    testCase.verifyEmpty(strfind(consoleText, 'Network burst: electrode'));
    testCase.verifyFalse(isfile(fullfile(makeOutputFolder(testCase, 'current'), 'export_log.txt')), ...
        'Valid data must not create a log file');
    testCase.verifyGreaterThan(newOutput.sheet{1, 2}, 0, 'Synthetic data should contain network bursts');
end

function testSampleFileUnchanged(testCase)
    % Real 15-minute recording from the repo's sample output
    sampleFile = fullfile(testCase.TestData.repoRoot, 'sample files', 'SAMS Output Files', 'burst_info_all.mat');
    loaded = load(sampleFile, 'raster_raw', 'maxTime', 'sorting_results');
    params = defaultParams();
    params.maxTime = loaded.maxTime;

    [newOutput, consoleText] = runCurrent(testCase, loaded.raster_raw, params, loaded.sorting_results);
    legacyOutput = runLegacy(testCase, loaded.raster_raw, params, loaded.sorting_results);

    synchronyDifference = verifyOutputsEqual(testCase, newOutput, legacyOutput, 'in sample file', 0.015);
    fprintf('Sample file: largest synchrony difference %.4f\n', synchronyDifference);
    testCase.verifyEmpty(strfind(consoleText, 'Network burst: electrode'));
end

function testUnitAnalysisDropsTimesAtOrBeforeZero(testCase)
    % A manual edit that adds unsorted waveforms at or before t = 0 to a unit must give the
    % same unit and electrode statistics as the unit without them. A real spike 1 ms after
    % t = 0 checks that the filter runs before the refractory step.
    fs = 12500;
    waveformLength = 38;
    maxTime = 2;
    randomStream = RandStream('mt19937ar', 'Seed', 5);
    cutoutStarts = [0.001, sort(rand(randomStream, 1, 59) * 1.8 + 0.05), -11 / fs, 0];
    times = (0:waveformLength - 1)' / fs + round(cutoutStarts * fs) / fs;
    voltages = repmat(-sin(linspace(0, pi, waveformLength))' * 50e-6, 1, numel(cutoutStarts)) ...
        + 1e-6 * randn(randomStream, waveformLength, numel(cutoutStarts));
    allData = cell(1, 1, 1, 2);
    allData{1, 1, 1, 2} = MockSpikeData(times, voltages);

    editedResults = cell(1, 1, 1, 2, 2);
    editedResults{1, 1, 1, 2, 1} = {(1:62)'};
    automaticResults = editedResults;
    automaticResults{1, 1, 1, 2, 1} = {(1:60)'};

    outputFolder = makeOutputFolder(testCase, 'current');
    [editedUnits, editedElectrodes] = runUnitAnalysis(editedResults, allData, maxTime, outputFolder);
    [automaticUnits, automaticElectrodes] = runUnitAnalysis(automaticResults, allData, maxTime, outputFolder);

    testCase.verifyEqual(editedUnits, automaticUnits);
    testCase.verifyEqual(editedElectrodes, automaticElectrodes);
end

% ---------------------------------------------------------------------------------------------
% Step 2: sparse spike storage and synchrony from 1 ms bins
% ---------------------------------------------------------------------------------------------

function testRandomLayoutsMatchLegacy(testCase)
    % Random wells, electrode layouts, spike patterns and burst parameters. Inputs for the
    % legacy function are cleaned of the times the current function drops. These recordings
    % are short and sparse (tens of spikes per electrode), where the synchrony index moves by
    % several hundredths with sub-millisecond timing changes under either definition, so it
    % only gets a loose check here (seeds 1-40 differ by up to 0.07; over seeds 1-600 one
    % reached 0.11). testSampleFileUnchanged checks it tightly on real data
    maxSynchronyDifference = 0;
    for seed = 1:40
        [raster_raw, params] = makeRandomRaster(seed);
        [newOutput, ~] = runCurrent(testCase, raster_raw, params);
        legacyOutput = runLegacy(testCase, cleanForLegacy(raster_raw, params.maxTime), params);
        maxSynchronyDifference = max(maxSynchronyDifference, ...
            verifyOutputsEqual(testCase, newOutput, legacyOutput, sprintf('seed %d', seed), 0.1));
        rmdir(testCase.TestData.outputRoot, 's');
        mkdir(testCase.TestData.outputRoot);
    end
    fprintf('Largest synchrony difference over random layouts: %.4f\n', maxSynchronyDifference);
end

function testBurstingElectrodeIncludesLastSample(testCase)
    % The legacy loop scanned burstStart:burstStart + duration*fs, and the float end can fall
    % just short of the last sample. That only happens for bursts that start early relative
    % to their own length (start below about duration/6, e.g. a 1 s burst in the first
    % 0.2 s); otherwise adding burstStart rounds the end back to an integer. An electrode
    % that fires only on such a last sample was not counted as bursting; the current
    % function counts it
    params = defaultParams();
    params.maxTime = 5;
    burstStart = 3;
    duration = 500;
    while numel(burstStart:burstStart + (duration / params.fs) * params.fs) == duration + 1
        duration = duration + 1;
        testCase.assertLessThan(duration, 5000, 'No duration where the legacy range drops the last sample');
    end
    leadSamples = burstStart + round(linspace(0, duration - 2, 30));
    raster_raw = cell(3, 2, 1);
    raster_raw(:, 1, 1) = {leadSamples / params.fs; (leadSamples + 1) / params.fs; ...
        (burstStart + duration) / params.fs};
    raster_raw(:, 2, 1) = {[1, 1, 1, 1, 1]; [1, 1, 2, 1, 1]; [1, 1, 3, 1, 1]};

    [newOutput, ~] = runCurrent(testCase, raster_raw, params);
    legacyOutput = runLegacy(testCase, raster_raw, params);

    isBurstingRow = strcmp(newOutput.sheet{:, 1}, 'Number of Bursting Electrodes');
    testCase.verifyEqual(newOutput.sheet{isBurstingRow, 2}, 3);
    testCase.verifyEqual(legacyOutput.sheet{isBurstingRow, 2}, 2);
    newOutput.sheet(isBurstingRow, :) = [];
    legacyOutput.sheet(isBurstingRow, :) = [];
    verifyOutputsEqual(testCase, newOutput, legacyOutput);
end

function testSingleValidSpikeInWell(testCase)
    % Two active electrodes but only one valid spike in the well, e.g. after an edit left an
    % electrode with nothing but a time at t = 0: no bursts, and no error
    params = defaultParams();
    params.minSpikesN = 2;
    raster_raw = cell(2, 2, 2);
    raster_raw(:, 2, 1) = {[1, 1, 1, 1, 1]; [1, 1, 2, 1, 1]};
    raster_raw(:, 1, 1) = {0.5; 0};
    raster_raw(:, :, 2) = makeBurstyRaster(2, 10, 7);  % a normal second well
    raster_raw{1, 2, 2}(1) = 2;
    raster_raw{2, 2, 2}(1) = 2;

    [newOutput, ~] = runCurrent(testCase, raster_raw, params);
    legacyOutput = runLegacy(testCase, cleanForLegacy(raster_raw, params.maxTime), params);

    % The second well is small (2 electrodes, 10 s), so synchrony only gets the loose check
    verifyOutputsEqual(testCase, newOutput, legacyOutput, '', 0.1);
    testCase.verifyEqual(newOutput.sheet{1, 2:3}, [0, legacyOutput.sheet{1, 3}]);
    testCase.verifyGreaterThan(newOutput.sheet{1, 3}, 0, 'The second well should still be analyzed');
    [networkBurstInfo, isBurstingElectrode] = get_network_spike_participation(12500, {5; zeros(0, 1)}, ...
        [1; 2], 0.35, 100, 2);
    testCase.verifyEqual(networkBurstInfo, zeros(6, 0));
    testCase.verifyEqual(isBurstingElectrode, [false; false]);
end

function testSynchronyFailureReturnsNaN(testCase)
    % An out-of-range sample makes the binning fail; the index must be NaN, not 0
    synchronyIndex = testCase.verifyWarning( ...
        @() calculate_multivariate_synchrony({[1; 5], [2; 100]}, 50, 12500), 'SAMS:synchronyFailed');
    testCase.verifyTrue(isnan(synchronyIndex));
end

function testSynchronyFailureReportedInExport(testCase)
    % Simulate hilbert running out of memory: the sheet gets NaN and the log says so
    stubFolder = fullfile(testCase.TestData.outputRoot, 'stub');
    mkdir(stubFolder);
    stubFile = fopen(fullfile(stubFolder, 'hilbert.m'), 'w');
    fprintf(stubFile, 'function x = hilbert(varargin)\nerror(''MATLAB:nomem'', ''Out of memory.'');\nend\n');
    fclose(stubFile);
    addpath(stubFolder);
    testCase.addTeardown(@() rmpath(stubFolder));

    params = defaultParams();
    warningState = warning('off', 'SAMS:synchronyFailed');
    testCase.addTeardown(@() warning(warningState));
    [newOutput, ~] = runCurrent(testCase, makeBurstyRaster(4, 10, 6), params);

    synchronyRow = strcmp(newOutput.sheet{:, 1}, 'Synchrony index');
    testCase.verifyTrue(isnan(newOutput.sheet{synchronyRow, 2}));
    testCase.verifySubstring(fileread(fullfile(makeOutputFolder(testCase, 'current'), 'export_log.txt')), ...
        'well 1: synchrony index could not be calculated (Out of memory.) and is reported as NaN');
end

% ---------------------------------------------------------------------------------------------
% Helpers
% ---------------------------------------------------------------------------------------------

function params = defaultParams()
    % Manual app defaults (manual_sorting_03252025.mlapp startup values)
    params.fs = 12500;
    params.maxTime = 10;
    params.threshold = 0.35;
    params.minSpikesE = 5;
    params.maxISIE = 100;
    params.minSpikesN = 50;
    params.maxISIN = 100;
end

function raster_raw = makeBurstyRaster(numElectrodes, maxTime, seed)
    % One well (A1), electrodes A1_11, A1_21, ...: ~2 Hz background plus a
    % 20-spike, 80 ms burst on every electrode every 2 s
    randomStream = RandStream('mt19937ar', 'Seed', seed);
    raster_raw = cell(numElectrodes, 2, 1);
    for electrodeNum = 1:numElectrodes
        spikeTimes = rand(randomStream, 1, round(2 * maxTime)) * maxTime;
        for burstStart = 1:2:maxTime - 1
            spikeTimes = [spikeTimes, burstStart + 0.002 * electrodeNum + rand(randomStream, 1, 20) * 0.08]; %#ok<AGROW>
        end
        raster_raw{electrodeNum, 1, 1} = sort(spikeTimes);
        raster_raw{electrodeNum, 2, 1} = [1, 1, electrodeNum, 1, 1];
    end
end

function [raster_raw, params] = makeRandomRaster(seed)
    % Up to 3 wells (one of them without data) with 0-6 electrodes each. Electrodes may be
    % emptied by an edit, hold only invalid times, share samples between units, fire past
    % the end of the recording, or carry times at or before t = 0 and NaN
    randomStream = RandStream('mt19937ar', 'Seed', seed);
    pick = @(values) values(randi(randomStream, numel(values)));
    params = defaultParams();
    params.fs = pick([12500, 20000, 12499.7]);
    params.maxTime = 5 + 15 * rand(randomStream);
    params.threshold = pick([0, 0.2, 0.35, 0.5, 0.9]);
    params.minSpikesN = pick([2, 10, 50]);
    params.maxISIN = pick([0.05, 5, 50, 100, 500]);

    numRows = 6;
    numWells = randi(randomStream, 3);
    raster_raw = cell(numRows, 2, numWells + 1);  % last well has no data
    for wellIndex = 1:numWells
        numElectrodes = randi(randomStream, [0, numRows]);
        burstTimes = sort(rand(randomStream, 1, randi(randomStream, [0, 8])) * params.maxTime);
        for electrodeNum = 1:numElectrodes
            raster_raw{electrodeNum, 2, wellIndex} = [wellIndex, 2, electrodeNum, 1, 1];
            background = rand(randomStream, 1, randi(randomStream, [1, 60])) * params.maxTime;
            bursts = [];
            for burstTime = burstTimes
                if rand(randomStream) < 0.8
                    bursts = [bursts, burstTime + rand(randomStream, 1, randi(randomStream, 25)) * 0.1]; %#ok<AGROW>
                end
            end
            spikeTimes = [background, bursts];
            % A second unit on the same sample, and one within half a sample
            spikeTimes = [spikeTimes, spikeTimes(1), spikeTimes(end) + 0.3 / params.fs];
            if rand(randomStream) < 0.3
                spikeTimes = [spikeTimes, params.maxTime + 0.01, 0.2 / params.fs];  %#ok<AGROW>
            end
            if rand(randomStream) < 0.3
                spikeTimes = [spikeTimes, -0.0008, 0, NaN];  %#ok<AGROW>
            end
            switch randi(randomStream, 10)
                case 1
                    spikeTimes = [];  % all units removed in the manual app
                case 2
                    spikeTimes = [-0.0004, NaN];  % only invalid times
            end
            raster_raw{electrodeNum, 1, wellIndex} = sort(spikeTimes);
        end
    end
end

function raster_raw = cleanForLegacy(raster_raw, maxTime)
    % Remove times at or before t = 0 and non-finite times, which the current function drops
    % (the legacy function crashed on most of them and put t = 0 on sample 1). An electrode
    % left without spikes gets one past the end of the recording, so it stays active with no
    % spikes, as it does in the current function
    for cellIndex = 1:numel(raster_raw(:, 1, :))
        [electrodeNum, ~, wellIndex] = ind2sub(size(raster_raw(:, 1, :)), cellIndex);
        spikeTimes = raster_raw{electrodeNum, 1, wellIndex};
        if isempty(spikeTimes)
            continue;
        end
        spikeTimes = spikeTimes(isfinite(spikeTimes) & spikeTimes > 0);
        if isempty(spikeTimes)
            spikeTimes = maxTime + 1;
        end
        raster_raw{electrodeNum, 1, wellIndex} = spikeTimes;
    end
end

function outputFolder = makeOutputFolder(testCase, name)
    outputFolder = fullfile(testCase.TestData.outputRoot, name);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
end

function [output, consoleText] = runCurrent(testCase, raster_raw, params, sorting_results)
    if nargin < 4
        sorting_results = {};
    end
    outputFolder = makeOutputFolder(testCase, 'current');
    consoleText = evalc(['get_network_burst_info(raster_raw, params.maxTime, params.fs, params.threshold, ', ...
        'params.minSpikesE, params.maxISIE, params.minSpikesN, params.maxISIN, outputFolder, sorting_results)']);
    output = readOutputs(outputFolder);
end

function output = runLegacy(testCase, raster_raw, params, sorting_results)
    if nargin < 4
        sorting_results = {};
    end
    outputFolder = makeOutputFolder(testCase, 'legacy');
    evalc(['get_network_burst_info_legacy(raster_raw, params.maxTime, params.fs, params.threshold, ', ...
        'params.minSpikesE, params.maxISIE, params.minSpikesN, params.maxISIN, outputFolder, sorting_results)']);
    output = readOutputs(outputFolder);
end

function [T_unit, T_electrode] = runUnitAnalysis(sorting_results, allData, maxTime, outputFolder)
    % Manual app defaults: cutoff 0.1 Hz, burst max ISI 100 ms, 5 spikes, refractory 1.5 ms
    evalc(['[T_unit, T_electrode] = get_individual_unit_analysis(sorting_results, allData, maxTime, ', ...
        '0.1, 100, 5, 0.0015, outputFolder);']);
end

function output = readOutputs(outputFolder)
    output.sheet = readtable(fullfile(outputFolder, 'spike_sorting.xlsx'), 'Sheet', 'network burst', ...
        'VariableNamingRule', 'preserve');
    loaded = load(fullfile(outputFolder, 'burst_info_all.mat'), 'burst_info_all');
    output.burst_info_all = loaded.burst_info_all;
end

function synchronyDifference = verifyOutputsEqual(testCase, actual, expected, label, synchronyTolerance)
    % Everything must match exactly except the synchrony index, which must agree within
    % synchronyTolerance (default 0.03)
    if nargin < 4
        label = '';
    end
    if nargin < 5
        synchronyTolerance = 0.03;
    end
    isSynchronyRow = strcmp(expected.sheet{:, 1}, 'Synchrony index');
    testCase.verifyEqual(actual.sheet(~isSynchronyRow, :), expected.sheet(~isSynchronyRow, :), ...
        sprintf("'network burst' sheet differs %s", label));
    actualSynchrony = actual.sheet{isSynchronyRow, 2:end};
    expectedSynchrony = expected.sheet{isSynchronyRow, 2:end};
    testCase.verifyTrue(all(isfinite(actualSynchrony)), sprintf('synchrony index not finite %s', label));
    testCase.verifyEqual(actualSynchrony, expectedSynchrony, 'AbsTol', synchronyTolerance, ...
        sprintf('synchrony index differs %s', label));
    synchronyDifference = max([abs(actualSynchrony - expectedSynchrony), 0]);
    testCase.verifyEqual(actual.burst_info_all, expected.burst_info_all, ...
        sprintf('burst_info_all differs %s', label));
end
