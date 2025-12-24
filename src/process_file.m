function process_file(currentFile, file_folder, parentFolder, ...
    params, file_idx, total_files, selectedWells)
% PROCESS_FILE_WITH_SELECTED_WELLS_LOG - Processes selected wells with comprehensive logging
%
% INPUTS:
%   currentFile - File information from dir()
%   file_folder - Path to the folder containing the file
%   parentFolder - Path to the output folder
%   params - Parameter structure
%   file_idx - Current file index
%   total_files - Total number of files
%   selectedWells - Cell array of well indices to process

% Start log file
logPath = fullfile(parentFolder, 'processing_log.txt');
logFile = fopen(logPath, 'a');
% Ensure log file is closed even if an error occurs
logFileCleanup = onCleanup(@() fclose(logFile));
fprintf(logFile, 'Processing file: %s\n', currentFile.name);

% Extract file information
fileName = currentFile.name;
baseFileName = fileName(1:end-4); % Remove .spk extension

% Create output folder for this file
outputFolder = fullfile(parentFolder, baseFileName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end


% Log processing
fprintf(logFile,'Processing file %d/%d: %s\n', file_idx, total_files, fileName);

% Load spike data
filePath = fullfile(file_folder, fileName);
try
    allData = AxisFile(filePath).SpikeData.LoadData;
    SPK_dataset = AxisFile(filePath).SpikeData;
    fr = SPK_dataset.SamplingFrequency; % Hz
    fprintf(logFile,'Sampling frequency: %d Hz\n', fr);

    % Get data dimensions
    [nwr, nwc, nec, ner] = size(allData);
    fprintf(logFile,'Data dimensions: [%d, %d, %d, %d]\n', nwr, nwc, nec, ner);
catch ME
    fprintf(logFile,'Error loading file: %s\nMessage: %s\n', filePath, ME.message);
    fclose(logFile);
    return;
end

% Initialize analysis tables
[T, T_electrode, T_parameters] = initialize_tables(params);
writetable(T_parameters, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'parameter list');

% Initialize data structures
raster_raw = {};
sorting_results = {};
num_well = 1;
total_num_cell = 1;

% Initialize DTW results collection
all_DTW_tables = {};
DTW_table_count = 0;

% Get recording properties
[maxTime, maxSpikeAmp, minSpikeAmp] = get_recording_properties(allData);
fprintf(logFile,'Recording duration: %.2f seconds\n', maxTime);

% Initialize statistics
num_total_electrodes_with_spikes = 0;
num_inactive_units = 0;
num_total_detected_units = 0;
fr_inactive_units = [];
fr_active_units = [];

% Initialize PowerPoint
try
    pptx = exportToPPTX('SAMS1.pptx');
    fprintf(logFile,'PowerPoint initialized successfully.\n');
catch ME
    fprintf(logFile,'Error initializing PowerPoint: %s\n', ME.message);
    pptx = [];
end

% Initialize lists for potentially problematic electrodes
possibleMUAList = {};
flagTooMuchList = {};
flagUnitsChangedList = {};

% Process wells (filter by selectedWells if provided)
if isempty(selectedWells)
    fprintf(logFile,'Processing all wells.\n');
else
    fprintf(logFile,'Processing selected wells: %s\n', strjoin(selectedWells, ', '));
end

for i = 1:nwr
    for j = 1:nwc
        % Build well name (e.g., 'A1', 'B2', etc.)
        wellName = [char(i+'A'-1), num2str(j)];
        
        % Skip this well if selectedWells is not empty and well is not in the list
        if ~isempty(selectedWells) && ~ismember(wellName, selectedWells)
            continue;
        end
        
        num_electrode = 1;  % Reset electrode counter for each well
        for m = 1:nec
            for n = 1:ner
                % num_units_manual_check = M_checklist{check_electrode,2};
                % Update progress
                frac4 = n/ner;
                frac3 = ((m-1)+frac4)/nec;
                frac2 = ((j-1)+frac3)/nwc;
                frac1 = ((i-1)+frac2)/nwr;
                frac0 = file_idx/total_files;

                try
                    progressbar(frac0, frac1);
                catch ME
                    fprintf(logFile,'Error opening progressbar: %s', ME.message);
                end

                % Get current electrode data
                Spikes = allData{i, j, m, n}(:);

                % Skip if no spikes
                if isempty(Spikes)
                    continue;
                end

                % Get spike times and waveforms
                [Times, spikes] = Spikes.GetTimeVoltageVector;

                % Skip if total spikes are below cutoff frequency
                if size(spikes,2) < maxTime*params.cutoff_frequency
                    continue;
                end

                electrodeName = [char(i+'A'-1), num2str(j), '_', num2str(m), num2str(n)];
                fprintf('----------Processing electrode %s: %d spikes---------------\n', electrodeName, size(spikes, 2));
                fprintf(logFile, '----------Processing electrode %s: %d spikes----------\n', electrodeName, size(spikes, 2));

                % Perform PCA
                try
                    [~, score, ~] = pca(spikes');
                    % score = -score;
                catch ME
                    fprintf('Error in PCA: %s\n', ME.message);
                    fprintf(logFile, 'Error in PCA for electrode %s: %s\n', electrodeName, ME.message);
                    continue;
                end

                % Process electrode
                [electrode_results, electrode_stats, spikes_aligned] = process_electrode(logFile, spikes, Times, i, j, m, n, ...
                    params, maxTime, maxSpikeAmp, minSpikeAmp, score);
         

                % Skip if no units detected
                if isempty(electrode_results.new_idx_list)
                    continue;
                end

                % Log unit detection
                fprintf(logFile, 'Electrode %s: %d units detected\n', electrodeName, length(electrode_results.new_idx_list));

                [T, T_electrode, ...
                    raster_raw, sorting_results, num_electrode, total_num_cell, ...
                    num_total_electrodes_with_spikes, num_inactive_units, num_total_detected_units, ...
                    fr_inactive_units, fr_active_units, possibleMUAList, flagTooMuchList] = ...
                    update_results(params.flag_threshold, electrode_results, electrode_stats, ...
                    T, T_electrode, ...
                    raster_raw, ...
                    sorting_results, num_electrode, total_num_cell, num_well, i, j, m, n, ...
                    num_total_electrodes_with_spikes, num_inactive_units, num_total_detected_units, ...
                    fr_inactive_units, fr_active_units, possibleMUAList, flagTooMuchList, pptx, ...
                    maxTime, baseFileName, spikes_aligned, Times, score);

            end
        end
        % Increment well counter after processing all electrodes in this well
        num_well = num_well + 1;
    end
end

% Update progress bar with completion message for this file
h = [];  % Initialize handle for cleanup
try
    % Custom progress bar with completion message
    h = figure('Position', [400, 400, 400, 100], 'Name', 'Processing Status', ...
        'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off');
    % Ensure figure is closed even if an error occurs later
    figureCleanup = onCleanup(@() closeIfValid(h));
    uicontrol('Style', 'text', 'Position', [20, 50, 360, 30], ...
        'String', ['File ' num2str(file_idx) '/' num2str(total_files) ' completed. Writing outputs...'], ...
        'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    uicontrol('Style', 'text', 'Position', [20, 20, 360, 20], ...
        'String', 'Please wait. Do not close MATLAB.', ...
        'FontSize', 10, 'HorizontalAlignment', 'center');
    drawnow;
catch
    % If figure creation fails, just print a message
    fprintf(logFile,'\nFile %d/%d completed. Writing outputs...\n', file_idx, total_files);
    fprintf(logFile,'Please wait. Do not close MATLAB.\n');
end


% Save PowerPoint
if ~isempty(pptx)
    try
        pptx.save([outputFolder, '\', baseFileName]);
        fprintf(logFile,'PowerPoint presentation saved successfully.\n');
    catch ME
        fprintf(logFile,'Error saving PowerPoint: %s\n', ME.message);
    end
end

% Save burst info
fprintf(logFile, 'Saving burst analysis data...\n');
save([outputFolder, '\burst_info_all.mat'], 'raster_raw', 'maxTime', 'sorting_results', '-v7.3');

% Perform network burst analysis
fprintf(logFile,'Performing network burst analysis...\n');
try
    get_network_burst_info(raster_raw, maxTime, fr, params.network_participation_threshold, ...
        params.min_spikes_E, params.max_ISI_E, params.min_spikes_N, params.max_ISI_N, ...
        outputFolder, sorting_results);
    fprintf(logFile, 'Network burst analysis completed successfully.\n');
catch ME
    fprintf(logFile, 'Error in network burst analysis: %s\n', ME.message);
end

% Save results to Excel
fprintf(logFile,'Writing results to Excel...\n');
try
    writetable(T, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'individual unit');
    writetable(T_electrode, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'electrode statistics');

    fprintf(logFile, 'Main results written to Excel successfully.\n');
catch ME
    fprintf(logFile, 'Error writing main results to Excel: %s\n', ME.message);
end

% Save check lists
if isempty(possibleMUAList)
    possibleMUAList{end+1} = 'no check list';
end
if isempty(flagTooMuchList)
    flagTooMuchList{end+1} = 'no check list';
end

% Make both tables have the same length
maxLength = max(length(possibleMUAList), length(flagTooMuchList));
while length(possibleMUAList) < maxLength
    possibleMUAList{end+1} = '';
end
while length(flagTooMuchList) < maxLength
    flagTooMuchList{end+1} = '';
end

% Create tables with equal numbers of rows
try
    T_pmua = cell2table(possibleMUAList', 'VariableNames', "Possible MultiUnit");
    T_flag = cell2table(flagTooMuchList', 'VariableNames', "over-exlcuded unit");
    T_checklist = [T_pmua T_flag];
    writetable(T_checklist, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'check list (active)');

    fprintf(logFile, 'Check lists written to Excel successfully.\n');
catch ME
    fprintf(logFile, 'Error writing check lists to Excel: %s\n', ME.message);
end

% Print and log final summary
summary_text = sprintf('\n=== PROCESSING SUMMARY ===\n');
summary_text = [summary_text, sprintf('File: %s\n', fileName)];
summary_text = [summary_text, sprintf('Processing mode: All wells\n')];
summary_text = [summary_text, sprintf('Total electrodes with spikes: %d\n', num_total_electrodes_with_spikes)];
summary_text = [summary_text, sprintf('Total units detected: %d\n', num_total_detected_units)];
summary_text = [summary_text, sprintf('Active units: %d\n', num_total_detected_units - num_inactive_units)];
summary_text = [summary_text, sprintf('Inactive units: %d\n', num_inactive_units)];
if DTW_table_count > 0
    summary_text = [summary_text, sprintf('Electrodes with DTW analysis: %d\n', DTW_table_count)];
end
summary_text = [summary_text, sprintf('Possible MUA electrodes: %d\n', length(possibleMUAList) - double(strcmp(possibleMUAList{end}, 'no check list') || strcmp(possibleMUAList{end}, '')))];
summary_text = [summary_text, sprintf('Over-excluded electrodes: %d\n', length(flagTooMuchList) - double(strcmp(flagTooMuchList{end}, 'no check list') || strcmp(flagTooMuchList{end}, '')))];
summary_text = [summary_text, sprintf('Recording duration: %.2f seconds\n', maxTime)];
summary_text = [summary_text, sprintf('Output folder: %s\n', outputFolder)];
summary_text = [summary_text, sprintf('Log file: %s\n', logPath)];
summary_text = [summary_text, sprintf('=========================================================\n\n')];

fprintf('%s', summary_text);
fprintf(logFile, '%s', summary_text);

% Close the notification window now that all outputs are written
try
    if exist('h', 'var') && ishandle(h)
        fprintf(logFile,'All outputs written successfully. Closing notification...\n');
        close(h);
    end
catch ME
    fprintf(logFile,'Error closing notification window: %s\n', ME.message);
end

fprintf('File processing completed: %s\n', fileName);
fprintf(logFile,'File processing completed: %s\n', fileName);
fprintf(logFile,'========================================\n\n');

% Log file will be automatically closed by onCleanup when function exits
% (No need to manually call fclose - handled by logFileCleanup)
end

% Helper function to safely close figure handles
function closeIfValid(h)
    if ~isempty(h) && isvalid(h)
        close(h);
    end
end