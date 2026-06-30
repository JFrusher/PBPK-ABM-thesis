function dosingRegimen = readDosingRegimen(inputFile)
% Read and parse the dosing regimen CSV file
% Same format as original simple PK model

    % Check if file exists
    if ~isfile(inputFile)
        errorMsg = sprintf('❌ DOSING FILE NOT FOUND\n\nFile: %s\n\n', inputFile);
        errorMsg = [errorMsg 'Troubleshooting steps:\n'];
        errorMsg = [errorMsg '1. Verify the file path is correct\n'];
        errorMsg = [errorMsg '2. Check file exists in the specified directory\n'];
        errorMsg = [errorMsg '3. Ensure file has .csv extension\n\n'];
        errorMsg = [errorMsg 'Example dosing file format:\n'];
        errorMsg = [errorMsg 'start_time_min,end_time_min,dosing_type,dose_amount\n'];
        errorMsg = [errorMsg '0,60,Bolus,1600\n\n'];
        errorMsg = [errorMsg 'Current working directory: ' pwd];
        error(errorMsg);
    end

    % Try to read the file with error handling
    try
        opts = detectImportOptions(inputFile);
        dosingRegimen = readtable(inputFile, opts);
    catch ME
        errorMsg = sprintf('❌ FAILED TO READ DOSING FILE\n\nFile: %s\n\n', inputFile);
        errorMsg = [errorMsg 'Error details: ' ME.message '\n\n'];
        errorMsg = [errorMsg 'Possible causes:\n'];
        errorMsg = [errorMsg '1. CSV file is corrupted or in wrong format\n'];
        errorMsg = [errorMsg '2. File uses incorrect delimiter (should be comma)\n'];
        errorMsg = [errorMsg '3. File contains special characters causing parse errors\n\n'];
        errorMsg = [errorMsg 'Required columns: start_time_min, end_time_min, dosing_type\n'];
        errorMsg = [errorMsg 'Optional columns: dose_amount, infusion_rate\n\n'];
        errorMsg = [errorMsg 'Original error: ' ME.message];
        error(errorMsg);
    end

    % Validate required columns (support aliases for dosing type)
    typeCandidates = {'dosing_type','dose_type','doseType','type'};
    foundType = '';
    for i = 1:length(typeCandidates)
        if any(strcmp(dosingRegimen.Properties.VariableNames, typeCandidates{i}))
            foundType = typeCandidates{i};
            break;
        end
    end

    requiredCols = {'start_time_min', 'end_time_min'};
    missingCols = {};
    for i = 1:length(requiredCols)
        if ~any(strcmp(dosingRegimen.Properties.VariableNames, requiredCols{i}))
            missingCols{end+1} = requiredCols{i};
        end
    end
    if isempty(foundType)
        missingCols{end+1} = 'dosing_type';
    end

    if ~isempty(missingCols)
        errorMsg = sprintf('❌ MISSING REQUIRED COLUMNS IN DOSING FILE\n\nFile: %s\n\n', inputFile);
        errorMsg = [errorMsg 'Missing columns: ' strjoin(missingCols, ', ') '\n\n'];
        errorMsg = [errorMsg 'Found columns: ' strjoin(dosingRegimen.Properties.VariableNames, ', ') '\n\n'];
        errorMsg = [errorMsg 'Example correct format:\n'];
        errorMsg = [errorMsg 'start_time_min,end_time_min,dosing_type,dose_amount\n'];
        errorMsg = [errorMsg '0,60,Bolus,1600\n'];
        error(errorMsg);
    end

    % Standardize dosing type column name
    if ~strcmp(foundType, 'dosing_type')
        dosingRegimen.dosing_type = dosingRegimen.(foundType);
    end

    % Standardize dose amount aliases
    if any(strcmp(dosingRegimen.Properties.VariableNames, 'dose_mg')) && ~any(strcmp(dosingRegimen.Properties.VariableNames, 'dose_amount'))
        dosingRegimen.dose_amount = dosingRegimen.dose_mg;
    end

    % Add missing optional columns with NaN / empty string
    optionalCols = {'dose_amount', 'infusion_rate', 'mean_rate', 'amplitude', 'frequency_per_min', 'rate_mg_per_min'};
    for i = 1:length(optionalCols)
        if ~any(strcmp(dosingRegimen.Properties.VariableNames, optionalCols{i}))
            dosingRegimen.(optionalCols{i}) = nan(height(dosingRegimen), 1);
        end
    end
    if ~any(strcmp(dosingRegimen.Properties.VariableNames, 'custom_function'))
        dosingRegimen.custom_function = repmat({''}, height(dosingRegimen), 1);
    end

    % Coerce numeric columns in case the CSV was read as strings
    numericCols = {'start_time_min', 'end_time_min', 'dose_amount', 'infusion_rate', ...
                   'mean_rate', 'amplitude', 'frequency_per_min', 'rate_mg_per_min'};
    for i = 1:length(numericCols)
        colName = numericCols{i};
        if any(strcmp(dosingRegimen.Properties.VariableNames, colName))
            col = dosingRegimen.(colName);
            if iscell(col) || isstring(col) || iscategorical(col)
                col = str2double(string(col));
            end
            dosingRegimen.(colName) = col;
        end
    end

    % Normalize dosing_type to lowercase trimmed string
    dosingRegimen.dosing_type = lower(strtrim(string(dosingRegimen.dosing_type)));

    % Sort by start time for deterministic overlap behavior
    [~, order] = sort(dosingRegimen.start_time_min, 'ascend');
    dosingRegimen = dosingRegimen(order, :);

    % Validate time values and create effective dosing windows.
    n = height(dosingRegimen);
    effective_end = dosingRegimen.end_time_min;
    for i = 1:n
        s = dosingRegimen.start_time_min(i);
        e = dosingRegimen.end_time_min(i);

        if ~isfinite(s)
            error('Invalid start_time_min at row %d (non-finite value).', i);
        end
        if ~isfinite(e)
            e = s;
        end

        if e < s
            warning('Row %d has end_time_min < start_time_min. Swapping values.', i);
            tmp = s; s = e; e = tmp;
            dosingRegimen.start_time_min(i) = s;
            dosingRegimen.end_time_min(i) = e;
        end

        type_i = string(dosingRegimen.dosing_type(i));
        is_bolus = any(contains(type_i, "bolus"));
        is_nonbolus = any(contains(type_i, ["continuous","constant","infusion","sinusoidal","chrono","custom","step"]));

        % Robust handling for point-style rows (start == end):
        % - bolus: keep short pulse behavior later in calculateDosingRate
        % - non-bolus: treat as step that continues until next schedule change
        if e <= s && is_nonbolus && ~is_bolus
            nextStart = NaN;
            if i < n
                futureStarts = dosingRegimen.start_time_min((i+1):end);
                futureStarts = futureStarts(futureStarts > s);
                if ~isempty(futureStarts)
                    nextStart = min(futureStarts);
                end
            end

            if isfinite(nextStart)
                e = nextStart;
            else
                e = s + 1.0; % safe fallback hold duration (1 minute)
            end
        end

        effective_end(i) = max(e, s);
    end

    dosingRegimen.effective_end_time_min = effective_end;
    dosingRegimen.effective_duration_min = max(dosingRegimen.effective_end_time_min - dosingRegimen.start_time_min, 0);

end
