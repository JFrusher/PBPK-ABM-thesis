function diagnose_file_loading()
    % DIAGNOSTIC SCRIPT FOR FILE LOADING AND ORDERING ISSUES
    % Mirrors the exact file pickup logic from plot_simulation_results.m
    % Helps identify spike issues at start of plots

    % Configuration - MUST MATCH plot_simulation_results.m
    sim_base_dir = 'simulation_results\sim_20251211_205200';
    column_to_analyse = 8; % current_phase
    interval_min = 60;

    fprintf('\n====================================================================\n');
    fprintf('FILE LOADING DIAGNOSTIC REPORT\n');
    fprintf('====================================================================\n\n');

    fprintf('Base Directory: %s\n', sim_base_dir);
    fprintf('Column to Analyse: %d\n', column_to_analyse);
    fprintf('Interval: %d minutes\n\n', interval_min);

    % Check if directory exists
    if ~isfolder(sim_base_dir)
        error('Base directory not found: %s', sim_base_dir);
    end

    % Get all sim_* folders
    sim_folders = dir(fullfile(sim_base_dir, 'sim_*'));
    sim_folders = sim_folders([sim_folders.isdir]);

    if isempty(sim_folders)
        error('No sim_* folders found');
    end

    fprintf('Found %d simulation folders\n\n', length(sim_folders));

    % Diagnostic data collection
    for s = 1:min(length(sim_folders), 3)  % Check first 3 simulations for diagnosis
        sim_path = fullfile(sim_base_dir, sim_folders(s).name);
        files = dir(fullfile(sim_path, '*_cells.mat'));

        fprintf('\n--- Simulation: %s (%d files) ---\n', sim_folders(s).name, length(files));
        fprintf('Files as discovered by dir():\n');

        % Show files in order discovered (may not be sorted!)
        for f = 1:min(length(files), 10)  % Show first 10 files
            fprintf('  [%d] %s\n', f, files(f).name);
        end

        if length(files) > 10
            fprintf('  ... and %d more files\n', length(files) - 10);
        end

        fprintf('\nPARSING RESULTS:\n');
        fprintf('File Index | Filename | Extracted Time | Converted Time (hours) | Cell Data Sample\n');
        fprintf('%-10s | %-30s | %-14s | %-22s | %s\n', ...
            '---', '---', '---', '---', '---');

        sim_times = [];

        for f = 1:length(files)
            fname = files(f).name;
            filepath = fullfile(sim_path, fname);

            try
                % Load file
                load_data = load(filepath, 'cells');
                cells = load_data.cells;

                % Extract column data
                if column_to_analyse <= size(cells, 1)
                    col_data = cells(column_to_analyse, :);
                    col_data = col_data(isfinite(col_data));
                    col_data = round(col_data);

                    % Extract by regex
                    numMatch = regexp(fname, 'output(\d+)_', 'tokens');
                    if ~isempty(numMatch)
                        t_regex = str2double(numMatch{1}{1});
                    else
                        t_regex = NaN;
                    end

                    % Convert to hours
                    t_hours_regex = t_regex * interval_min / 60;

                    % Show data
                    sample_val = col_data(min(10, length(col_data)));  % 10th value or last

                    fprintf('%-10d | %-30s | %6.0f | %6.2f hours | %d\n', ...
                        f, fname, t_regex, t_hours_regex, sample_val);

                    sim_times = [sim_times; t_regex, t_hours_regex, f];

                    if f <= 3  % Show details for first 3 files
                        fprintf('      >> File %d details: %d cells with phase data\n', f, length(col_data));
                        fprintf('      >> Min: %d, Max: %d, Mean: %.1f\n', ...
                            min(col_data), max(col_data), mean(col_data));
                    end
                end

            catch ME
                fprintf('%-10d | %-30s | ERROR: %s\n', f, fname, ME.message);
            end
        end

        % ANALYSIS: Check for ordering issues
        fprintf('\n*** ORDERING ANALYSIS ***\n');

        if ~isempty(sim_times)
            % Check if times are monotonically increasing
            time_diffs = diff(sim_times(:, 1));

            if all(time_diffs >= 0)
                fprintf('OK: Times are MONOTONIC (properly ordered)\n');
            else
                fprintf('WARNING: Times not monotonic - files may be out of order!\n');
                bad_indices = find(time_diffs < 0);
                if ~isempty(bad_indices)
                    fprintf('  First disorder at transition %d -> %d\n', bad_indices(1), bad_indices(1)+1);
                    fprintf('  Time goes from %.0f to %.0f minutes\n', ...
                        sim_times(bad_indices(1), 1), sim_times(bad_indices(1)+1, 1));
                end
            end

            % Check for gaps
            fprintf('\nTime point gaps analysis:\n');
            for idx = 1:min(5, size(sim_times, 1))
                if idx < size(sim_times, 1)
                    gap = sim_times(idx+1, 1) - sim_times(idx, 1);
                else
                    gap = 0;
                end
                fprintf('  File %d (t=%d min): gap = %.0f min to next\n', ...
                    idx, sim_times(idx, 1), gap);
            end
        end

        fprintf('\n');
    end

    % RECOMMENDATIONS
    fprintf('\n====================================================================\n');
    fprintf('DIAGNOSTIC RECOMMENDATIONS\n');
    fprintf('====================================================================\n\n');

    fprintf('COMMON CAUSES OF START-OF-PLOT SPIKES:\n');
    fprintf('1. FILES NOT SORTED - dir() returns files in filesystem order\n');
    fprintf('   FIX: Sort files by extracted timepoint\n\n');

    fprintf('2. FIRST TIMEPOINT MISMATCH\n');
    fprintf('   Check: Does first file start at output00000000 or output00000001?\n\n');

    fprintf('3. TIME EXTRACTION MISMATCH\n');
    fprintf('   Verify regex extraction matches your actual filenames\n\n');

    fprintf('4. INTERPOLATION AT BOUNDARIES\n');
    fprintf('   First timepoint may interpolate poorly with sparse data\n\n');

    fprintf('NEXT STEPS:\n');
    fprintf('a) Note the first 3-5 filenames and times from output above\n');
    fprintf('b) Check if times are monotonic (should always increase)\n');
    fprintf('c) If NOT monotonic, files need sorting before processing\n');
    fprintf('d) If times jump unexpectedly, check time conversion logic\n');
    fprintf('\n');
end