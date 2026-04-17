function diagnose_mat_files()
    % COMPREHENSIVE .MAT FILE DIAGNOSTIC SCRIPT
    % Scans all sim_* folders in simulation_results/
    % Identifies where errors start in each sim folder
    % Tests a sample of bad files and generates detailed report
    
    clear; clc;
    
    % Configuration
    base_dir = 'simulation_results';
    report_file = 'mat_diagnostic_report.txt';
    sample_size = 5; % Test up to 5 bad files per sim
    
    % Open report file
    fid = fopen(report_file, 'w');
    fprintf(fid, '=== COMPREHENSIVE .MAT FILE DIAGNOSTIC REPORT ===\n');
    fprintf(fid, 'Generated: %s\n', datetime('now'));
    fprintf(fid, 'Base directory: %s\n\n', base_dir);
    
    % Find all sim_* folders
    sim_folders = dir(fullfile(base_dir, 'sim_*'));
    sim_folders = sim_folders([sim_folders.isdir]);
    
    if isempty(sim_folders)
        fprintf(fid, 'ERROR: No sim_* folders found in %s\n', base_dir);
        fclose(fid);
        return;
    end
    
    fprintf('Found %d simulation folders\n', length(sim_folders));
    fprintf(fid, 'Found %d simulation folders\n\n', length(sim_folders));
    
    % Process each sim folder
    all_stats = [];
    
    for sim_idx = 1:length(sim_folders)
        sim_name = sim_folders(sim_idx).name;
        sim_path = fullfile(base_dir, sim_name);
        
        fprintf('\n=== Processing %s ===\n', sim_name);
        fprintf(fid, '\n==============================================\n');
        fprintf(fid, 'SIMULATION: %s\n', sim_name);
        fprintf(fid, '==============================================\n\n');
        
        % Get all *_cells.mat files
        files = dir(fullfile(sim_path, '*_cells.mat'));
        
        if isempty(files)
            fprintf('  No .mat files found\n');
            fprintf(fid, 'No .mat files found in this folder\n\n');
            continue;
        end
        
        fprintf('  Total .mat files: %d\n', length(files));
        fprintf(fid, 'Total .mat files: %d\n\n', length(files));
        
        % Classify files as OK or FAILED
        file_status = {};
        error_start_idx = [];
        first_error = true;
        
        for f = 1:length(files)
            fname = files(f).name;
            filepath = fullfile(sim_path, fname);
            
            % Try to load and check
            [is_ok, error_msg, data_info] = check_file(filepath);
            
            file_status{f, 1} = fname;
            file_status{f, 2} = is_ok;
            file_status{f, 3} = error_msg;
            file_status{f, 4} = data_info;
            
            % Track where errors start
            if ~is_ok && first_error
                error_start_idx = f;
                first_error = false;
            end
        end
        
        % Count OK vs FAILED
        num_ok = sum([file_status{:, 2}]);
        num_failed = length(files) - num_ok;
        
        fprintf('  Status: %d OK, %d FAILED\n', num_ok, num_failed);
        fprintf(fid, 'STATUS SUMMARY:\n');
        fprintf(fid, '  OK files: %d\n', num_ok);
        fprintf(fid, '  FAILED files: %d\n', num_failed);
        
        % Find contiguous regions
        if ~isempty(error_start_idx)
            fprintf(fid, '  Error pattern starts at: File %d (%s)\n', ...
                error_start_idx, file_status{error_start_idx, 1});
            fprintf('  Errors start at file %d\n', error_start_idx);
            
            % Check if all files from error_start_idx onwards are bad
            errors_from_here = all(~[file_status{error_start_idx:end, 2}]);
            if errors_from_here
                fprintf(fid, '  Pattern: ALL files from #%d onwards are corrupted\n', error_start_idx);
                fprintf('  All files corrupted from #%d onwards\n', error_start_idx);
            else
                fprintf(fid, '  Pattern: Sporadic errors (not all files from start point)\n');
                fprintf('  Pattern: Mix of OK and FAILED files\n');
            end
        else
            fprintf(fid, '  All files loaded successfully!\n');
        end
        
        fprintf(fid, '\nDETAILED FILE ANALYSIS:\n');
        fprintf(fid, '------------------------\n');
        
        % Sample bad files for detailed analysis
        if num_failed > 0
            bad_indices = find(~[file_status{:, 2}]);
            sample_indices = bad_indices(1:min(sample_size, length(bad_indices)));
            
            fprintf('\n  Testing %d bad files in detail...\n', length(sample_indices));
            fprintf(fid, '\nDETAILED ANALYSIS OF %d SAMPLE BAD FILES:\n', length(sample_indices));
            fprintf(fid, '\n');
            
            for i = 1:length(sample_indices)
                idx = sample_indices(i);
                fname = file_status{idx, 1};
                error_msg = file_status{idx, 3};
                data_info = file_status{idx, 4};
                
                fprintf(fid, '--- File %d: %s ---\n', idx, fname);
                fprintf(fid, 'Error: %s\n', error_msg);
                
                if ~isempty(data_info)
                    fprintf(fid, 'Data Info:\n');
                    fprintf(fid, '%s\n', data_info);
                end
                fprintf(fid, '\n');
            end
        end
        
        % Test one OK file for reference
        ok_indices = find([file_status{:, 2}]);
        if ~isempty(ok_indices)
            idx = ok_indices(1);
            fname = file_status{idx, 1};
            data_info = file_status{idx, 4};
            
            fprintf(fid, '\nREFERENCE (GOOD FILE):\n');
            fprintf(fid, '--- File %d: %s ---\n', idx, fname);
            fprintf(fid, 'Status: OK\n');
            if ~isempty(data_info)
                fprintf(fid, 'Data Info:\n');
                fprintf(fid, '%s\n', data_info);
            end
            fprintf(fid, '\n');
        end
        
        % Store stats
        all_stats = [all_stats; struct('sim', sim_name, 'total', length(files), ...
                                       'ok', num_ok, 'failed', num_failed, ...
                                       'error_starts_at', error_start_idx)];
    end
    
    % Summary section
    fprintf(fid, '\n\n==============================================\n');
    fprintf(fid, 'OVERALL SUMMARY\n');
    fprintf(fid, '==============================================\n\n');
    
    total_files = sum([all_stats.total]);
    total_ok = sum([all_stats.ok]);
    total_failed = sum([all_stats.failed]);
    
    fprintf(fid, 'Across all simulations:\n');
    fprintf(fid, '  Total files scanned: %d\n', total_files);
    fprintf(fid, '  OK files: %d (%.1f%%)\n', total_ok, 100*total_ok/total_files);
    fprintf(fid, '  FAILED files: %d (%.1f%%)\n', total_failed, 100*total_failed/total_files);
    
    fprintf(fid, '\nPER-SIMULATION BREAKDOWN:\n');
    for i = 1:length(all_stats)
        s = all_stats(i);
        fprintf(fid, '  %s: %d OK, %d FAILED', s.sim, s.ok, s.failed);
        if ~isempty(s.error_starts_at)
            fprintf(fid, ' [errors from file %d]', s.error_starts_at);
        end
        fprintf(fid, '\n');
    end
    
    % Recommendations
    fprintf(fid, '\n\nRECOMMENDATIONS:\n');
    fprintf(fid, '1. Common causes of corrupted files:\n');
    fprintf(fid, '   - Simulation crashed or stopped prematurely\n');
    fprintf(fid, '   - Disk full or write permission error\n');
    fprintf(fid, '   - Incomplete file write (partial data)\n');
    fprintf(fid, '   - PhysiCell version mismatch\n');
    fprintf(fid, '\n2. Actions to take:\n');
    fprintf(fid, '   - If errors are at the END of a sim: Files are incomplete, safe to exclude\n');
    fprintf(fid, '   - If errors are THROUGHOUT a sim: Rerun entire simulation\n');
    fprintf(fid, '   - If only a few sims affected: Rerun just those sims\n');
    fprintf(fid, '   - Keep backup of working files (e.g., sim_6) before rerunning\n');
    
    fclose(fid);
    
    fprintf('\n=== DIAGNOSTIC COMPLETE ===\n');
    fprintf('Report written to: %s\n', report_file);
    fprintf('Open it to see detailed analysis.\n');
end

% ===== HELPER FUNCTION =====
function [is_ok, error_msg, data_info] = check_file(filepath)
    % Check if file can be loaded and provide detailed info
    
    is_ok = false;
    error_msg = '';
    data_info = '';
    
    % Check file existence
    if ~isfile(filepath)
        error_msg = 'File not found';
        return;
    end
    
    % Get file size
    file_info = dir(filepath);
    file_size_kb = file_info.bytes / 1024;
    
    % Try to load
    try
        load_data = load(filepath);
        
        % Check if cells field exists
        if ~isfield(load_data, 'cells')
            is_ok = false;
            error_msg = 'No cells field in file';
            data_info = sprintf('Fields present: %s', strjoin(fieldnames(load_data), ', '));
            return;
        end
        
        cells = load_data.cells;
        
        % Check data integrity
        if isempty(cells)
            is_ok = false;
            error_msg = 'Cells array is empty';
            return;
        end
        
        % Check for NaN or Inf
        num_nan = sum(isnan(cells(:)));
        num_inf = sum(isinf(cells(:)));
        
        if num_nan > 0 || num_inf > 0
            is_ok = false;
            error_msg = sprintf('Data corruption detected: %d NaN, %d Inf values', num_nan, num_inf);
            data_info = sprintf('Array shape: %s, File size: %.1f KB', ...
                mat2str(size(cells)), file_size_kb);
            return;
        end
        
        % File is OK
        is_ok = true;
        error_msg = 'OK';
        data_info = sprintf('Array shape: %s, File size: %.1f KB, Values range: [%.2e, %.2e]', ...
            mat2str(size(cells)), file_size_kb, min(cells(:)), max(cells(:)));
        
    catch ME
        is_ok = false;
        error_msg = ME.message;
        data_info = sprintf('File size: %.1f KB', file_size_kb);
    end
end
