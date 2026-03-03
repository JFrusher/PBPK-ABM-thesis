function diagnose_and_fix_corruption()
    % Analyzes failure patterns and attempts recovery
    
    sim_base_dir = 'simulation_results';
    sim_folders = dir(fullfile(sim_base_dir, 'sim_*'));
    sim_folders = sim_folders([sim_folders.isdir]);
    
    fprintf('=== CORRUPTION PATTERN ANALYSIS ===\n\n');
    
    for s = 1:length(sim_folders)
        sim_path = fullfile(sim_base_dir, sim_folders(s).name);
        files = dir(fullfile(sim_path, '*_cells.mat'));
        
        % Sort files by number
        [~, idx] = sort(str2double(regexp({files.name}, '\d+', 'match', 'once')));
        files = files(idx);
        
        ok_count = 0;
        failed_start = -1;
        file_sizes = [];
        ok_file_sizes = [];
        
        for f = 1:length(files)
            filepath = fullfile(sim_path, files(f).name);
            file_size = files(f).bytes;
            file_sizes(f) = file_size;
            
            try
                load(filepath, 'cells');
                ok_count = ok_count + 1;
                ok_file_sizes(f) = file_size;
            catch
                if failed_start == -1
                    failed_start = f;
                end
            end
        end
        
        if ok_count < length(files)
            fprintf('%s:\n', sim_folders(s).name);
            fprintf('  OK files: %d/%d\n', ok_count, length(files));
            fprintf('  Failure starts at file: %d (output%08d_cells.mat)\n', ...
                failed_start, failed_start - 3);
            
            % Calculate file sizes
            if ~isempty(ok_file_sizes)
                ok_sizes = ok_file_sizes(ok_file_sizes > 0);
                mean_ok_size = mean(ok_sizes);
                max_ok_size = max(ok_sizes);
                
                fprintf('  Mean OK file size: %.1f KB\n', mean_ok_size / 1024);
                fprintf('  Max OK file size: %.1f KB\n', max_ok_size / 1024);
                
                % Estimate when corruption started
                last_ok = ok_count;
                first_fail = failed_start;
                fprintf('  Last OK file: %d, First failed: %d\n', last_ok, first_fail);
                fprintf('  Corruption threshold crossed ~%.1f%% through simulation\n', ...
                    (last_ok / length(files)) * 100);
            end
            fprintf('\n');
        end
    end
    
    % Print recommendations
    fprintf('\n=== RECOMMENDATIONS ===\n\n');
    fprintf('ROOT CAUSES (in order of likelihood):\n');
    fprintf('1. MAT file writer has a 2GB limit (max array size)\n');
    fprintf('   Solution: Reduce cell array size in PhysiCell config\n\n');
    
    fprintf('2. Memory allocation fails partway through simulation\n');
    fprintf('   Solution: Reduce MAX_CELLS_PER_VOXEL or domain size\n\n');
    
    fprintf('3. Variable name length buffer overflow\n');
    fprintf('   Solution: Update PhysiCell to latest version\n\n');
    
    fprintf('4. HDF5 library version mismatch\n');
    fprintf('   Solution: Recompile with matching HDF5 version\n\n');
end
