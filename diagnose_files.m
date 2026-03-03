function diagnose_file()
    
    sim_path = 'simulation_results/sim_20251129_162959_05'; % Change to your sim folder
    
    % Get all *_cells.mat files
    files = dir(fullfile(sim_path, '*_cells.mat'));
    
    fprintf('Total files found: %d\n\n', length(files));
    
    for f = 1:length(files)
        fname = files(f).name;
        filepath = fullfile(sim_path, fname);
        file_size = files(f).bytes;
        
        % Extract file number
        numMatch = regexp(fname, 'output(\d+)_', 'tokens');
        if ~isempty(numMatch)
            file_num = str2double(numMatch{1}{1});
        else
            file_num = NaN;
        end
        
        fprintf('File %s (%.2f KB): ', fname, file_size/1024);
        
        try
            load_data = load(filepath);
            
            if isfield(load_data, 'cells')
                cells = load_data.cells;
                fprintf('OK - cells array [%d x %d]\n', size(cells, 1), size(cells, 2));
            else
                fprintf('ERROR - No cells field. Fields: %s\n', strjoin(fieldnames(load_data), ', '));
            end
            
        catch ME
            fprintf('FAILED - %s\n', ME.message);
        end
    end
    
end
