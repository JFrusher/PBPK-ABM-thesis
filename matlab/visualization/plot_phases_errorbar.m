function plot_simulation_results()
    % Configuration - EASY TO CHANGE
    column_to_analyse = 8; % current_phase (change this to analyse different columns)
    interval_min = 60; % Time interval in minutes
    sim_base_dir = 'simulation_results'; % Base directory containing sim_* folders
    
    % Phase and cell type definitions
    phase_values = [4, 10, 12, 13, 17, 18, 100, 101, 102];
    phase_labels = {'G1', 'S', 'G2', 'M', 'Apoptosis', 'Necrosis', 'Dead1', 'Dead2', 'Dead3'};
    
    % Alive phases (G1, S, G2, M)
    alive_phases = [4, 10, 12, 13];
    alive_labels = {'G1', 'S', 'G2', 'M'};
    
    % Dead phases (Dead1, Dead2, Dead3)
    dead_phases = [100, 101, 102];
    dead_labels = {'Dead1', 'Dead2', 'Dead3'};
    
    % Check if base directory exists
    if ~isfolder(sim_base_dir)
        error('Base directory not found: %s', sim_base_dir);
    end
    
    % Get all sim_* folders
    sim_folders = dir(fullfile(sim_base_dir, 'sim_*'));
    sim_folders = sim_folders([sim_folders.isdir]);
    
    if isempty(sim_folders)
        error('No sim_* folders found in %s', sim_base_dir);
    end
    
    fprintf('Found %d simulation folders\n', length(sim_folders));
    
    % Store all data for each simulation
    all_sim_data = struct();
    unique_times = [];
    valid_sims = 0;
    
    % Load all data from all simulations
    for s = 1:length(sim_folders)
        sim_path = fullfile(sim_base_dir, sim_folders(s).name);
        files = dir(fullfile(sim_path, '*_cells.mat'));
        
        if isempty(files)
            fprintf('Warning: No *_cells.mat files in %s (skipping)\n', sim_folders(s).name);
            continue;
        end
        
        fprintf('Processing %s (%d files)...\n', sim_folders(s).name, length(files));
        
        sim_times = [];
        sim_data = [];
        file_count = 0;
        
        for f = 1:length(files)
            fname = files(f).name;
            filepath = fullfile(sim_path, fname);
            
            try
                % Load with error handling
                load_data = load(filepath, 'cells');
                
                if ~isfield(load_data, 'cells')
                    fprintf('  Warning: No cells field in %s (skipping)\n', fname);
                    continue;
                end
                
                cells = load_data.cells;
                
                % Validate cells array
                if ~isnumeric(cells) || isempty(cells)
                    fprintf('  Warning: Invalid cells array in %s (skipping)\n', fname);
                    continue;
                end
                
                % Check if column exists
                if column_to_analyse > size(cells, 1)
                    fprintf('  Warning: Column %d does not exist in %s (file has %d rows, skipping)\n', ...
                        column_to_analyse, fname, size(cells, 1));
                    continue;
                end
                
                % Extract column data safely
                col_data = cells(column_to_analyse, :);
                col_data = col_data(isfinite(col_data)); % Remove NaN and Inf
                
                if isempty(col_data)
                    fprintf('  Warning: No valid data in column %d of %s (skipping)\n', column_to_analyse, fname);
                    continue;
                end
                
                col_data = round(col_data);
                
                % Get time from filename
                numMatch = regexp(fname, 'output(\d+)_', 'tokens');
                if ~isempty(numMatch)
                    t = str2double(numMatch{1}{1});
                else
                    t = f - 1;
                end
                
                sim_times(end+1) = t;
                
                % Count occurrences of each phase/value
                counts = zeros(1, length(phase_values));
                for k = 1:length(phase_values)
                    counts(k) = sum(col_data == phase_values(k));
                end
                
                sim_data(end+1, :) = counts;
                file_count = file_count + 1;
                
            catch ME
                fprintf('  Error processing %s: %s (skipping)\n', fname, ME.message);
                continue;
            end
        end
        
        % Only store if we have valid data
        if file_count > 0
            sim_times = sim_times * interval_min / 60; % Convert to hours
            unique_times = union(unique_times, sim_times);
            
            valid_sims = valid_sims + 1;
            all_sim_data(valid_sims).times = sim_times;
            all_sim_data(valid_sims).data = sim_data;
            all_sim_data(valid_sims).name = sim_folders(s).name;
            fprintf('  Successfully loaded %d timepoints\n', file_count);
        else
            fprintf('  No valid data from %s\n', sim_folders(s).name);
        end
    end
    
    % Validate we have data
    if valid_sims == 0
        error('No valid simulation data found');
    end
    
    if isempty(unique_times)
        error('No valid timepoints found');
    end
    
    unique_times = sort(unique_times);
    fprintf('\nProcessing complete: %d simulations with valid data\n', valid_sims);
    fprintf('Time range: %.2f to %.2f hours\n', min(unique_times), max(unique_times));
    
    % Interpolate all simulations to common timepoints and calculate stats
    num_phases = length(phase_values);
    mean_data = zeros(length(unique_times), num_phases);
    std_data = zeros(length(unique_times), num_phases);
    
    for t_idx = 1:length(unique_times)
        t = unique_times(t_idx);
        phase_counts = [];
        
        for s = 1:valid_sims
            times = all_sim_data(s).times;
            data = all_sim_data(s).data;
            
            % Find closest timepoint
            [~, closest_idx] = min(abs(times - t));
            phase_counts(s, :) = data(closest_idx, :);
        end
        
        mean_data(t_idx, :) = mean(phase_counts, 1);
        std_data(t_idx, :) = std(phase_counts, 1);
    end
    
    % Plot 1: Alive Phases
    try
        fig1 = figure('Position', [100 100 1000 500]);
        plot_with_errorbars(unique_times, mean_data, std_data, alive_phases, alive_labels, 'Alive Phases Over Time');
    catch ME
        fprintf('Error creating alive phases plot: %s\n', ME.message);
    end
    
    % Plot 2: Dead Phases
    try
        fig2 = figure('Position', [100 650 1000 500]);
        plot_with_errorbars(unique_times, mean_data, std_data, dead_phases, dead_labels, 'Dead Phases Over Time');
    catch ME
        fprintf('Error creating dead phases plot: %s\n', ME.message);
    end
    
    % Plot 3: Cell Type
    try
        all_unique_vals = [];
        for s = 1:valid_sims
            all_unique_vals = [all_unique_vals; all_sim_data(s).data(:)];
        end
        unique_cell_types = unique(all_unique_vals);
        
        fig3 = figure('Position', [1150 100 1000 500]);
        plot_with_errorbars(unique_times, mean_data, std_data, unique_cell_types, ...
            arrayfun(@num2str, unique_cell_types, 'UniformOutput', false), ...
            sprintf('Column %d Values Over Time', column_to_analyse));
    catch ME
        fprintf('Error creating cell type plot: %s\n', ME.message);
    end
    
end

% Helper function to plot with error bars
function plot_with_errorbars(times, mean_data, std_data, phase_indices, phase_labels, title_str)
    
    % Validate inputs
    if isempty(times) || isempty(mean_data) || isempty(std_data)
        warning('Skipping plot %s: empty data', title_str);
        return;
    end
    
    if length(phase_indices) == 0
        warning('Skipping plot %s: no phases to plot', title_str);
        return;
    end
    
    % Find indices of phases to plot
    indices_to_plot = [];
    labels_to_plot = {};
    all_phases = [4, 10, 12, 13, 17, 18, 100, 101, 102];
    
    for p = 1:length(phase_indices)
        phase_val = phase_indices(p);
        
        % Handle both numeric phase values and arbitrary column values
        if ismember(phase_val, all_phases)
            idx = find(all_phases == phase_val, 1);
        else
            % For arbitrary column values, map directly
            idx = p;
        end
        
        if ~isempty(idx) && idx <= size(mean_data, 2)
            indices_to_plot = [indices_to_plot, idx];
            if iscell(phase_labels) && p <= length(phase_labels)
                labels_to_plot = [labels_to_plot, phase_labels(p)];
            else
                labels_to_plot = [labels_to_plot, {num2str(phase_val)}];
            end
        end
    end
    
    if isempty(indices_to_plot)
        warning('Skipping plot %s: no valid indices to plot', title_str);
        return;
    end
    
    % Plot each phase with error bars
    hold on;
    num_plots = length(indices_to_plot);
    colors = lines(max(num_plots, 3));
    
    for i = 1:num_plots
        try
            idx = indices_to_plot(i);
            means = mean_data(:, idx);
            stds = std_data(:, idx);
            
            % Validate data
            if any(~isfinite(means))
                fprintf('Warning: Non-finite values in means for %s\n', labels_to_plot{i});
                means(~isfinite(means)) = 0;
            end
            if any(~isfinite(stds))
                stds(~isfinite(stds)) = 0;
            end
            
            errorbar(times, means, stds, '-', ...
                'Color', colors(i, :), ...
                'LineWidth', 2, ...
                'MarkerSize', 0.05, ...
                'CapSize', 3, ...
                'DisplayName', labels_to_plot{i});
        catch ME
            fprintf('Error plotting %s: %s\n', labels_to_plot{i}, ME.message);
            continue;
        end
    end
    
    xlabel('Time (hours)', 'FontSize', 12);
    ylabel('Cell Count', 'FontSize', 12);
    title(title_str, 'FontSize', 14);
    
    if ~isempty(labels_to_plot)
        legend(labels_to_plot, 'Location', 'Best', 'FontSize', 11);
    end
    
    grid on;
    hold off;
    
end
 
