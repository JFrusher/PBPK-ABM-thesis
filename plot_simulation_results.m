function plot_simulation_results()

% Configuration - EASY TO CHANGE
column_to_analyse = 8; % current_phase (change this to analyse different columns)
interval_min = 60; % Time interval in minutes
sim_base_dir = 'simulation_results/sim_20251213_220411'; % Base directory containing sim_* folders

% Phase and cell type definitions
phase_values = [4, 10, 12, 13, 17, 18, 100, 101, 102];
phase_labels = {'G1', 'S', 'G2', 'M', 'Apoptosis', 'Necrosis', 'Dead1', 'Dead2', 'Dead3'};

% Alive phases (G1, S, G2, M)
alive_phases = [4, 10, 12, 13];
alive_labels = {'G1', 'S', 'G2', 'M'};

% Dead phases (Apoptosis, Necrosis, Dead1, Dead2, Dead3)
dead_phases = [17, 18, 100, 101, 102];
dead_labels = {'Apoptosis', 'Necrosis', 'Dead1', 'Dead2', 'Dead3'};

% Cell types 
column_celltype = 6; % Label index 5 = Row 6 in MATLAB
celltype_values = [0, 1, 2, 3, 4, 5];
celltype_labels = {'Cancer_Stem', 'Cancer_Prolif', 'Cancer_Diff', 'CAF', 'CD8_T', 'M2_Macro'};

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
sim_file_counts = [];

% Load all data from all simulations
for s = 1:length(sim_folders)
    sim_path = fullfile(sim_base_dir, sim_folders(s).name);
    files = dir(fullfile(sim_path, 'output*_cells.mat'));
    
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
                fprintf(' Warning: No cells field in %s (skipping)\n', fname);
                continue;
            end
            
            cells = load_data.cells;
            
            % Validate cells array
            if ~isnumeric(cells) || isempty(cells)
                fprintf(' Warning: Invalid cells array in %s (skipping)\n', fname);
                continue;
            end
            
            % Check if column exists
            if column_to_analyse > size(cells, 1)
                fprintf(' Warning: Column %d does not exist in %s (file has %d rows, skipping)\n', ...
                    column_to_analyse, fname, size(cells, 1));
                continue;
            end
            
            % Extract column data safely
            col_data = cells(column_to_analyse, :);
            col_data = col_data(isfinite(col_data)); % Remove NaN and Inf
            
            if isempty(col_data)
                fprintf(' Warning: No valid data in column %d of %s (skipping)\n', column_to_analyse, fname);
                continue;
            end
            
            col_data = round(col_data);
            
            % Get time from filename
            % Get time from filename with hard-coded first file
            if strcmp(fname, 'output00000001_cells.mat')
                t = 0;  % First timepoint is always at t=0
            else
                numMatch = regexp(fname, 'output(\d+)_', 'tokens');
                if ~isempty(numMatch)
                    t = str2double(numMatch{1}{1});
                else
                    t = f - 1;  % fallback
                end
            end
            
            sim_times(end+1) = t;

            
            % Count occurrences of each phase/value
            counts = zeros(1, length(phase_values));
            for k = 1:length(phase_values)
                counts(k) = sum(col_data == phase_values(k));
            end
            
            % ALSO extract cell type column (row 6) for same file
            celltype_counts = zeros(1, length(celltype_values));
            if column_celltype <= size(cells, 1)
                celltype_data = cells(column_celltype, :);
                celltype_data = celltype_data(isfinite(celltype_data));
                celltype_data = round(celltype_data);
                for k = 1:length(celltype_values)
                    celltype_counts(k) = sum(celltype_data == celltype_values(k));
                end
            end
            
            % Combine phase and cell type counts for storage
            sim_data(end+1, :) = counts;
            
            % Note: cell type data extracted above for Plot 3
            if f == 1
                celltype_sim_data = []; % Initialize on first file
            end
            
            celltype_sim_data(end+1, :) = celltype_counts;
            file_count = file_count + 1;
            
        catch ME
            fprintf(' Error processing %s: %s (skipping)\n', fname, ME.message);
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
        all_sim_data(valid_sims).celltype_data = celltype_sim_data; % ADD THIS LINE
        all_sim_data(valid_sims).name = sim_folders(s).name;
        
        sim_file_counts(valid_sims) = file_count;
        fprintf(' Successfully loaded %d timepoints\n', file_count);
    else
        fprintf(' No valid data from %s\n', sim_folders(s).name);
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
sem_data = zeros(length(unique_times), num_phases); % Standard error of mean
cv_data = zeros(length(unique_times), num_phases); % Coefficient of variation

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
    sem_data(t_idx, :) = std_data(t_idx, :) / sqrt(valid_sims);
    
    % Coefficient of variation (CV = std/mean * 100%)
    for p = 1:num_phases
        if mean_data(t_idx, p) > 0
            cv_data(t_idx, p) = (std_data(t_idx, p) / mean_data(t_idx, p)) * 100;
        else
            cv_data(t_idx, p) = 0;
        end
    end
end

% ===== PLOT 1: ALIVE PHASES (G1, S, G2, M) =====
try
    fig1 = figure('Position', [100 100 1200 600]);
    plot_with_errorbars(unique_times, mean_data, std_data, alive_phases, alive_labels, 'Alive Cell Phases (G1, S, G2, M)', valid_sims);
    %sgtitle('Plot 1: Alive Phases Over Time', 'FontSize', 14, 'FontWeight', 'bold');
catch ME
    fprintf('Error creating alive phases plot: %s\n', ME.message);
end

% ===== PLOT 2: DEAD PHASES (Apoptosis, Necrosis, Dead1-3) =====
try
    fig2 = figure('Position', [100 750 1200 600]);
    plot_with_errorbars(unique_times, mean_data, std_data, dead_phases, dead_labels, 'Dead Cell Phases (Apoptosis, Necrosis, Dead1-3)', valid_sims);
    %sgtitle('Plot 2: Dead Phases Over Time', 'FontSize', 14, 'FontWeight', 'bold');
catch ME
    fprintf('Error creating dead phases plot: %s\n', ME.message);
end

% ===== PLOT 3: CELL TYPES (6 Types from Your Data) =====
try
    fig3 = figure('Position', [1350 100 1200 600]);
    fprintf('\nProcessing cell type data (row %d)...\n', column_celltype);
    
    celltype_mean_data = zeros(length(unique_times), length(celltype_values));
    celltype_std_data = zeros(length(unique_times), length(celltype_values));
    
    % Collect cell type data for each timepoint across all simulations
    for t_idx = 1:length(unique_times)
        t = unique_times(t_idx);
        celltype_counts_across_sims = [];
        
        for s = 1:valid_sims
            times = all_sim_data(s).times;
            
            % Get cell type data for this simulation (NOT phase data!)
            if isfield(all_sim_data(s), 'celltype_data') && ~isempty(all_sim_data(s).celltype_data)
                celltype_data = all_sim_data(s).celltype_data; % ← Correct: celltype_data
                
                % Find closest timepoint
                [~, closest_idx] = min(abs(times - t));
                
                % Get cell type counts for this timepoint
                if closest_idx <= size(celltype_data, 1)
                    celltype_counts_across_sims(s, :) = celltype_data(closest_idx, :);
                else
                    celltype_counts_across_sims(s, :) = celltype_data(end, :);
                end
            else
                fprintf('Warning: No celltype_data for sim %d at t=%.2f\n', s, t);
            end
        end
        
        % Calculate mean and std across simulations
        if ~isempty(celltype_counts_across_sims)
            celltype_mean_data(t_idx, :) = mean(celltype_counts_across_sims, 1);
            celltype_std_data(t_idx, :) = std(celltype_counts_across_sims, 1);
        end
    end
    
    % Verify data was collected
    if max(max(celltype_mean_data)) == 0
        fprintf('WARNING: celltype_mean_data is all zeros! Checking data structure...\n');
        fprintf(' Number of simulations: %d\n', valid_sims);
        for s = 1:valid_sims
            has_celltype = isfield(all_sim_data(s), 'celltype_data');
            if has_celltype
                celltype_size = size(all_sim_data(s).celltype_data);
                fprintf(' Sim %d: celltype_data exists, size: %d x %d\n', s, celltype_size(1), celltype_size(2));
            else
                fprintf(' Sim %d: NO celltype_data field!\n', s);
            end
        end
    end
    
    % Plot cell types directly (columns 1-6)
    hold on;
    colors = lines(max(length(celltype_values), 3));
    for i = 1:length(celltype_values)
        means = celltype_mean_data(:, i);
        stds = celltype_std_data(:, i);
        
        % Validate data
        if any(~isfinite(means))
            means(~isfinite(means)) = 0;
        end
        if any(~isfinite(stds))
            stds(~isfinite(stds)) = 0;
        end
        
        % Plot line
        plot(unique_times, means, '-', 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', celltype_labels{i});
        hold on;
        
        % Plot error bars
        errorbar(unique_times, means, 1.96 * stds / sqrt(valid_sims), 'k', 'LineStyle', 'none', 'MarkerSize', 0.1, 'CapSize', 3, 'HandleVisibility', 'off');
    end
    
    xlabel('Time (hours)', 'FontSize', 12);
    ylabel('Cell Count by Type', 'FontSize', 12);
    title('Cell Type Distribution Over Time', 'FontSize', 14);
    legend(celltype_labels, 'Location', 'Best', 'FontSize', 11);
    grid on;
    hold off;
    sgtitle('Plot 3: Cell Type Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    
catch ME
    fprintf('Error creating cell type plot: %s\n', ME.message);
    fprintf('Stack: %s\n', ME.stack(1).name);
end

% ===== ENHANCED STATISTICAL ANALYSIS AND REPORT GENERATION =====
generate_enhanced_statistics_report(unique_times, mean_data, std_data, sem_data, cv_data, ...
    phase_labels, phase_values, valid_sims, sim_file_counts, column_to_analyse);

end

% Helper function to plot with error bars
function plot_with_errorbars(times, mean_data, std_data, phase_indices, phase_labels, title_str, valid_sims)

% All 9 phase values in order (matches columns of mean_data)
all_phase_values = [4, 10, 12, 13, 17, 18, 100, 101, 102];

% Validate inputs
if isempty(times) || isempty(mean_data) || isempty(std_data)
    warning('Skipping plot %s: empty data', title_str);
    return;
end

if isempty(phase_indices)
    warning('Skipping plot %s: no phases to plot', title_str);
    return;
end

% Convert phase VALUES to column INDICES
% phase_indices contains phase values like [17, 18, 100, 101, 102]
% We need to find which COLUMN each phase value is in

indices_to_plot = [];
labels_to_plot = {};

for p = 1:length(phase_indices)
    phase_val = phase_indices(p); % e.g. 17 (Apoptosis)
    col_idx = find(all_phase_values == phase_val, 1); % Find which column has this phase
    
    if ~isempty(col_idx) && col_idx <= size(mean_data, 2)
        indices_to_plot = [indices_to_plot, col_idx];
        
        if iscell(phase_labels) && p <= length(phase_labels)
            labels_to_plot = [labels_to_plot, phase_labels(p)];
        else
            labels_to_plot = [labels_to_plot, {num2str(phase_val)}];
        end
    end
end

if isempty(indices_to_plot)
    warning('Skipping plot %s: no valid indices to plot (searched for phase values: %s)', title_str, num2str(phase_indices));
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
            means(~isfinite(means)) = 0;
        end
        if any(~isfinite(stds))
            stds(~isfinite(stds)) = 0;
        end
        
        % Plot line with colour
        plot(times, means, '-', 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', labels_to_plot{i});
        hold on;
        
        % Plot error bars in black
        errorbar(times, means, 1.96 * stds / sqrt(valid_sims), 'k', 'LineStyle', 'none', 'MarkerSize', 0.025, 'CapSize', 3, 'HandleVisibility', 'off');
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

% =========================================================================
% ENHANCED STATISTICS REPORT FUNCTION
% Combines methodology from plot_phases_enhanced.m with plot_simulation_results context
% =========================================================================
function generate_enhanced_statistics_report(times, mean_data, std_data, sem_data, cv_data, ...
    phase_labels, phase_values, num_sims, ~, column_analysed)

% Create filename with timestamp
timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
filename = sprintf('simulation_statistics_%s.txt', timestamp);
fid = fopen(filename, 'w');

% ===== REPORT HEADER =====
fprintf(fid, '================================================================================\n');
fprintf(fid, ' TUMOUR CELL PHASE ANALYSIS REPORT FOR DOSING STRATEGY EVALUATION\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'Report Generated: %s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
fprintf(fid, 'Column Analysed: %d (current_phase)\n\n', column_analysed);

% ===== ANALYSIS METADATA =====
fprintf(fid, '--- ANALYSIS METADATA ---\n');
fprintf(fid, 'Number of Independent Simulations: %d\n', num_sims);
fprintf(fid, 'Time Points Sampled: %d\n', length(times));
fprintf(fid, 'Sampling Interval: %d minutes (%.2f hours)\n', 60, 60/60);
fprintf(fid, 'Total Time Coverage: %.2f hours (%.2f days)\n', times(end), times(end)/24);
fprintf(fid, 'Time Points Range: %.2f - %.2f hours\n\n', times(1), times(end));

fprintf(fid, 'Phases Tracked (9 total):\n');
for i = 1:length(phase_values)
    fprintf(fid, ' [%d] %s (ID: %d)\n', i, phase_labels{i}, phase_values(i));
end
fprintf(fid, '\n');

% ===== GLOBAL SUMMARY STATISTICS =====
fprintf(fid, '================================================================================\n');
fprintf(fid, 'GLOBAL SUMMARY STATISTICS\n');
fprintf(fid, '================================================================================\n\n');

% Calculate total cell population at each timepoint
total_cells_each_time = sum(mean_data, 2);
total_cells_std = sqrt(sum(std_data.^2, 2)); % Sum of variances
total_cells_sem = total_cells_std / sqrt(num_sims);

fprintf(fid, 'TOTAL CELL POPULATION DYNAMICS:\n');
fprintf(fid, ' Initial Population (t=0): %.1f ± %.1f cells\n', total_cells_each_time(1), total_cells_sem(1));
fprintf(fid, ' Final Population (t=%.2f h): %.1f ± %.1f cells\n', times(end), total_cells_each_time(end), total_cells_sem(end));

if total_cells_each_time(1) > 0
    percent_change = ((total_cells_each_time(end) - total_cells_each_time(1)) / total_cells_each_time(1)) * 100;
else
    percent_change = 0;
end

fprintf(fid, ' Net Change: %+.2f%% [%+.1f cells]\n\n', percent_change, total_cells_each_time(end) - total_cells_each_time(1));

% Find minimum and maximum total population
[min_pop, min_idx] = min(total_cells_each_time);
[max_pop, max_idx] = max(total_cells_each_time);

fprintf(fid, ' Minimum Population: %.1f cells at t = %.2f hours\n', min_pop, times(min_idx));
fprintf(fid, ' Maximum Population: %.1f cells at t = %.2f hours\n', max_pop, times(max_idx));

if mean(total_cells_each_time) > 0
    fprintf(fid, ' Population Variability (CV): %.2f%%\n\n', (std(total_cells_each_time)/mean(total_cells_each_time))*100);
end

% ===== PHASE-SPECIFIC STATISTICS =====
fprintf(fid, '================================================================================\n');
fprintf(fid, 'PHASE-SPECIFIC STATISTICS AND ONCOLOGICAL METRICS\n');
fprintf(fid, '================================================================================\n\n');

% Initialize stats table for output
stats_table = table();
stats_table.Phase = cell(length(phase_values), 1);
stats_table.Phase_ID = phase_values';
stats_table.Initial_Count = zeros(length(phase_values), 1);
stats_table.Initial_SE = zeros(length(phase_values), 1);
stats_table.Peak_Count = zeros(length(phase_values), 1);
stats_table.Peak_Time_h = zeros(length(phase_values), 1);
stats_table.Final_Count = zeros(length(phase_values), 1);
stats_table.Final_SE = zeros(length(phase_values), 1);
stats_table.Mean_Count = zeros(length(phase_values), 1);
stats_table.Std_Count = zeros(length(phase_values), 1);
stats_table.AUC_cells_h = zeros(length(phase_values), 1);
stats_table.Growth_Rate_per_h = zeros(length(phase_values), 1);
stats_table.Doubling_Time_h = zeros(length(phase_values), 1);

for phase_idx = 1:length(phase_values)
    phase_val = phase_values(phase_idx);
    phase_name = phase_labels{phase_idx};
    
    % Extract statistics for this phase across all timepoints
    phase_means = mean_data(:, phase_idx);
    phase_stds = std_data(:, phase_idx);
    phase_ses = sem_data(:, phase_idx);
    
    fprintf(fid, '--- %s (Cell Phase ID: %d) ---\n', phase_name, phase_val);
    
    % Kinetic parameters
    init_count = phase_means(1);
    final_count = phase_means(end);
    [peak_count, peak_idx] = max(phase_means);
    peak_time = times(peak_idx);
    
    fprintf(fid, '\nPopulation Dynamics:\n');
    fprintf(fid, ' Initial Count (t=0): %.1f ± %.1f cells\n', init_count, phase_ses(1));
    fprintf(fid, ' Final Count (t=%.2f h): %.1f ± %.1f cells\n', times(end), final_count, phase_ses(end));
    fprintf(fid, ' Peak Count: %.1f cells at t = %.2f hours\n', peak_count, peak_time);
    fprintf(fid, ' Min-Max Range: %.1f - %.1f cells\n', min(phase_means), max(phase_means));
    
    % Calculate growth rate
    if init_count > 0 && final_count > 0 && final_count ~= init_count
        growth_rate = log(final_count / init_count) / (times(end) - times(1));
        doubling_time = log(2) / abs(growth_rate);
        fprintf(fid, ' Net Growth Rate: %.4f per hour\n', growth_rate);
        
        if growth_rate > 0
            fprintf(fid, ' Doubling Time: %.2f hours\n', doubling_time);
        elseif growth_rate < 0
            fprintf(fid, ' Half-life (decline): %.2f hours\n', log(2)/abs(growth_rate));
        end
    else
        growth_rate = 0;
        doubling_time = NaN;
        fprintf(fid, ' Net Growth Rate: Not calculable (zero or steady count)\n');
    end
    
    % Calculate AUC using trapezoidal integration
    auc = trapz(times, phase_means);
    fprintf(fid, '\nAccumulation Metrics:\n');
    fprintf(fid, ' Area Under Curve (AUC): %.1f cell-hours\n', auc);
    fprintf(fid, ' Mean Count Across Time: %.1f ± %.1f cells\n', mean(phase_means), std(phase_means));
    
    if mean(phase_means) > 0
        fprintf(fid, ' Coefficient of Variation: %.2f%%\n', (std(phase_means)/mean(phase_means))*100);
    end
    
    % Phase-specific dosage response metrics
    percent_init = (init_count / total_cells_each_time(1)) * 100;
    percent_final = (final_count / total_cells_each_time(end)) * 100;
    
    fprintf(fid, '\nDosage Response Metrics:\n');
    fprintf(fid, ' Proportion of Initial Population: %.1f%%\n', percent_init);
    fprintf(fid, ' Proportion of Final Population: %.1f%%\n', percent_final);
    fprintf(fid, ' Change in Relative Proportion: %+.1f%%\n', percent_final - percent_init);
    
    % Store in stats table
    stats_table.Phase{phase_idx} = phase_name;
    stats_table.Initial_Count(phase_idx) = init_count;
    stats_table.Initial_SE(phase_idx) = phase_ses(1);
    stats_table.Peak_Count(phase_idx) = peak_count;
    stats_table.Peak_Time_h(phase_idx) = peak_time;
    stats_table.Final_Count(phase_idx) = final_count;
    stats_table.Final_SE(phase_idx) = phase_ses(end);
    stats_table.Mean_Count(phase_idx) = mean(phase_means);
    stats_table.Std_Count(phase_idx) = std(phase_means);
    stats_table.AUC_cells_h(phase_idx) = auc;
    stats_table.Growth_Rate_per_h(phase_idx) = growth_rate;
    stats_table.Doubling_Time_h(phase_idx) = doubling_time;
    
    fprintf(fid, '\n');
end

% ===== CYCLE-SPECIFIC METRICS (ALIVE vs DEAD) =====
fprintf(fid, '================================================================================\n');
fprintf(fid, 'ONCOLOGICAL DOSING STRATEGY EVALUATION METRICS\n');
fprintf(fid, '================================================================================\n\n');

% Identify column indices for cell cycle phases
g1_idx = find(phase_values == 4);
s_idx = find(phase_values == 10);
g2_idx = find(phase_values == 12);
m_idx = find(phase_values == 13);

cycling_cells = zeros(length(times), 1);
if ~isempty(g1_idx)
    cycling_cells = cycling_cells + mean_data(:, g1_idx);
end
if ~isempty(s_idx)
    cycling_cells = cycling_cells + mean_data(:, s_idx);
end
if ~isempty(g2_idx)
    cycling_cells = cycling_cells + mean_data(:, g2_idx);
end
if ~isempty(m_idx)
    cycling_cells = cycling_cells + mean_data(:, m_idx);
end

apoptosis_idx = find(phase_values == 17);
necrotic_idx = find(phase_values == 18);
dead1_idx = find(phase_values == 100);
dead2_idx = find(phase_values == 101);
dead3_idx = find(phase_values == 102);

apoptotic_cells = zeros(length(times), 1);
if ~isempty(apoptosis_idx)
    apoptotic_cells = mean_data(:, apoptosis_idx);
end

necrotic_cells = zeros(length(times), 1);
if ~isempty(necrotic_idx)
    necrotic_cells = mean_data(:, necrotic_idx);
end

dead_cells = zeros(length(times), 1);
if ~isempty(dead1_idx)
    dead_cells = dead_cells + mean_data(:, dead1_idx);
end
if ~isempty(dead2_idx)
    dead_cells = dead_cells + mean_data(:, dead2_idx);
end
if ~isempty(dead3_idx)
    dead_cells = dead_cells + mean_data(:, dead3_idx);
end

fprintf(fid, 'CELL CYCLE KINETICS (G1 + S + G2 + M):\n');
fprintf(fid, ' Initial Cycling Cells: %.1f cells (%.1f%% of population)\n', ...
    cycling_cells(1), (cycling_cells(1)/total_cells_each_time(1))*100);
fprintf(fid, ' Final Cycling Cells: %.1f cells (%.1f%% of population)\n', ...
    cycling_cells(end), (cycling_cells(end)/total_cells_each_time(end))*100);
fprintf(fid, ' Mean Cycling Cells: %.1f ± %.1f cells\n', mean(cycling_cells), std(cycling_cells));

fprintf(fid, '\nAPOPTOSIS & NECROSIS (Drug Response Indicators):\n');
fprintf(fid, ' Total Apoptotic Cells: %.1f cells (AUC: %.1f cell-hours)\n', ...
    apoptotic_cells(end), trapz(times, apoptotic_cells));
fprintf(fid, ' Total Necrotic Cells: %.1f cells (AUC: %.1f cell-hours)\n', ...
    necrotic_cells(end), trapz(times, necrotic_cells));
fprintf(fid, ' Total Dead Cells: %.1f cells\n', dead_cells(end));

total_death = apoptotic_cells(end) + necrotic_cells(end);
if times(end) > 0
    fprintf(fid, ' Cell Death Rate (apoptosis+necrosis): %.2f cells/hour\n', total_death / times(end));
end

fprintf(fid, '\nTHERAPEUTIC EFFICACY INDICATORS:\n');
fprintf(fid, ' Total Cell Elimination: %.1f cells\n', total_death);

if total_cells_each_time(1) > 0
    efficacy_index = (total_death / total_cells_each_time(1)) * 100;
    fprintf(fid, ' Efficacy Index (% of initial pop): %.1f%%\n', efficacy_index);
end

if total_cells_each_time(end) < total_cells_each_time(1)
    fprintf(fid, ' Tumour Control Status: RESPONSIVE (net reduction)\n');
else
    fprintf(fid, ' Tumour Control Status: PROGRESSIVE (net growth)\n');
end

% ===== DETAILED TIME-SERIES TABLE =====
fprintf(fid, '\n================================================================================\n');
fprintf(fid, 'TIME-SERIES STATISTICS TABLE\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'Time(h)\t');
for p = 1:length(phase_labels)
    fprintf(fid, '%s_Mean\t%s_SE\t%s_CV(%%))\t', phase_labels{p}, phase_labels{p}, phase_labels{p});
end
fprintf(fid, '\n');

for t = 1:length(times)
    fprintf(fid, '%.2f\t', times(t));
    for p = 1:length(phase_labels)
        fprintf(fid, '%.1f\t%.3f\t%.2f\t', mean_data(t, p), sem_data(t, p), cv_data(t, p));
    end
    fprintf(fid, '\n');
end

% ===== CONCLUSIONS =====
fprintf(fid, '\n================================================================================\n');
fprintf(fid, 'CONCLUSIONS AND RECOMMENDATIONS\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'SUMMARY OF KEY FINDINGS:\n\n');

fprintf(fid, '1. TUMOUR CONTROL:\n');
if total_cells_each_time(end) < total_cells_each_time(1) * 0.5
    fprintf(fid, ' - Excellent tumour control achieved (>50%% reduction in cell population)\n');
elseif total_cells_each_time(end) < total_cells_each_time(1)
    fprintf(fid, ' - Moderate tumour control achieved (net reduction in cell population)\n');
else
    fprintf(fid, ' - Limited tumour control (cell population growth observed)\n');
end

fprintf(fid, '\n2. CELL CYCLE DISRUPTION:\n');
if cycling_cells(end) < cycling_cells(1) * 0.3
    fprintf(fid, ' - Strong cell cycle arrest achieved\n');
elseif cycling_cells(end) < cycling_cells(1) * 0.7
    fprintf(fid, ' - Partial cell cycle arrest observed\n');
else
    fprintf(fid, ' - Minimal cell cycle disruption\n');
end

fprintf(fid, '\n3. CELL DEATH INDUCTION:\n');
if total_cells_each_time(1) > 0 && total_death > total_cells_each_time(1) * 0.3
    fprintf(fid, ' - Strong apoptotic/necrotic response (>30%% cell death)\n');
elseif total_cells_each_time(1) > 0 && total_death > total_cells_each_time(1) * 0.1
    fprintf(fid, ' - Moderate cell death induction (10-30%%)\n');
else
    fprintf(fid, ' - Limited cell death response\n');
end

fprintf(fid, '\n4. DOSING RECOMMENDATIONS:\n');
if total_cells_each_time(end) < total_cells_each_time(1)
    fprintf(fid, ' - Current dosing shows INHIBITORY response profile\n');
else
    fprintf(fid, ' - Current dosing shows STIMULATORY response profile\n');
end

% Find dominant phases
[~, rank_idx] = sort(stats_table.AUC_cells_h(1:4), 'descend');
if ~isempty(rank_idx)
    fprintf(fid, ' - Most responsive phases (by AUC): %s, %s\n', ...
        stats_table.Phase{rank_idx(1)}, stats_table.Phase{rank_idx(2)});
    fprintf(fid, ' - Focus on targeting cells in %s phase(s) for optimisation\n', ...
        stats_table.Phase{rank_idx(1)});
end

fprintf(fid, '\n================================================================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '================================================================================\n');

fclose(fid);

% Display to console
fprintf('\n*** STATISTICAL ANALYSIS REPORT GENERATED ***\n');
fprintf('Report saved to: %s\n', filename);
fprintf('File size: %.1f KB\n', dir(filename).bytes / 1024);

fprintf('\n=== QUICK SUMMARY ===\n');
fprintf('Initial population: %.1f cells\n', total_cells_each_time(1));
fprintf('Final population: %.1f cells\n', total_cells_each_time(end));
fprintf('Total cell deaths (apoptosis+necrosis): %.1f cells\n', total_death);
if total_cells_each_time(1) > 0
    fprintf('Efficacy index: %.1f%% of initial population eliminated\n', (total_death/total_cells_each_time(1))*100);
end

end