function plot_simulation_single_run()

% =========================================================================
% SINGLE RUN ANALYSIS - Adapted from plot_simulation_results.m
% Processes ONE simulation output folder and generates statistics
% CSV EXPORT: Outputs all plotted data to CSV file for further analysis
% =========================================================================

% Configuration
column_to_analyse = 8; % current_phase
interval_min = 60; % Time interval in minutes
output_folder = 'output'; % Single output folder (not sim_*)

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
column_celltype = 6;
celltype_values = [0, 1, 2, 3, 4, 5];
celltype_labels = {'Cancer_Stem', 'Cancer_Prolif', 'Cancer_Diff', 'CAF', 'CD8_T_Cell', 'M2_Macro'};

% Check if folder exists
if ~isfolder(output_folder)
    error('Output folder not found: %s', output_folder);
end

% Get all output*_cells.mat files
files = dir(fullfile(output_folder, 'output*_cells.mat'));
if isempty(files)
    error('No *_cells.mat files found in %s', output_folder);
end

fprintf('Found %d output files in %s\n', length(files), output_folder);

% Load all data
times = [];
sim_data = [];
celltype_sim_data = [];

for f = 1:length(files)
    fname = files(f).name;
    filepath = fullfile(output_folder, fname);
    
    try
        load_data = load(filepath, 'cells');
        if ~isfield(load_data, 'cells')
            fprintf('Warning: No cells field in %s (skipping)\n', fname);
            continue;
        end
        
        cells = load_data.cells;
        
        % Validate cells array
        if ~isnumeric(cells) || isempty(cells)
            fprintf('Warning: Invalid cells array in %s (skipping)\n', fname);
            continue;
        end
        
        % Check if column exists
        if column_to_analyse > size(cells, 1)
            fprintf('Warning: Column %d does not exist in %s (skipping)\n', column_to_analyse, fname);
            continue;
        end
        
        % Extract phase data
        col_data = cells(column_to_analyse, :);
        col_data = col_data(isfinite(col_data));
        if isempty(col_data)
            fprintf('Warning: No valid data in column %d of %s (skipping)\n', column_to_analyse, fname);
            continue;
        end
        
        col_data = round(col_data);
        
        % Get time from filename
        if strcmp(fname, 'output00000001_cells.mat')
            t = 0;
        else
            numMatch = regexp(fname, 'output(\d+)_', 'tokens');
            if ~isempty(numMatch)
                t = str2double(numMatch{1}{1});
            else
                t = f - 1;
            end
        end
        
        times(f) = t;
        
        % Count occurrences of each phase
        counts = zeros(1, length(phase_values));
        for k = 1:length(phase_values)
            counts(k) = sum(col_data == phase_values(k));
        end
        
        sim_data(f, :) = counts;
        
        % Extract cell type data
        celltype_counts = zeros(1, length(celltype_values));
        if column_celltype <= size(cells, 1)
            celltype_data = cells(column_celltype, :);
            celltype_data = celltype_data(isfinite(celltype_data));
            celltype_data = round(celltype_data);
            for k = 1:length(celltype_values)
                celltype_counts(k) = sum(celltype_data == celltype_values(k));
            end
        end
        
        celltype_sim_data(f, :) = celltype_counts;
        
    catch ME
        fprintf('Error processing %s: %s (skipping)\n', fname, ME.message);
        continue;
    end
end

% Validate we have data
if isempty(times)
    error('No valid data loaded');
end

% Convert times to hours
times = times * interval_min / 60;
times = times(1:length(sim_data));

fprintf('Successfully loaded %d timepoints\n', length(times));
fprintf('Time range: %.2f to %.2f hours\n\n', times(1), times(end));

% ===== EXPORT DATA TO CSV =====
export_simulation_data_to_csv(times, sim_data, celltype_sim_data, phase_labels, celltype_labels);

% ===== PLOT 1: ALIVE PHASES =====
try
    fig1 = figure('Position', [100 100 1200 600]);
    plot_single_run(times, sim_data, alive_phases, alive_labels, 'Alive Cell Phases (G1, S, G2, M)');
catch ME
    fprintf('Error creating alive phases plot: %s\n', ME.message);
end

% ===== PLOT 2: DEAD PHASES =====
try
    fig2 = figure('Position', [100 750 1200 600]);
    plot_single_run(times, sim_data, dead_phases, dead_labels, 'Dead Cell Phases (Apoptosis, Necrosis, Dead1-3)');
catch ME
    fprintf('Error creating dead phases plot: %s\n', ME.message);
end

% ===== PLOT 3: CELL TYPES =====
try
    fig3 = figure('Position', [1350 100 1200 600]);
    hold on;
    colors = lines(length(celltype_values));
    for i = 1:length(celltype_values)
        plot(times, celltype_sim_data(:, i), '-o', 'Color', colors(i, :), 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', celltype_labels{i});
        hold on;
    end
    xlabel('Time (hours)', 'FontSize', 12);
    ylabel('Cell Count by Type', 'FontSize', 12);
    title('Cell Type Distribution Over Time', 'FontSize', 14);
    legend(celltype_labels, 'Location', 'Best', 'FontSize', 11);
    grid on;
    hold off;
catch ME
    fprintf('Error creating cell type plot: %s\n', ME.message);
end

% ===== GENERATE STATISTICS REPORT =====
generate_statistics_report_single(times, sim_data, celltype_sim_data, phase_labels, phase_values, column_to_analyse);

end

% =========================================================================
% CSV EXPORT FUNCTION: Export all plotted data per timepoint
% =========================================================================
function export_simulation_data_to_csv(times, sim_data, celltype_sim_data, phase_labels, celltype_labels)

timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
csv_filename = sprintf('simulation_data_export_%s.csv', timestamp);

% Open CSV file for writing
fid = fopen(csv_filename, 'w');

% Write header
fprintf(fid, 'Time_Hours');
for i = 1:length(phase_labels)
    fprintf(fid, ',%s', phase_labels{i});
end
for i = 1:length(celltype_labels)
    fprintf(fid, ',%s', celltype_labels{i});
end
fprintf(fid, ',Total_Cells\n');

% Calculate total cells per timepoint
total_cells_each_time = sum(sim_data, 2);

% Write data for each timepoint
for t = 1:length(times)
    fprintf(fid, '%.2f', times(t));
    
    % Write phase data
    for p = 1:length(phase_labels)
        fprintf(fid, ',%.0f', sim_data(t, p));
    end
    
    % Write cell type data
    for c = 1:length(celltype_labels)
        fprintf(fid, ',%.0f', celltype_sim_data(t, c));
    end
    
    % Write total cells
    fprintf(fid, ',%.0f\n', total_cells_each_time(t));
end

fclose(fid);

fprintf('\n*** CSV EXPORT COMPLETED ***\n');
fprintf('Data exported to: %s\n', csv_filename);
fprintf('Timepoints: %d\n', length(times));
fprintf('Variables: %d cell phases + %d cell types + 1 total\n', length(phase_labels), length(celltype_labels));

end

% =========================================================================
% HELPER FUNCTION: Plot single run (no error bars, just line plot)
% =========================================================================
function plot_single_run(times, sim_data, phase_indices, phase_labels, title_str)

all_phase_values = [4, 10, 12, 13, 17, 18, 100, 101, 102];

% Convert phase values to column indices
indices_to_plot = [];
labels_to_plot = {};

for p = 1:length(phase_indices)
    phase_val = phase_indices(p);
    col_idx = find(all_phase_values == phase_val, 1);
    if ~isempty(col_idx) && col_idx <= size(sim_data, 2)
        indices_to_plot = [indices_to_plot, col_idx];
        if iscell(phase_labels) && p <= length(phase_labels)
            labels_to_plot = [labels_to_plot, phase_labels(p)];
        else
            labels_to_plot = [labels_to_plot, {num2str(phase_val)}];
        end
    end
end

if isempty(indices_to_plot)
    warning('No valid indices to plot');
    return;
end

% Plot each phase (no error bars)
hold on;
num_plots = length(indices_to_plot);
colors = lines(max(num_plots, 3));

for i = 1:num_plots
    try
        idx = indices_to_plot(i);
        vals = sim_data(:, idx);
        
        % Validate data
        if any(~isfinite(vals))
            vals(~isfinite(vals)) = 0;
        end
        
        % Plot line with markers
        plot(times, vals, '-o', 'Color', colors(i, :), 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', labels_to_plot{i});
        hold on;
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
% STATISTICS REPORT GENERATION FOR SINGLE RUN
% =========================================================================
function generate_statistics_report_single(times, sim_data, celltype_sim_data, phase_labels, phase_values, column_analysed)

timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
filename = sprintf('simulation_statistics_%s.txt', timestamp);

fid = fopen(filename, 'w');

% ===== REPORT HEADER =====
fprintf(fid, '================================================================================\n');
fprintf(fid, ' TUMOUR CELL PHASE ANALYSIS REPORT - SINGLE RUN\n');
fprintf(fid, '================================================================================\n\n');
fprintf(fid, 'Report Generated: %s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
fprintf(fid, 'Column Analysed: %d (current_phase)\n\n', column_analysed);

% ===== ANALYSIS METADATA =====
fprintf(fid, '--- ANALYSIS METADATA ---\n');
fprintf(fid, 'Analysis Type: Single Run (No Replicates)\n');
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

% Total population at each timepoint
total_cells_each_time = sum(sim_data, 2);

fprintf(fid, 'TOTAL CELL POPULATION DYNAMICS:\n');
fprintf(fid, ' Initial Population (t=0): %.1f cells\n', total_cells_each_time(1));
fprintf(fid, ' Final Population (t=%.2f h): %.1f cells\n', times(end), total_cells_each_time(end));

if total_cells_each_time(1) > 0
    percent_change = ((total_cells_each_time(end) - total_cells_each_time(1)) / total_cells_each_time(1)) * 100;
else
    percent_change = 0;
end

fprintf(fid, ' Net Change: %+.2f%% [%+.1f cells]\n\n', percent_change, total_cells_each_time(end) - total_cells_each_time(1));

% Min/max population
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

% Create stats table
stats_table = table();
stats_table.Phase = cell(length(phase_values), 1);
stats_table.Phase_ID = phase_values';
stats_table.Initial_Count = zeros(length(phase_values), 1);
stats_table.Peak_Count = zeros(length(phase_values), 1);
stats_table.Peak_Time_h = zeros(length(phase_values), 1);
stats_table.Final_Count = zeros(length(phase_values), 1);
stats_table.Mean_Count = zeros(length(phase_values), 1);
stats_table.Min_Count = zeros(length(phase_values), 1);
stats_table.Max_Count = zeros(length(phase_values), 1);
stats_table.AUC_cells_h = zeros(length(phase_values), 1);
stats_table.Growth_Rate_per_h = zeros(length(phase_values), 1);
stats_table.Doubling_Time_h = zeros(length(phase_values), 1);

for phase_idx = 1:length(phase_values)
    phase_val = phase_values(phase_idx);
    phase_name = phase_labels{phase_idx};
    
    % Extract data for this phase
    phase_vals = sim_data(:, phase_idx);
    
    fprintf(fid, '--- %s (Cell Phase ID: %d) ---\n', phase_name, phase_val);
    
    % Kinetic parameters
    init_count = phase_vals(1);
    final_count = phase_vals(end);
    [peak_count, peak_idx] = max(phase_vals);
    peak_time = times(peak_idx);
    min_count = min(phase_vals);
    max_count = max(phase_vals);
    mean_count = mean(phase_vals);
    
    fprintf(fid, '\nPopulation Dynamics:\n');
    fprintf(fid, ' Initial Count (t=0): %.1f cells\n', init_count);
    fprintf(fid, ' Final Count (t=%.2f h): %.1f cells\n', times(end), final_count);
    fprintf(fid, ' Peak Count: %.1f cells at t = %.2f hours\n', peak_count, peak_time);
    fprintf(fid, ' Min-Max Range: %.1f - %.1f cells\n', min_count, max_count);
    fprintf(fid, ' Mean Count: %.1f cells\n', mean_count);
    
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
    auc = trapz(times, phase_vals);
    fprintf(fid, ' Area Under Curve (AUC): %.1f cell-hours\n', auc);
    
    % Phase proportion
    percent_init = (init_count / total_cells_each_time(1)) * 100;
    percent_final = (final_count / total_cells_each_time(end)) * 100;
    
    fprintf(fid, '\nDosage Response Metrics:\n');
    fprintf(fid, ' Initial Proportion: %.1f%% of population\n', percent_init);
    fprintf(fid, ' Final Proportion: %.1f%% of population\n', percent_final);
    fprintf(fid, ' Change in Proportion: %+.1f%%\n', percent_final - percent_init);
    
    % Store in table
    stats_table.Phase{phase_idx} = phase_name;
    stats_table.Initial_Count(phase_idx) = init_count;
    stats_table.Peak_Count(phase_idx) = peak_count;
    stats_table.Peak_Time_h(phase_idx) = peak_time;
    stats_table.Final_Count(phase_idx) = final_count;
    stats_table.Mean_Count(phase_idx) = mean_count;
    stats_table.Min_Count(phase_idx) = min_count;
    stats_table.Max_Count(phase_idx) = max_count;
    stats_table.AUC_cells_h(phase_idx) = auc;
    stats_table.Growth_Rate_per_h(phase_idx) = growth_rate;
    stats_table.Doubling_Time_h(phase_idx) = doubling_time;
    
    fprintf(fid, '\n');
end

% ===== CELL CYCLE & DEATH ANALYSIS =====
fprintf(fid, '================================================================================\n');
fprintf(fid, 'ONCOLOGICAL DOSING STRATEGY EVALUATION METRICS\n');
fprintf(fid, '================================================================================\n\n');

% Find indices
g1_idx = find(phase_values == 4);
s_idx = find(phase_values == 10);
g2_idx = find(phase_values == 12);
m_idx = find(phase_values == 13);

cycling_cells = zeros(length(times), 1);
if ~isempty(g1_idx), cycling_cells = cycling_cells + sim_data(:, g1_idx); end
if ~isempty(s_idx), cycling_cells = cycling_cells + sim_data(:, s_idx); end
if ~isempty(g2_idx), cycling_cells = cycling_cells + sim_data(:, g2_idx); end
if ~isempty(m_idx), cycling_cells = cycling_cells + sim_data(:, m_idx); end

apoptosis_idx = find(phase_values == 17);
necrotic_idx = find(phase_values == 18);
dead1_idx = find(phase_values == 100);
dead2_idx = find(phase_values == 101);
dead3_idx = find(phase_values == 102);

apoptotic_cells = zeros(length(times), 1);
if ~isempty(apoptosis_idx), apoptotic_cells = sim_data(:, apoptosis_idx); end

necrotic_cells = zeros(length(times), 1);
if ~isempty(necrotic_idx), necrotic_cells = sim_data(:, necrotic_idx); end

dead_cells = zeros(length(times), 1);
if ~isempty(dead1_idx), dead_cells = dead_cells + sim_data(:, dead1_idx); end
if ~isempty(dead2_idx), dead_cells = dead_cells + sim_data(:, dead2_idx); end
if ~isempty(dead3_idx), dead_cells = dead_cells + sim_data(:, dead3_idx); end

fprintf(fid, 'CELL CYCLE KINETICS (G1 + S + G2 + M):\n');
fprintf(fid, ' Initial Cycling Cells: %.1f cells (%.1f%% of population)\n', cycling_cells(1), (cycling_cells(1)/total_cells_each_time(1))*100);
fprintf(fid, ' Final Cycling Cells: %.1f cells (%.1f%% of population)\n', cycling_cells(end), (cycling_cells(end)/total_cells_each_time(end))*100);
fprintf(fid, ' Mean Cycling Cells: %.1f cells\n', mean(cycling_cells));
fprintf(fid, ' Peak Cycling Cells: %.1f cells at t = %.2f hours\n', max(cycling_cells), times(find(cycling_cells == max(cycling_cells), 1)));

fprintf(fid, '\nAPOPTOSIS & NECROSIS (Drug Response Indicators):\n');
fprintf(fid, ' Total Apoptotic Cells (final): %.1f cells\n', apoptotic_cells(end));
fprintf(fid, ' Total Apoptotic AUC: %.1f cell-hours\n', trapz(times, apoptotic_cells));
fprintf(fid, ' Peak Apoptosis: %.1f cells at t = %.2f hours\n', max(apoptotic_cells), times(find(apoptotic_cells == max(apoptotic_cells), 1)));
fprintf(fid, ' Total Necrotic Cells (final): %.1f cells\n', necrotic_cells(end));
fprintf(fid, ' Total Necrotic AUC: %.1f cell-hours\n', trapz(times, necrotic_cells));
fprintf(fid, ' Total Dead Cells: %.1f cells\n', dead_cells(end));

total_death = apoptotic_cells(end) + necrotic_cells(end);

fprintf(fid, '\nTHERAPEUTIC EFFICACY INDICATORS:\n');
fprintf(fid, ' Total Cell Deaths (apoptosis+necrosis): %.1f cells\n', total_death);

if total_cells_each_time(1) > 0
    efficacy_index = (total_death / total_cells_each_time(1)) * 100;
    fprintf(fid, ' Efficacy Index (cell death as %% of initial): %.1f%%\n', efficacy_index);
end

if times(end) > 0
    fprintf(fid, ' Cell Death Rate: %.4f cells/hour\n', total_death / times(end));
end

if total_cells_each_time(end) < total_cells_each_time(1)
    fprintf(fid, ' Tumour Control Status: RESPONSIVE (net reduction in population)\n');
else
    fprintf(fid, ' Tumour Control Status: PROGRESSIVE (net growth in population)\n');
end

% ===== TIME-SERIES TABLE =====
fprintf(fid, '\n================================================================================\n');
fprintf(fid, 'TIME-SERIES DATA TABLE\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'Time(h)\t');
for p = 1:length(phase_labels)
    fprintf(fid, '%s\t', phase_labels{p});
end
fprintf(fid, '\n');

for t = 1:length(times)
    fprintf(fid, '%.2f\t', times(t));
    for p = 1:length(phase_labels)
        fprintf(fid, '%.0f\t', sim_data(t, p));
    end
    fprintf(fid, '\n');
end

% ===== CELL TYPE SUMMARY =====
if ~isempty(celltype_sim_data)
    fprintf(fid, '\n================================================================================\n');
    fprintf(fid, 'CELL TYPE DISTRIBUTION\n');
    fprintf(fid, '================================================================================\n\n');
    
    celltype_labels = {'Cancer_Stem', 'Cancer_Prolif', 'Cancer_Diff', 'CAF', 'CD8_T', 'M2_Macro'};
    
    fprintf(fid, 'Initial Cell Type Composition:\n');
    for i = 1:size(celltype_sim_data, 2)
        if i <= length(celltype_labels)
            pct = (celltype_sim_data(1, i) / sum(celltype_sim_data(1, :))) * 100;
            fprintf(fid, ' %s: %.1f cells (%.1f%%)\n', celltype_labels{i}, celltype_sim_data(1, i), pct);
        end
    end
    
    fprintf(fid, '\nFinal Cell Type Composition:\n');
    for i = 1:size(celltype_sim_data, 2)
        if i <= length(celltype_labels)
            pct = (celltype_sim_data(end, i) / sum(celltype_sim_data(end, :))) * 100;
            fprintf(fid, ' %s: %.1f cells (%.1f%%)\n', celltype_labels{i}, celltype_sim_data(end, i), pct);
        end
    end
end

% ===== CONCLUSIONS =====
fprintf(fid, '\n================================================================================\n');
fprintf(fid, 'CONCLUSIONS AND RECOMMENDATIONS\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'SUMMARY OF KEY FINDINGS:\n\n');

fprintf(fid, '1. TUMOUR CONTROL:\n');
if total_cells_each_time(end) < total_cells_each_time(1) * 0.5
    fprintf(fid, ' - Excellent tumour control achieved (>50%% reduction)\n');
elseif total_cells_each_time(end) < total_cells_each_time(1)
    fprintf(fid, ' - Moderate tumour control achieved (net reduction)\n');
else
    fprintf(fid, ' - Limited tumour control (population growth observed)\n');
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

fprintf(fid, '\n================================================================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '================================================================================\n');

fclose(fid);

% Display summary to console
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

fprintf('Report duration: %.2f hours (%.2f days)\n', times(end), times(end)/24);

end

