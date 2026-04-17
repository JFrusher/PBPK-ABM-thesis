function analyze_batch_simulations(batch_folder)
% =========================================================================
% BATCH ANALYSIS - Process multiple MC_run_XXXX directories
% Analyzes all output folders within a parent batch directory
% Aggregates results into a single CSV file
% =========================================================================
%
% Usage:
%   analyze_batch_simulations('./mc_results/my_batch')
%   analyze_batch_simulations('D:\PhysiCell\mc_results\batch_2025_01_11')
%
% Input:
%   batch_folder - Path to folder containing MC_run_XXXX subdirectories
%
% Output:
%   batch_summary_TIMESTAMP.csv - Aggregated results from all runs
%
% =========================================================================

% Handle input
if nargin < 1
    batch_folder = pwd;
    fprintf('No batch folder specified. Using current directory: %s\n\n', batch_folder);
end

% Ensure folder exists
if ~isfolder(batch_folder)
    error('Batch folder not found: %s', batch_folder);
end

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('PhysiCell Batch Analysis - Process Multiple Runs\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Configuration (same as single run)
column_to_analyse = 8; % current_phase
interval_min = 60; % Time interval in minutes
phase_values = [4, 10, 12, 13, 17, 18, 100, 101, 102];
phase_labels = {'G1', 'S', 'G2', 'M', 'Apoptosis', 'Necrosis', 'Dead1', 'Dead2', 'Dead3'};
alive_phases = [4, 10, 12, 13];
dead_phases = [17, 18, 100, 101, 102];
celltype_values = [0, 1, 2, 3, 4, 5];
celltype_labels = {'Cancer_Stem', 'Cancer_Prolif', 'Cancer_Diff', 'CAF', 'CD8_T_Cell', 'M2_Macro'};

% Find all MC_run_XXXX directories
run_dirs = dir(fullfile(batch_folder, 'MC_run_*'));
run_dirs = run_dirs([run_dirs.isdir]);

if isempty(run_dirs)
    error('No MC_run_* directories found in: %s', batch_folder);
end

% Sort by numeric index
[~, idx] = sort_run_numbers({run_dirs.name});
run_dirs = run_dirs(idx);

fprintf('Found %d MC_run directories:\n', length(run_dirs));
for i = 1:length(run_dirs)
    fprintf('  [%2d] %s\n', i, run_dirs(i).name);
end
fprintf('\n');

% =========================================================================
% PROCESS EACH RUN
% =========================================================================

all_results = [];
batch_summary_data = {};
batch_idx = 1;

for run_num = 1:length(run_dirs)
    
    run_name = run_dirs(run_num).name;
    run_path = fullfile(batch_folder, run_name);
    
    % Extract numeric index from MC_run_XXX
    idx_match = regexp(run_name, 'MC_run_(\d+)', 'tokens');
    if ~isempty(idx_match)
        run_idx = str2double(idx_match{1}{1});
    else
        run_idx = run_num;
    end
    
    % Look for output folder or physicell_output folder
    output_path = fullfile(run_path, 'output');
    if ~isfolder(output_path)
        output_path = fullfile(run_path, 'physicell_output');
    end
    
    if ~isfolder(output_path)
        fprintf('[%2d/%2d] %s - ⚠ No output folder found (SKIP)\n', run_num, length(run_dirs), run_name);
        continue;
    end
    
    fprintf('[%2d/%2d] %s - Processing... ', run_num, length(run_dirs), run_name);
    
    % Get all output*_cells.mat files
    files = dir(fullfile(output_path, 'output*_cells.mat'));
    
    if isempty(files)
        fprintf('✗ No data files found\n');
        continue;
    end
    
    % Load and analyze this run
    try
        [times, sim_data, celltype_sim_data, summary_stats] = ...
            load_and_analyze_single_run(output_path, files, column_to_analyse, ...
            interval_min, phase_values, celltype_values);
        
        fprintf('✓ Processed %d timepoints\n', length(times));
        
        % Store results
        batch_summary_data{batch_idx, 1} = run_idx;
        batch_summary_data{batch_idx, 2} = run_name;
        batch_summary_data{batch_idx, 3} = length(times);
        batch_summary_data{batch_idx, 4} = times(end);  % Final time
        batch_summary_data{batch_idx, 5} = summary_stats.final_alive;
        batch_summary_data{batch_idx, 6} = summary_stats.final_dead;
        batch_summary_data{batch_idx, 7} = summary_stats.final_total;
        batch_summary_data{batch_idx, 8} = summary_stats.max_alive;
        batch_summary_data{batch_idx, 9} = summary_stats.min_alive;
        
        batch_idx = batch_idx + 1;
        
        % Store detailed results
        all_results(run_num).run_idx = run_idx;
        all_results(run_num).run_name = run_name;
        all_results(run_num).times = times;
        all_results(run_num).sim_data = sim_data;
        all_results(run_num).celltype_sim_data = celltype_sim_data;
        all_results(run_num).summary_stats = summary_stats;
        
    catch ME
        fprintf('✗ Error: %s\n', ME.message);
        continue;
    end
    
end

fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('Processing Complete: %d runs analyzed\n', batch_idx - 1);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% =========================================================================
% EXPORT SUMMARY CSV
% =========================================================================

if ~isempty(batch_summary_data)
    
    % Create summary table
    summary_table = cell2table(batch_summary_data);
    summary_table.Properties.VariableNames = {...
        'Run_Index', 'Run_Name', 'Num_Timepoints', 'Final_Time_Hours', ...
        'Final_Alive_Cells', 'Final_Dead_Cells', 'Final_Total_Cells', ...
        'Max_Alive_Cells', 'Min_Alive_Cells'};
    
    % Export to CSV
    timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
    csv_filename = fullfile(batch_folder, sprintf('batch_summary_%s.csv', timestamp));
    
    writetable(summary_table, csv_filename);
    fprintf('✓ Summary CSV exported: %s\n', csv_filename);
    
    % Display summary table
    fprintf('\n');
    disp(summary_table);
    fprintf('\n');
end

% =========================================================================
% CREATE COMPARISON PLOTS
% =========================================================================

if ~isempty(all_results)
    create_batch_comparison_plots(all_results, phase_labels, celltype_labels);
end

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('Analysis Complete!\n');
fprintf('═══════════════════════════════════════════════════════════════\n');

end

% =========================================================================
% HELPER FUNCTION: Load and analyze single run
% =========================================================================
function [times, sim_data, celltype_sim_data, summary_stats] = ...
    load_and_analyze_single_run(output_path, files, column_to_analyse, ...
    interval_min, phase_values, celltype_values)

times = [];
sim_data = [];
celltype_sim_data = [];
column_celltype = 6;

for f = 1:length(files)
    fname = files(f).name;
    filepath = fullfile(output_path, fname);
    
    try
        load_data = load(filepath, 'cells');
        if ~isfield(load_data, 'cells')
            continue;
        end
        
        cells = load_data.cells;
        
        if ~isnumeric(cells) || isempty(cells)
            continue;
        end
        
        if column_to_analyse > size(cells, 1)
            continue;
        end
        
        col_data = cells(column_to_analyse, :);
        col_data = col_data(isfinite(col_data));
        
        if isempty(col_data)
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
        % Skip this file silently
        continue;
    end
end

% Convert times to hours
if ~isempty(times)
    times = times * interval_min / 60;
    times = times(1:length(sim_data));
end

% Calculate summary statistics
alive_phases = [4, 10, 12, 13];
dead_phases = [17, 18, 100, 101, 102];

summary_stats.final_alive = sum(sim_data(end, ismember(phase_values, alive_phases)));
summary_stats.final_dead = sum(sim_data(end, ismember(phase_values, dead_phases)));
summary_stats.final_total = sum(sim_data(end, :));
summary_stats.max_alive = max(sum(sim_data(:, ismember(phase_values, alive_phases)), 2));
summary_stats.min_alive = min(sum(sim_data(:, ismember(phase_values, alive_phases)), 2));

end

% =========================================================================
% HELPER FUNCTION: Sort run numbers naturally
% =========================================================================
function [sorted_names, idx] = sort_run_numbers(names)

numbers = [];
for i = 1:length(names)
    match = regexp(names{i}, 'MC_run_(\d+)', 'tokens');
    if ~isempty(match)
        numbers(i) = str2double(match{1}{1});
    else
        numbers(i) = i;
    end
end

[~, idx] = sort(numbers);
sorted_names = names(idx);

end

% =========================================================================
% HELPER FUNCTION: Create batch comparison plots
% =========================================================================
function create_batch_comparison_plots(all_results, phase_labels, celltype_labels)

fprintf('\nGenerating comparison plots...\n');

% Extract final statistics for each run
run_names = {all_results.run_name};
num_runs = length(all_results);
final_alive = [];
final_total = [];

for i = 1:num_runs
    final_alive(i) = all_results(i).summary_stats.final_alive;
    final_total(i) = all_results(i).summary_stats.final_total;
end

% Create comparison figure
fig = figure('Position', [100 100 1200 600]);

% Plot 1: Final alive cells per run
subplot(1, 2, 1);
bar(1:num_runs, final_alive, 'FaceColor', [0.2 0.8 0.2]);
xlabel('Run Number', 'FontSize', 12);
ylabel('Final Alive Cells', 'FontSize', 12);
title('Final Alive Cell Count Across Runs', 'FontSize', 14);
grid on;

% Plot 2: Final total cells per run
subplot(1, 2, 2);
bar(1:num_runs, final_total, 'FaceColor', [0.2 0.4 0.8]);
xlabel('Run Number', 'FontSize', 12);
ylabel('Final Total Cells', 'FontSize', 12);
title('Final Total Cell Count Across Runs', 'FontSize', 14);
grid on;

% Save figure
timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
savefig(sprintf('batch_comparison_%s.fig', timestamp));
saveas(fig, sprintf('batch_comparison_%s.png', timestamp));

fprintf('✓ Comparison plots saved\n');

end
