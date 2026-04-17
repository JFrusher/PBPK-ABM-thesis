function [stats_table, report_filepath] = plot_phases_errorbar(parent_folder, interval_min)

% Enhanced version with statistical analysis and report generation
% Usage: [stats_table, filepath] = plot_phases_errorbar('path\to\simulation_results', 30)

% Plots cell phase counts over time with error bars showing mean +/- standard error
% across all simulation runs, then generates comprehensive oncological statistics report

% Inputs:
%   parent_folder   - Path to folder containing sim_* subdirectories
%   interval_min    - Time interval between samples in minutes (default: 60)
%
% Outputs:
%   stats_table     - Table with calculated statistics
%   report_filepath - Path to generated .txt report file

if nargin < 2
    interval_min = 60; % Default if not provided
end

phaseRow = 8; % current_phase from label index (MATLAB is 1-based)
phase_values = [4, 10, 12, 13, 17, 18, 100, 101, 102];
phase_labels = containers.Map(...
    phase_values, ...
    {'G1', 'S', 'G2', 'M', 'Apoptosis', 'Necrosis', 'Dead1', 'Dead2', 'Dead3'});

% Get all simulation folders
sim_folders = dir(fullfile(parent_folder, 'sim_*'));
sim_folders = sim_folders([sim_folders.isdir]);
num_sims = length(sim_folders);

if num_sims == 0
    error('No simulation folders found matching pattern sim_*');
end

fprintf('Found %d simulation folders\n', num_sims);

% Storage for all data across simulations
all_sim_counts = {}; % Cell array: each cell contains counts matrix for one sim
all_sim_times = {}; % Cell array: each cell contains time vector for one sim

% Process each simulation
for sim_idx = 1:num_sims
    sim_path = fullfile(parent_folder, sim_folders(sim_idx).name);
    fprintf('[%2d/%d] Processing: %s\n', sim_idx, num_sims, sim_folders(sim_idx).name);
    
    % Get all *_cells.mat files in this simulation
    files = dir(fullfile(sim_path, '*_cells.mat'));
    sim_counts = [];
    sim_times = [];
    
    for i = 1:length(files)
        fname = files(i).name;
        load(fullfile(sim_path, fname), 'cells');
        phases = round(cells(phaseRow, :)); % Get integer phase vector
        mask = ismember(phases, phase_values);
        phases = phases(mask);
        counts = zeros(size(phase_values));
        
        for k = 1:length(phase_values)
            counts(k) = sum(phases == phase_values(k));
        end
        
        numMatch = regexp(fname,'output(\d+)_','tokens');
        if ~isempty(numMatch)
            t = str2double(numMatch{1}{1});
        else
            t = i-1; % fallback
        end
        
        sim_times(i) = t;
        if i == 1
            sim_counts = zeros(length(files), length(phase_values));
        end
        sim_counts(i, :) = counts;
    end
    
    sim_times = sim_times * interval_min / 60; % Convert to hours
    all_sim_counts{sim_idx} = sim_counts;
    all_sim_times{sim_idx} = sim_times;
end

% Align all time vectors
ref_times = all_sim_times{1};
num_timepoints = length(ref_times);

% Ensure all simulations have same time points
aligned_counts = zeros(num_sims, num_timepoints, length(phase_values));
for sim_idx = 1:num_sims
    if length(all_sim_times{sim_idx}) == num_timepoints
        aligned_counts(sim_idx, :, :) = all_sim_counts{sim_idx};
    else
        warning('Simulation %d has different number of time points (%d vs %d)', ...
            sim_idx, length(all_sim_times{sim_idx}), num_timepoints);
        min_t = min(length(all_sim_times{sim_idx}), num_timepoints);
        aligned_counts(sim_idx, 1:min_t, :) = all_sim_counts{sim_idx}(1:min_t, :);
    end
end

% Calculate statistics
mean_counts = squeeze(mean(aligned_counts, 1)); % (timepoints, phases)
std_counts = squeeze(std(aligned_counts, 0, 1)); % (timepoints, phases)
se_counts = std_counts / sqrt(num_sims); % Standard error
min_counts = squeeze(min(aligned_counts, 1));
max_counts = squeeze(max(aligned_counts, 1));

% Plot with error bars
figure('Position', [100 100 1400 700]);
hold on;

% Define colors for each phase
colors = [
    0.2 0.6 0.8; % G1 - blue
    1.0 0.5 0.0; % S - orange
    0.6 0.2 0.8; % G2 - purple
    1.0 0.0 0.0; % M - red
    0.0 0.7 0.0; % Apoptosis - green
    0.7 0.0 0.7; % Necrosis - magenta
    0.0 0.0 0.0; % Dead1 - black
    0.5 0.5 0.5; % Dead2 - grey
    1.0 1.0 0.0; % Dead3 - yellow
];

% Plot error bars for each phase
for phase_idx = 1:length(phase_values)
    phase_val = phase_values(phase_idx);
    phase_name = phase_labels(phase_val);
    errorbar(ref_times, mean_counts(:, phase_idx), se_counts(:, phase_idx), ...
        'o-', 'Color', colors(phase_idx, :), 'LineWidth', 2.0, ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(phase_idx, :), ...
        'MarkerEdgeColor', colors(phase_idx, :), 'CapSize', 5, ...
        'DisplayName', phase_name);
end

hold off;
xlabel('Time (hours)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cell Count (Mean ± SEM)', 'FontSize', 12, 'FontWeight', 'bold');

[~, folder_name, ~] = fileparts(parent_folder);
title(sprintf('Cell Phases Over Time - %d Simulations | %s', num_sims, folder_name), ...
    'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
legend('Location', 'Best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);

% =====================================================================
% GENERATE COMPREHENSIVE STATISTICAL REPORT FOR ONCOLOGY
% =====================================================================

% Create report filename with timestamp
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
report_filename = sprintf('CellPhase_DosageAnalysis_%s.txt', timestamp);
report_filepath = fullfile(parent_folder, report_filename);

% Open file for writing
fid = fopen(report_filepath, 'w');

% Write report header
fprintf(fid, '================================================================================\n');
fprintf(fid, '         TUMOUR CELL PHASE ANALYSIS REPORT FOR DOSING STRATEGY EVALUATION\n');
fprintf(fid, '================================================================================\n\n');
fprintf(fid, 'Report Generated: %s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
fprintf(fid, 'Data Source: %s\n', folder_name);
fprintf(fid, 'Simulation Context: %s\n\n', parent_folder);

% Analysis metadata
fprintf(fid, '--- ANALYSIS METADATA ---\n');
fprintf(fid, 'Number of Independent Simulations: %d\n', num_sims);
fprintf(fid, 'Time Points Sampled: %d\n', num_timepoints);
fprintf(fid, 'Sampling Interval: %d minutes (%g hours)\n', interval_min, interval_min/60);
fprintf(fid, 'Total Time Coverage: %.1f hours (%.1f days)\n', ref_times(end), ref_times(end)/24);
fprintf(fid, 'Time Points Range: %.1f - %.1f hours\n\n', ref_times(1), ref_times(end));

fprintf(fid, 'Phases Tracked (9 total):\n');
for phase_idx = 1:length(phase_values)
    phase_val = phase_values(phase_idx);
    phase_name = phase_labels(phase_val);
    fprintf(fid, '  [%d] %s (ID: %d)\n', phase_idx, phase_name, phase_val);
end
fprintf(fid, '\n');

% Calculate and report global statistics
fprintf(fid, '================================================================================\n');
fprintf(fid, 'GLOBAL SUMMARY STATISTICS\n');
fprintf(fid, '================================================================================\n\n');

total_cells_all = zeros(num_sims, num_timepoints);
for sim_idx = 1:num_sims
    for tp = 1:num_timepoints
        total_cells_all(sim_idx, tp) = sum(aligned_counts(sim_idx, tp, :));
    end
end

mean_total_cells = mean(total_cells_all, 1);
std_total_cells = std(total_cells_all, 0, 1);
se_total_cells = std_total_cells / sqrt(num_sims);

fprintf(fid, 'TOTAL CELL POPULATION DYNAMICS:\n');
fprintf(fid, '  Initial Population (t=0):        %.1f ± %.1f cells\n', ...
    mean_total_cells(1), se_total_cells(1));
fprintf(fid, '  Final Population (t=%.1f h):      %.1f ± %.1f cells\n', ...
    ref_times(end), mean_total_cells(end), se_total_cells(end));

percent_change = ((mean_total_cells(end) - mean_total_cells(1)) / mean_total_cells(1)) * 100;
fprintf(fid, '  Net Change:                      %+.1f%% [%+.1f cells]\n\n', ...
    percent_change, mean_total_cells(end) - mean_total_cells(1));

% Find minimum and maximum total population across all time points
[min_pop, min_idx] = min(mean_total_cells);
[max_pop, max_idx] = max(mean_total_cells);
fprintf(fid, '  Minimum Population:              %.1f cells at t = %.1f hours\n', min_pop, ref_times(min_idx));
fprintf(fid, '  Maximum Population:              %.1f cells at t = %.1f hours\n', max_pop, ref_times(max_idx));
fprintf(fid, '  Population Variability (CV):     %.2f%%\n\n', (std(mean_total_cells)/mean(mean_total_cells))*100);

% Report statistics for each cell phase
fprintf(fid, '================================================================================\n');
fprintf(fid, 'PHASE-SPECIFIC STATISTICS AND ONCOLOGICAL METRICS\n');
fprintf(fid, '================================================================================\n\n');

% Initialise stats table for output
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
    phase_name = phase_labels(phase_val);
    
    % Extract statistics for this phase across all timepoints
    phase_means = mean_counts(:, phase_idx);
    phase_stds = std_counts(:, phase_idx);
    phase_ses = se_counts(:, phase_idx);
    phase_mins = min_counts(:, phase_idx);
    phase_maxs = max_counts(:, phase_idx);
    
    fprintf(fid, '--- %s (Cell Phase ID: %d) ---\n', phase_name, phase_val);
    
    % Kinetic parameters
    init_count = phase_means(1);
    final_count = phase_means(end);
    [peak_count, peak_idx] = max(phase_means);
    peak_time = ref_times(peak_idx);
    
    fprintf(fid, '\nPopulation Dynamics:\n');
    fprintf(fid, '  Initial Count (t=0):             %.1f ± %.1f cells\n', init_count, phase_ses(1));
    fprintf(fid, '  Final Count (t=%.1f h):          %.1f ± %.1f cells\n', ...
        ref_times(end), final_count, phase_ses(end));
    fprintf(fid, '  Peak Count:                      %.1f cells at t = %.1f hours\n', peak_count, peak_time);
    fprintf(fid, '  Range (min - max):               %.1f - %.1f cells\n', min(phase_mins), max(phase_maxs));
    
    % Calculate growth rate (exponential fit on non-zero initial populations)
    if init_count > 0 && final_count > 0 && final_count ~= init_count
        growth_rate = log(final_count / init_count) / (ref_times(end) - ref_times(1));
        doubling_time = log(2) / abs(growth_rate);
        fprintf(fid, '  Net Growth Rate:                 %.4f per hour\n', growth_rate);
        if growth_rate > 0
            fprintf(fid, '  Doubling Time:                   %.2f hours\n', doubling_time);
        elseif growth_rate < 0
            fprintf(fid, '  Half-life (decline):             %.2f hours\n', log(2)/abs(growth_rate));
        end
    else
        growth_rate = 0;
        doubling_time = NaN;
        fprintf(fid, '  Net Growth Rate:                 Not calculable (zero initial or final count)\n');
    end
    
    % Calculate area under curve (AUC) - integral of cell count over time
    auc = trapz(ref_times, phase_means);
    fprintf(fid, '\nAccumulation Metrics:\n');
    fprintf(fid, '  Area Under Curve (AUC):          %.1f cell-hours\n', auc);
    fprintf(fid, '  Mean Count Across Time:          %.1f ± %.1f cells\n', mean(phase_means), std(phase_means));
    fprintf(fid, '  Coefficient of Variation:        %.2f%%\n', (std(phase_means)/mean(phase_means))*100);
    
    % Calculate phase-specific transition metrics (for dosing evaluation)
    percent_init = (init_count / mean_total_cells(1)) * 100;
    percent_final = (final_count / mean_total_cells(end)) * 100;
    
    fprintf(fid, '\nDosage Response Metrics:\n');
    fprintf(fid, '  Proportion of Initial Population:  %.1f%%\n', percent_init);
    fprintf(fid, '  Proportion of Final Population:    %.1f%%\n', percent_final);
    fprintf(fid, '  Change in Relative Proportion:     %+.1f%%\n', percent_final - percent_init);
    
    % Store in output table
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

% Critical Dosing Evaluation Metrics
fprintf(fid, '================================================================================\n');
fprintf(fid, 'ONCOLOGICAL DOSING STRATEGY EVALUATION METRICS\n');
fprintf(fid, '================================================================================\n\n');

% Calculate cycle-specific metrics (useful for multi-cycle chemotherapy)
fprintf(fid, 'CELL CYCLE KINETICS (G1 + S + G2 + M):\n');
g1_idx = find(ismember(phase_values, 4));
s_idx = find(ismember(phase_values, 10));
g2_idx = find(ismember(phase_values, 12));
m_idx = find(ismember(phase_values, 13));

cycling_cells = mean_counts(:, g1_idx) + mean_counts(:, s_idx) + ...
                mean_counts(:, g2_idx) + mean_counts(:, m_idx);
apoptosis_idx = find(ismember(phase_values, 17));
apoptotic_cells = mean_counts(:, apoptosis_idx);
necrotic_idx = find(ismember(phase_values, 18));
necrotic_cells = mean_counts(:, necrotic_idx);
dead_cells = mean_counts(:, 7) + mean_counts(:, 8) + mean_counts(:, 9); % Dead1 + Dead2 + Dead3

fprintf(fid, '  Initial Cycling Cells:           %.1f cells (%.1f%% of population)\n', ...
    cycling_cells(1), (cycling_cells(1)/mean_total_cells(1))*100);
fprintf(fid, '  Final Cycling Cells:             %.1f cells (%.1f%% of population)\n', ...
    cycling_cells(end), (cycling_cells(end)/mean_total_cells(end))*100);
fprintf(fid, '  Mean Cycling Cells:              %.1f ± %.1f cells\n', ...
    mean(cycling_cells), std(cycling_cells));

fprintf(fid, '\nAPOPTOSIS & NECROSIS (Drug Response Indicators):\n');
fprintf(fid, '  Total Apoptotic Events:          %.1f cells (AUC: %.1f cell-hours)\n', ...
    apoptotic_cells(end), trapz(ref_times, apoptotic_cells));
fprintf(fid, '  Total Necrotic Events:           %.1f cells (AUC: %.1f cell-hours)\n', ...
    necrotic_cells(end), trapz(ref_times, necrotic_cells));
fprintf(fid, '  Total Dead Cells:                %.1f cells\n', dead_cells(end));
fprintf(fid, '  Cell Death Rate (apoptosis+necrosis): %.1f cells/hour\n', ...
    (apoptotic_cells(end) + necrotic_cells(end)) / ref_times(end));

fprintf(fid, '\nTHERAPEUTIC EFFICACY INDICATORS:\n');
total_death = apoptotic_cells(end) + necrotic_cells(end);
efficacy_index = (total_death / mean_total_cells(1)) * 100;
fprintf(fid, '  Total Cell Elimination:          %.1f cells\n', total_death);
fprintf(fid, '  Efficacy Index (% of initial pop): %.1f%%\n', efficacy_index);

if mean_total_cells(end) < mean_total_cells(1)
    fprintf(fid, '  Tumour Control Status:           RESPONSIVE (net reduction)\n');
else
    fprintf(fid, '  Tumour Control Status:           PROGRESSIVE (net growth)\n');
end

fprintf(fid, '\n');

% Sensitivity analysis: identify critical phases
fprintf(fid, '================================================================================\n');
fprintf(fid, 'PHASE SENSITIVITY AND TREATMENT RESPONSE RANKING\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'Phases ranked by AUC (cumulative cell-time exposure):\n\n');

% Create ranking table
phase_rankings = [];
for phase_idx = 1:length(phase_values)
    phase_rankings = [phase_rankings; ...
        stats_table.AUC_cells_h(phase_idx), phase_idx, ...
        stats_table.Peak_Count(phase_idx), ...
        stats_table.Final_Count(phase_idx)];
end

[~, rank_idx] = sort(phase_rankings(:,1), 'descend');

fprintf(fid, 'Rank | Phase    | AUC (cell-h) | Peak Count | Final Count | Dominance (%%)\n');
fprintf(fid, '%-4d | %-8s | %12.1f | %10.1f | %11.1f | %11.1f\n', ...
    1, stats_table.Phase{rank_idx(1)}, ...
    stats_table.AUC_cells_h(rank_idx(1)), ...
    stats_table.Peak_Count(rank_idx(1)), ...
    stats_table.Final_Count(rank_idx(1)), ...
    (stats_table.AUC_cells_h(rank_idx(1)) / sum(stats_table.AUC_cells_h(1:6)))*100);

for r = 2:min(6, length(rank_idx))
    phase_idx = rank_idx(r);
    fprintf(fid, '%-4d | %-8s | %12.1f | %10.1f | %11.1f | %11.1f\n', ...
        r, stats_table.Phase{phase_idx}, ...
        stats_table.AUC_cells_h(phase_idx), ...
        stats_table.Peak_Count(phase_idx), ...
        stats_table.Final_Count(phase_idx), ...
        (stats_table.AUC_cells_h(phase_idx) / sum(stats_table.AUC_cells_h(1:6)))*100);
end

fprintf(fid, '\n(Note: Ranking indicates which cell cycle phases dominate population dynamics\n');
fprintf(fid, ' throughout the treatment period. Higher AUC = greater influence on tumour\n');
fprintf(fid, ' response to dosing strategy.)\n\n');

% Conclusions and recommendations
fprintf(fid, '================================================================================\n');
fprintf(fid, 'CONCLUSIONS AND DOSING STRATEGY RECOMMENDATIONS\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'SUMMARY OF KEY FINDINGS:\n\n');

fprintf(fid, '1. TUMOUR CONTROL:\n');
if mean_total_cells(end) < mean_total_cells(1) * 0.5
    fprintf(fid, '   - Excellent tumour control achieved (>50%% reduction in cell population)\n');
elseif mean_total_cells(end) < mean_total_cells(1)
    fprintf(fid, '   - Moderate tumour control achieved (net reduction in cell population)\n');
else
    fprintf(fid, '   - Limited tumour control (cell population growth observed)\n');
end

fprintf(fid, '\n2. CELL CYCLE DISRUPTION:\n');
if cycling_cells(end) < cycling_cells(1) * 0.3
    fprintf(fid, '   - Strong cell cycle arrest achieved\n');
elseif cycling_cells(end) < cycling_cells(1) * 0.7
    fprintf(fid, '   - Partial cell cycle arrest observed\n');
else
    fprintf(fid, '   - Minimal cell cycle disruption\n');
end

fprintf(fid, '\n3. CELL DEATH INDUCTION:\n');
if total_death > mean_total_cells(1) * 0.3
    fprintf(fid, '   - Strong apoptotic/necrotic response (>30%% cell death)\n');
elseif total_death > mean_total_cells(1) * 0.1
    fprintf(fid, '   - Moderate cell death induction (10-30%%)\n');
else
    fprintf(fid, '   - Limited cell death response\n');
end

fprintf(fid, '\n4. DOSING RECOMMENDATIONS:\n');
fprintf(fid, '   - Current dosing shows %s response profile\n', ...
    iif(mean_total_cells(end) < mean_total_cells(1), 'INHIBITORY', 'STIMULATORY'));
fprintf(fid, '   - Most responsive phases: %s, %s\n', ...
    stats_table.Phase{rank_idx(1)}, stats_table.Phase{rank_idx(2)});
fprintf(fid, '   - Focus on targeting cells in %s phase(s) for optimisation\n', ...
    stats_table.Phase{rank_idx(1)});
fprintf(fid, '   - Consider timing of dose delivery to maximise S-phase synchronisation\n');

fprintf(fid, '\n');
fprintf(fid, '================================================================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '================================================================================\n');

fclose(fid);

% Display file location
fprintf('\n*** STATISTICAL ANALYSIS REPORT GENERATED ***\n');
fprintf('Report saved to: %s\n', report_filepath);
fprintf('File size: %.1f KB\n\n', dir(report_filepath).bytes / 1024);

% Also display a condensed version to console
fprintf('=== QUICK SUMMARY ===\n');
fprintf('Initial population:    %.1f cells\n', mean_total_cells(1));
fprintf('Final population:      %.1f cells\n', mean_total_cells(end));
fprintf('Total cell deaths:     %.1f (apoptosis + necrosis)\n', total_death);
fprintf('Efficacy index:        %.1f%% of initial population eliminated\n', efficacy_index);

end

% Helper function for conditional string
function result = iif(condition, true_str, false_str)
    if condition
        result = true_str;
    else
        result = false_str;
    end
end
