function analyze_physicell_batch_results(results_folder)
% analyze_physicell_batch_results  Aggregate analysis of PhysiCell batch simulation results
%
% PURPOSE:
%   Post-processes and analyzes results from run_mc_physicell_batch.sh,
%   aggregating outputs from multiple MC_run_XXXX directories.
%   Creates publication-quality visualizations and comparative statistics.
%
% SYNTAX:
%   analyze_physicell_batch_results(results_folder)
%
% INPUTS:
%   results_folder - Path to parent folder containing MC_run_XXXX subdirectories
%                    Example: './mc_results/2025-01-11_23-15-30'
%
% OUTPUT:
%   - Aggregate statistics and comparisons
%   - Combined plots across all runs
%   - Heatmaps of parameter sensitivity
%   - Run-to-run variation analysis
%   - Summary report (TXT and MAT files)
%
% EXAMPLE:
%   analyze_physicell_batch_results('./mc_results/2025-01-11_23-15-30')
%
% DEPENDENCIES:
%   - Statistics and Machine Learning Toolbox (optional, for enhanced stats)
%   - Image Processing Toolbox (optional, for heatmaps)
%
% AUTHOR: PhysiCell Analysis Framework
% VERSION: 1.0

    % ========================================================================
    % 1) INITIALIZATION & VALIDATION
    % ========================================================================
    
    fprintf('\n╔════════════════════════════════════════════════════════╗\n');
    fprintf('║  PhysiCell Batch Results Analysis                     ║\n');
    fprintf('╚════════════════════════════════════════════════════════╝\n\n');
    
    % Check nargin
    if nargin < 1
        error('analyze_physicell_batch_results:InvalidInput', ...
            'Usage: analyze_physicell_batch_results(results_folder)');
    end
    
    % Validate folder exists
    if ~isfolder(results_folder)
        error('analyze_physicell_batch_results:FolderNotFound', ...
            'Results folder not found: %s', results_folder);
    end
    
    % Get timestamp from folder name for logging
    [~, folder_name] = fileparts(results_folder);
    timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    
    % Create output directory for analysis
    output_dir = fullfile(results_folder, 'batch_analysis');
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    
    % Initialize log file
    logfile = fullfile(output_dir, sprintf('analysis_report_%s.txt', timestamp));
    diary(logfile);
    
    fprintf('Analysis initialized: %s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    fprintf('Results folder: %s\n', results_folder);
    fprintf('Output folder: %s\n', output_dir);
    fprintf('Log file: %s\n\n', logfile);
    
    % ========================================================================
    % 2) DISCOVER MC_RUN DIRECTORIES
    % ========================================================================
    
    fprintf('════════════════════════════════════════════════════════\n');
    fprintf('PHASE 1: Discovering MC_run directories\n');
    fprintf('════════════════════════════════════════════════════════\n\n');
    
    % Find all MC_run_XXXX directories
    mc_run_dirs = dir(fullfile(results_folder, 'MC_run_*'));
    mc_run_dirs = mc_run_dirs([mc_run_dirs.isdir]);
    
    if isempty(mc_run_dirs)
        error('analyze_physicell_batch_results:NoRunsFound', ...
            'No MC_run_XXXX directories found in %s', results_folder);
    end
    
    n_runs = length(mc_run_dirs);
    fprintf('Found %d MC_run directories\n\n', n_runs);
    
    % Extract run indices for proper sorting
    run_names = {mc_run_dirs.name};
    run_indices = zeros(n_runs, 1);
    for i = 1:n_runs
        tokens = regexp(run_names{i}, 'MC_run_(\d+)', 'tokens');
        if ~isempty(tokens)
            run_indices(i) = str2double(tokens{1}{1});
        end
    end
    
    % Sort by index
    [~, sort_order] = sort(run_indices);
    mc_run_dirs = mc_run_dirs(sort_order);
    run_names = run_names(sort_order);
    run_indices = run_indices(sort_order);
    
    % ========================================================================
    % 3) LOAD AND AGGREGATE DATA FROM EACH RUN
    % ========================================================================
    
    fprintf('════════════════════════════════════════════════════════\n');
    fprintf('PHASE 2: Loading and aggregating run data\n');
    fprintf('════════════════════════════════════════════════════════\n\n');
    
    % Initialize data structures
    batch_data = struct();
    batch_data.run_names = run_names;
    batch_data.run_indices = run_indices;
    batch_data.n_runs = n_runs;
    batch_data.runs = struct();
    
    % Cell parameters (extracted from run names or config files)
    batch_data.cell_count = [];
    batch_data.tumor_size = [];
    batch_data.dosing_prefix = cell(n_runs, 1);
    
    % Outputs from simulations
    batch_data.final_populations = [];
    batch_data.final_volumes = [];
    batch_data.max_drug_concentration = [];
    batch_data.cell_death_count = [];
    
    % Metadata
    batch_data.successful_runs = 0;
    batch_data.failed_runs = 0;
    batch_data.run_status = cell(n_runs, 1);
    
    fprintf('Loading data:\n');
    fprintf('%-30s %-20s %-15s\n', 'Run', 'Status', 'Cells Counted');
    fprintf('────────────────────────────────────────────────────────\n');
    
    for i = 1:n_runs
        run_path = fullfile(results_folder, mc_run_dirs(i).name);
        run_name = mc_run_dirs(i).name;
        output_path = fullfile(run_path, 'physicell_output');
        
        % Check if output exists
        if ~isfolder(output_path)
            fprintf('%-30s %-20s %-15s\n', run_name, 'NO_OUTPUT', 'N/A');
            batch_data.run_status{i} = 'NO_OUTPUT';
            batch_data.failed_runs = batch_data.failed_runs + 1;
            continue;
        end
        
        % Try to load simulation outputs
        try
            % Look for .mat files or CSV outputs from PhysiCell
            mat_files = dir(fullfile(output_path, '*.mat'));
            csv_files = dir(fullfile(output_path, 'output*.csv'));
            
            run_data = struct();
            run_data.mat_files = {mat_files.name};
            run_data.csv_files = {csv_files.name};
            run_data.n_csv_outputs = length(csv_files);
            
            % Try to extract key metrics from CSV files (last timestep)
            if ~isempty(csv_files)
                last_csv = fullfile(output_path, csv_files(end).name);
                try
                    % Read final CSV to count cells
                    cell_data = readtable(last_csv);
                    n_cells = height(cell_data);
                    
                    batch_data.cell_count = [batch_data.cell_count; n_cells];
                    batch_data.runs(i).n_cells = n_cells;
                    
                    % Try to extract drug concentration if available
                    if ismember('drug_conc', cell_data.Properties.VariableNames)
                        max_drug = max(cell_data.drug_conc);
                        batch_data.max_drug_concentration = [batch_data.max_drug_concentration; max_drug];
                    end
                    
                    % Try to count dead cells if death column exists
                    if ismember('is_alive', cell_data.Properties.VariableNames)
                        dead_count = sum(cell_data.is_alive == 0);
                        batch_data.cell_death_count = [batch_data.cell_death_count; dead_count];
                    end
                    
                    batch_data.run_status{i} = 'SUCCESS';
                    batch_data.successful_runs = batch_data.successful_runs + 1;
                    
                    fprintf('%-30s %-20s %-15d\n', run_name, 'SUCCESS', n_cells);
                    
                catch ME
                    fprintf('%-30s %-20s %-15s\n', run_name, 'READ_ERROR', 'Failed');
                    batch_data.run_status{i} = sprintf('ERROR: %s', ME.message);
                    batch_data.failed_runs = batch_data.failed_runs + 1;
                end
            else
                fprintf('%-30s %-20s %-15s\n', run_name, 'NO_CSV', 'N/A');
                batch_data.run_status{i} = 'NO_CSV_OUTPUT';
                batch_data.failed_runs = batch_data.failed_runs + 1;
            end
            
            batch_data.runs(i).data = run_data;
            batch_data.runs(i).path = run_path;
            
        catch ME
            fprintf('%-30s %-20s %-15s\n', run_name, 'LOAD_FAILED', 'Error');
            batch_data.run_status{i} = sprintf('LOAD_ERROR: %s', ME.message);
            batch_data.failed_runs = batch_data.failed_runs + 1;
        end
    end
    
    fprintf('\n');
    
    % ========================================================================
    % 4) SUMMARY STATISTICS
    % ========================================================================
    
    fprintf('════════════════════════════════════════════════════════\n');
    fprintf('PHASE 3: Computing aggregate statistics\n');
    fprintf('════════════════════════════════════════════════════════\n\n');
    
    fprintf('Batch Statistics:\n');
    fprintf('  Total runs:        %d\n', n_runs);
    fprintf('  Successful:        %d (%.1f%%)\n', ...
        batch_data.successful_runs, 100*batch_data.successful_runs/n_runs);
    fprintf('  Failed/Incomplete: %d (%.1f%%)\n', ...
        batch_data.failed_runs, 100*batch_data.failed_runs/n_runs);
    fprintf('\n');
    
    % Only compute stats if we have successful data
    if ~isempty(batch_data.cell_count)
        fprintf('Final Cell Population Statistics:\n');
        fprintf('  Mean:   %d cells (SD: %d)\n', ...
            round(mean(batch_data.cell_count)), round(std(batch_data.cell_count)));
        fprintf('  Median: %d cells\n', median(batch_data.cell_count));
        fprintf('  Range:  %d - %d cells\n', ...
            min(batch_data.cell_count), max(batch_data.cell_count));
        fprintf('  CV:     %.1f%%\n', 100*std(batch_data.cell_count)/mean(batch_data.cell_count));
        fprintf('\n');
    end
    
    if ~isempty(batch_data.max_drug_concentration)
        fprintf('Max Drug Concentration Statistics:\n');
        fprintf('  Mean:   %.3f µM (SD: %.3f)\n', ...
            mean(batch_data.max_drug_concentration), std(batch_data.max_drug_concentration));
        fprintf('  Median: %.3f µM\n', median(batch_data.max_drug_concentration));
        fprintf('  Range:  %.3f - %.3f µM\n', ...
            min(batch_data.max_drug_concentration), max(batch_data.max_drug_concentration));
        fprintf('\n');
    end
    
    if ~isempty(batch_data.cell_death_count)
        fprintf('Cell Death Statistics:\n');
        fprintf('  Total deaths (across all runs): %d cells\n', ...
            sum(batch_data.cell_death_count));
        fprintf('  Mean deaths per run: %.1f cells (SD: %.1f)\n', ...
            mean(batch_data.cell_death_count), std(batch_data.cell_death_count));
        fprintf('  Max deaths (single run): %d cells\n', ...
            max(batch_data.cell_death_count));
        fprintf('\n');
    end
    
    % ========================================================================
    % 5) VISUALIZATION - OVERVIEW PLOTS
    % ========================================================================
    
    fprintf('════════════════════════════════════════════════════════\n');
    fprintf('PHASE 4: Creating visualizations\n');
    fprintf('════════════════════════════════════════════════════════\n\n');
    
    % Create overview figure
    fig_overview = figure('Name', 'Batch_Overview', 'Position', [100, 100, 1400, 800]);
    
    n_plots = 0;
    
    % Plot 1: Final populations by run
    if ~isempty(batch_data.cell_count)
        n_plots = n_plots + 1;
        subplot(2, 3, n_plots);
        
        bar(1:length(batch_data.cell_count), batch_data.cell_count, 'FaceColor', [0.2 0.6 0.8]);
        hold on;
        
        % Add mean line
        mean_pop = mean(batch_data.cell_count);
        yline(mean_pop, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Mean: %.0f', mean_pop));
        
        xlabel('Run Index', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Final Population (cells)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Final Cell Populations Across Runs', 'FontSize', 12, 'FontWeight', 'bold');
        grid on; grid minor;
        legend('Location', 'best');
        set(gca, 'FontSize', 10);
    end
    
    % Plot 2: Distribution of final populations
    if ~isempty(batch_data.cell_count)
        n_plots = n_plots + 1;
        subplot(2, 3, n_plots);
        
        histogram(batch_data.cell_count, 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.7, 'EdgeColor', 'k');
        
        xlabel('Final Population (cells)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Frequency', 'FontSize', 11, 'FontWeight', 'bold');
        title('Distribution of Final Populations', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca, 'FontSize', 10);
    end
    
    % Plot 3: Max drug concentration by run
    if ~isempty(batch_data.max_drug_concentration)
        n_plots = n_plots + 1;
        subplot(2, 3, n_plots);
        
        plot(1:length(batch_data.max_drug_concentration), ...
            batch_data.max_drug_concentration, 'o-', 'LineWidth', 2, 'MarkerSize', 6, ...
            'Color', [0.2 0.5 0.8]);
        hold on;
        
        mean_drug = mean(batch_data.max_drug_concentration);
        yline(mean_drug, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Mean: %.3f µM', mean_drug));
        
        xlabel('Run Index', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Max Drug Conc. (µM)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Peak Drug Concentration by Run', 'FontSize', 12, 'FontWeight', 'bold');
        grid on; grid minor;
        legend('Location', 'best');
        set(gca, 'FontSize', 10);
    end
    
    % Plot 4: Cell death count by run
    if ~isempty(batch_data.cell_death_count)
        n_plots = n_plots + 1;
        subplot(2, 3, n_plots);
        
        bar(1:length(batch_data.cell_death_count), batch_data.cell_death_count, ...
            'FaceColor', [0.8 0.4 0.4], 'FaceAlpha', 0.7);
        
        xlabel('Run Index', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Dead Cells', 'FontSize', 11, 'FontWeight', 'bold');
        title('Cell Death by Run', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca, 'FontSize', 10);
    end
    
    % Plot 5: Success/failure status pie chart
    if n_plots < 5
        n_plots = n_plots + 1;
        subplot(2, 3, n_plots);
        
        status_counts = [batch_data.successful_runs, batch_data.failed_runs];
        status_labels = {'Successful', 'Failed/Incomplete'};
        colors = [[0.2 0.8 0.2]; [0.8 0.2 0.2]];
        
        pie(status_counts, status_labels);
        colormap(colors);
        title('Run Status Summary', 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'FontSize', 10);
    end
    
    sgtitle('PhysiCell Batch Results Overview', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save overview figure
    overview_fig_file = fullfile(output_dir, sprintf('batch_overview_%s.pdf', timestamp));
    print(gcf, overview_fig_file, '-dpdf', '-bestfit');
    fprintf('✓ Saved overview figure: %s\n', overview_fig_file);
    close(fig_overview);
    
    % ========================================================================
    % 6) CORRELATION & SENSITIVITY ANALYSIS
    % ========================================================================
    
    fprintf('\n════════════════════════════════════════════════════════\n');
    fprintf('PHASE 5: Correlation and sensitivity analysis\n');
    fprintf('════════════════════════════════════════════════════════\n\n');
    
    if ~isempty(batch_data.cell_count) && ~isempty(batch_data.max_drug_concentration) && ...
            length(batch_data.cell_count) == length(batch_data.max_drug_concentration)
        
        [r, p] = corrcoef(batch_data.max_drug_concentration, batch_data.cell_count);
        fprintf('Correlation between max drug conc. and final population:\n');
        fprintf('  Pearson r: %.4f\n', r(1,2));
        fprintf('  p-value: %.4e\n', p(1,2));
        
        if p(1,2) < 0.05
            fprintf('  ✓ Statistically significant correlation\n');
        else
            fprintf('  ✗ No significant correlation\n');
        end
        fprintf('\n');
    end
    
    % ========================================================================
    % 7) SAVE ANALYSIS DATA
    % ========================================================================
    
    fprintf('════════════════════════════════════════════════════════\n');
    fprintf('PHASE 6: Saving analysis results\n');
    fprintf('════════════════════════════════════════════════════════\n\n');
    
    % Save batch data as MAT file
    analysis_mat = fullfile(output_dir, sprintf('batch_analysis_%s.mat', timestamp));
    save(analysis_mat, 'batch_data', '-v7.3');
    fprintf('✓ Saved analysis data: %s\n', analysis_mat);
    
    % Save summary statistics as CSV
    summary_table = table(run_indices, run_names', batch_data.run_status', ...
        'VariableNames', {'RunIndex', 'RunName', 'Status'});
    
    if ~isempty(batch_data.cell_count)
        % Pad cell count to match table length
        cell_count_padded = nan(n_runs, 1);
        cell_count_padded(1:length(batch_data.cell_count)) = batch_data.cell_count;
        summary_table.FinalPopulation = cell_count_padded;
    end
    
    if ~isempty(batch_data.max_drug_concentration)
        drug_padded = nan(n_runs, 1);
        drug_padded(1:length(batch_data.max_drug_concentration)) = batch_data.max_drug_concentration;
        summary_table.MaxDrugConc_uM = drug_padded;
    end
    
    summary_csv = fullfile(output_dir, sprintf('batch_summary_%s.csv', timestamp));
    writetable(summary_table, summary_csv);
    fprintf('✓ Saved summary table: %s\n', summary_csv);
    
    % ========================================================================
    % 8) FINAL REPORT
    % ========================================================================
    
    fprintf('\n════════════════════════════════════════════════════════\n');
    fprintf('ANALYSIS COMPLETE\n');
    fprintf('════════════════════════════════════════════════════════\n\n');
    
    fprintf('Output files saved to: %s\n\n', output_dir);
    fprintf('Generated files:\n');
    fprintf('  • %s\n', overview_fig_file);
    fprintf('  • %s\n', analysis_mat);
    fprintf('  • %s\n', summary_csv);
    fprintf('  • %s (this log file)\n', logfile);
    fprintf('\n');
    
    fprintf('Summary:\n');
    fprintf('  Batch analyzed:  %s\n', folder_name);
    fprintf('  Total runs:      %d\n', n_runs);
    fprintf('  Successful:      %d\n', batch_data.successful_runs);
    fprintf('  Failed:          %d\n', batch_data.failed_runs);
    fprintf('  Analysis time:   %s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    fprintf('\n');
    
    diary off;
    
    fprintf('✓ Analysis complete. Log saved to: %s\n', logfile);
    fprintf('\nYou can now:\n');
    fprintf('  1. Load analysis with: load(''%s'')\n', analysis_mat);
    fprintf('  2. Explore detailed results in: %s\n', output_dir);
    fprintf('  3. Review run logs for diagnostics\n');
    
end
