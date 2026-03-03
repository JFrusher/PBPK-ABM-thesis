% Check MC results for data quality issues
clear; clc;

% Load the most recent results
results_dir = 'D:\PhysiCell\MC_results';
matfiles = dir(fullfile(results_dir, 'MC_5FU_PK_sensitivity_*.mat'));
[~, idx] = sort([matfiles.datenum], 'descend');
latest = fullfile(results_dir, matfiles(idx(1)).name);

fprintf('Loading results from: %s\n', latest);
load(latest, 'results');

fprintf('\n=== DATA QUALITY CHECK ===\n');
fprintf('Total runs: %d\n', length(results.Cmax));
fprintf('Valid Cmax values: %d\n', sum(isfinite(results.Cmax)));
fprintf('NaN Cmax values: %d\n', sum(isnan(results.Cmax)));
fprintf('Inf Cmax values: %d\n', sum(isinf(results.Cmax)));
fprintf('Negative Cmax: %d\n', sum(results.Cmax < 0));

fprintf('\nCmax statistics (finite values only):\n');
valid_cmax = results.Cmax(isfinite(results.Cmax));
if ~isempty(valid_cmax)
    fprintf('  Mean: %.4f µM\n', mean(valid_cmax));
    fprintf('  Std:  %.4f µM\n', std(valid_cmax));
    fprintf('  Min:  %.4f µM\n', min(valid_cmax));
    fprintf('  Max:  %.4f µM\n', max(valid_cmax));
    fprintf('  Median: %.4f µM\n', median(valid_cmax));
else
    fprintf('  NO VALID DATA!\n');
end

fprintf('\nAUC statistics (finite values only):\n');
valid_auc = results.AUC(isfinite(results.AUC));
if ~isempty(valid_auc)
    fprintf('  Mean: %.4f mg·h/L\n', mean(valid_auc));
    fprintf('  Std:  %.4f mg·h/L\n', std(valid_auc));
    fprintf('  Min:  %.4f mg·h/L\n', min(valid_auc));
    fprintf('  Max:  %.4f mg·h/L\n', max(valid_auc));
    fprintf('  Median: %.4f mg·h/L\n', median(valid_auc));
else
    fprintf('  NO VALID DATA!\n');
end

% Check distribution
fprintf('\n=== DISTRIBUTION CHECK ===\n');
if ~isempty(valid_cmax)
    figure('Position', [100 100 1200 500]);
    
    subplot(1,2,1);
    histogram(valid_cmax, 15);
    xlabel('Cmax (µM)');
    ylabel('Frequency');
    title(sprintf('Cmax Distribution (n=%d)', length(valid_cmax)));
    grid on;
    
    subplot(1,2,2);
    boxplot(valid_cmax);
    ylabel('Cmax (µM)');
    title('Cmax Boxplot');
    grid on;
end

% Check first few simulation outputs
fprintf('\n=== FIRST 3 RUNS DETAILED CHECK ===\n');
mc_dir = 'D:\PhysiCell\MC_results';
for i = 1:min(3, length(results.Cmax))
    csv_file = fullfile(mc_dir, sprintf('MC_run_%04d_5FU_compartments.csv', i));
    if isfile(csv_file)
        T = readtable(csv_file);
        fprintf('\nRun %d:\n', i);
        fprintf('  CSV exists: yes\n');
        fprintf('  Rows: %d\n', height(T));
        fprintf('  Variables: %s\n', sprintf('%s ', T.Properties.VariableNames{:}));
        
        % Find central compartment
        central_cols = T.Properties.VariableNames(contains(T.Properties.VariableNames, 'entral', 'IgnoreCase', true));
        if ~isempty(central_cols)
            col = central_cols{1};
            C = T.(col);
            fprintf('  Central column: %s\n', col);
            fprintf('  Mean: %.4f, Max: %.4f, Min: %.4f\n', mean(C), max(C), min(C));
            fprintf('  Result Cmax: %.4f µM\n', results.Cmax(i));
        end
    else
        fprintf('Run %d: CSV NOT FOUND\n', i);
    end
end
