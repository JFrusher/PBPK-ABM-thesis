function results = analyze_MC10k_comprehensive(mergedCsvPath, outDir, options)
% analyze_MC10k_comprehensive
% Comprehensive post-merge analysis for chunked MC 10k outputs.
%
% Generates dissertation-level analysis outputs:
%   1) QC and completeness checks
%   2) Distribution plots (hist, ECDF, violin/box fallback)
%   3) Sensitivity tornado plots (Spearman rho)
%   4) Scatter plots and trend lines (top drivers)
%   5) Correlation heatmap
%   6) Target-attainment communication figures and tables
%   7) Task-level throughput/QC plots when chunk_task_id exists
%
% Usage:
%   analyze_MC10k_comprehensive
%   analyze_MC10k_comprehensive('MC_results/MC10k_chunks/MC_10k_merged_outputs.csv')
%   analyze_MC10k_comprehensive('MC_results/MC10k_chunks/MC_10k_merged_outputs.csv', 'MC_results/MC10k_chunks/analysis_pack')
%
% Output files are written into outDir.

if nargin < 1 || isempty(mergedCsvPath)
    mergedCsvPath = fullfile('MC_10k_merged_outputs_cont.csv');
end

if nargin < 2 || isempty(outDir)
    stamp = char(string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    outDir = fullfile(fileparts(mergedCsvPath), ['analysis_pack_' stamp]);
end

if nargin < 3 || isempty(options)
    options = struct();
end

options = apply_default_options(options);

if ~isfile(mergedCsvPath)
    error('Merged CSV not found: %s', mergedCsvPath);
end
if ~isfolder(outDir)
    mkdir(outDir);
end

fprintf('\n=== MC10k Comprehensive Analysis ===\n');
fprintf('Input:  %s\n', mergedCsvPath);
fprintf('Output: %s\n\n', outDir);

T = readtable(mergedCsvPath, 'VariableNamingRule', 'preserve');
if height(T) == 0
    error('Input table is empty: %s', mergedCsvPath);
end

required = {'AUC_mg_h_L', 'Cmax_uM'};
for i = 1:numel(required)
    if ~ismember(required{i}, T.Properties.VariableNames)
        error('Required column missing: %s', required{i});
    end
end

results = struct();
results.input_file = mergedCsvPath;
results.output_dir = outDir;
results.n_rows_raw = height(T);

% Convert key columns to numeric safely.
T.AUC_mg_h_L = to_numeric(T.AUC_mg_h_L);
T.Cmax_uM = to_numeric(T.Cmax_uM);

isGood = isfinite(T.AUC_mg_h_L) & isfinite(T.Cmax_uM);
Tclean = T(isGood, :);
results.n_rows_clean = height(Tclean);
results.n_rows_dropped = height(T) - height(Tclean);

fprintf('Rows raw:   %d\n', results.n_rows_raw);
fprintf('Rows clean: %d\n', results.n_rows_clean);
fprintf('Dropped:    %d\n\n', results.n_rows_dropped);

% Parameter columns = numeric columns excluding outcomes and metadata.
excludeCols = {'AUC_mg_h_L', 'Cmax_uM', 'sample_local_index', 'chunk_task_id'};
numericMask = varfun(@isnumeric, Tclean, 'OutputFormat', 'uniform');
numVars = Tclean.Properties.VariableNames(numericMask);
paramCols = setdiff(numVars, excludeCols, 'stable');
results.parameter_columns = paramCols;

if isempty(paramCols)
    warning('No parameter columns detected; sensitivity plots will be skipped.');
end

% Core descriptive stats and communication tables.
summaryTbl = build_summary_table(Tclean, options);
writetable(summaryTbl, fullfile(outDir, 'summary_statistics.csv'));

targetTbl = build_target_attainment_table(Tclean, options);
writetable(targetTbl, fullfile(outDir, 'target_attainment_summary.csv'));

% Per-task QC if chunk_task_id available.
if ismember('chunk_task_id', Tclean.Properties.VariableNames)
    Tclean.chunk_task_id = to_numeric(Tclean.chunk_task_id);
end

save_qc_report_txt(fullfile(outDir, 'analysis_report.txt'), Tclean, summaryTbl, targetTbl, results, options);

% ---- Plot set ----
figs = {};

figs{end+1} = plot_endpoint_distributions(Tclean, outDir, options); %#ok<AGROW>
figs{end+1} = plot_endpoint_ecdf(Tclean, outDir, options); %#ok<AGROW>
figs{end+1} = plot_violin_or_fallback(Tclean, outDir, options); %#ok<AGROW>
figs{end+1} = plot_auc_vs_cmax_scatter(Tclean, outDir, options); %#ok<AGROW>
figs{end+1} = plot_target_attainment_bars(targetTbl, outDir, options); %#ok<AGROW>

if ~isempty(paramCols)
    [sensAUC, sensCmax] = compute_sensitivity(Tclean, paramCols);
    writetable(sensAUC, fullfile(outDir, 'sensitivity_spearman_auc.csv'));
    writetable(sensCmax, fullfile(outDir, 'sensitivity_spearman_cmax.csv'));

    figs{end+1} = plot_tornado(sensAUC, sensCmax, outDir, options); %#ok<AGROW>
    figs{end+1} = plot_top_scatter_grid(Tclean, sensAUC, sensCmax, outDir, options); %#ok<AGROW>
    figs{end+1} = plot_correlation_heatmap(Tclean, paramCols, outDir, options); %#ok<AGROW>
end

if ismember('chunk_task_id', Tclean.Properties.VariableNames)
    if any(isfinite(Tclean.chunk_task_id))
        figs{end+1} = plot_task_level_qc(Tclean, outDir, options); %#ok<AGROW>
    end
end

results.figure_handles = figs;
save(fullfile(outDir, 'analysis_workspace.mat'), 'results', 'options', 'summaryTbl', 'targetTbl');

fprintf('\nAnalysis complete. Outputs written to:\n  %s\n', outDir);

end

function options = apply_default_options(options)
def = struct();
def.auc_target = [20 30];
def.cmax_target = [300 500];
def.auto_frame_cmax = true;
def.axis_margin_fraction = 0.08;
def.top_n_tornado = 10;
def.top_n_scatter = 6;
def.export_png_dpi = 400;
def.export_pdf = true;
def.pdf_content_type = 'image';
def.use_exportgraphics = false;
def.enable_violin = false;
def.hist_bins = 40;
def.hist_normalization = 'probability';
def.base_font_size = 10;
def.label_font_size = 10;
def.title_font_size = 11;
def.legend_font_size = 9;
def.axis_line_width = 0.9;
def.line_width = 1.6;
def.tick_label_max_chars = 24;
def.max_heatmap_vars = 16;
% Publication sizing policy (A4-aware): width is a fraction of A4 width,
% height is proportional to width (with minimum height guard).
def.a4_width_cm = 21.0;
def.figure_width_fraction_a4 = 0.8;
def.figure_height_to_width_ratio = 0.62;
def.figure_min_height_cm = 7.4;
def.fig_size_cm = []; % auto-computed from publication sizing policy unless explicitly set
def.color_auc = [0.80 0.25 0.25];
def.color_cmax = [0.20 0.45 0.80];
def.color_joint = [0.10 0.10 0.10];
def.alpha = 0.05;

f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(options, f{i}) || isempty(options.(f{i}))
        options.(f{i}) = def.(f{i});
    end
end

% Derive final figure size if not explicitly provided.
if ~isfield(options, 'fig_size_cm') || isempty(options.fig_size_cm)
    w = options.a4_width_cm * options.figure_width_fraction_a4;
    h = max(options.figure_min_height_cm, w * options.figure_height_to_width_ratio);
    options.fig_size_cm = [w h];
end
end

function x = to_numeric(v)
if isnumeric(v)
    x = double(v);
elseif islogical(v)
    x = double(v);
elseif isstring(v) || iscellstr(v) || iscategorical(v)
    x = str2double(string(v));
elseif iscell(v)
    x = str2double(string(v));
else
    try
        x = double(v);
    catch
        x = nan(size(v));
    end
end
x = double(x);
end

function summaryTbl = build_summary_table(T, options)
metrics = {'AUC_mg_h_L', 'Cmax_uM'};
rows = [];
for i = 1:numel(metrics)
    x = to_numeric(T.(metrics{i}));
    x = x(isfinite(x));
    rows = [rows; make_summary_row(metrics{i}, x)]; %#ok<AGROW>
end

summaryTbl = struct2table(rows);
summaryTbl.Properties.VariableNames = {'metric','n','mean','sd','median','q25','q75','p05','p95','min','max'};

% Add concise communication columns.
summaryTbl.cv_percent = 100 * summaryTbl.sd ./ max(summaryTbl.mean, eps);
summaryTbl.target_range = strings(height(summaryTbl), 1);
summaryTbl.target_in_range_percent = nan(height(summaryTbl), 1);
for i = 1:height(summaryTbl)
    switch string(summaryTbl.metric(i))
        case "AUC_mg_h_L"
            lo = options.auc_target(1); hi = options.auc_target(2);
            x = to_numeric(T.AUC_mg_h_L);
        otherwise
            lo = options.cmax_target(1); hi = options.cmax_target(2);
            x = to_numeric(T.Cmax_uM);
    end
    x = x(isfinite(x));
    summaryTbl.target_range(i) = sprintf('[%.3g, %.3g]', lo, hi);
    summaryTbl.target_in_range_percent(i) = 100 * mean(x >= lo & x <= hi);
end
end

function row = make_summary_row(name, x)
row = struct();
row.metric = string(name);
row.n = numel(x);
row.mean = mean(x);
row.sd = std(x);
row.median = median(x);
row.q25 = prctile(x, 25);
row.q75 = prctile(x, 75);
row.p05 = prctile(x, 5);
row.p95 = prctile(x, 95);
row.min = min(x);
row.max = max(x);
end

function targetTbl = build_target_attainment_table(T, options)
auc = to_numeric(T.AUC_mg_h_L);
cmax = to_numeric(T.Cmax_uM);
valid = isfinite(auc) & isfinite(cmax);
auc = auc(valid);
cmax = cmax(valid);

aucIn = auc >= options.auc_target(1) & auc <= options.auc_target(2);
cmaxIn = cmax >= options.cmax_target(1) & cmax <= options.cmax_target(2);
bothIn = aucIn & cmaxIn;

targetTbl = table();
targetTbl.metric = ["AUC"; "Cmax"; "Both"];
targetTbl.target_range = string([
    sprintf('[%.3g, %.3g]', options.auc_target(1), options.auc_target(2)); ...
    sprintf('[%.3g, %.3g]', options.cmax_target(1), options.cmax_target(2)); ...
    "AUC in range AND Cmax in range" ...
    ]);
targetTbl.percent_attained = [100*mean(aucIn); 100*mean(cmaxIn); 100*mean(bothIn)];
targetTbl.n = repmat(numel(auc), 3, 1);
end

function [sensAUC, sensCmax] = compute_sensitivity(T, paramCols)
n = numel(paramCols);
rhoA = nan(n,1); pA = nan(n,1);
rhoC = nan(n,1); pC = nan(n,1);
auc = to_numeric(T.AUC_mg_h_L);
cmax = to_numeric(T.Cmax_uM);

for i = 1:n
    x = to_numeric(T.(paramCols{i}));
    [rhoA(i), pA(i)] = corr(x, auc, 'Type', 'Spearman', 'Rows', 'complete');
    [rhoC(i), pC(i)] = corr(x, cmax, 'Type', 'Spearman', 'Rows', 'complete');
end

sensAUC = table(string(paramCols(:)), rhoA, abs(rhoA), pA, ...
    'VariableNames', {'parameter','rho','abs_rho','p_value'});
sensCmax = table(string(paramCols(:)), rhoC, abs(rhoC), pC, ...
    'VariableNames', {'parameter','rho','abs_rho','p_value'});

sensAUC = sortrows(sensAUC, 'abs_rho', 'descend');
sensCmax = sortrows(sensCmax, 'abs_rho', 'descend');
end

function h = plot_endpoint_distributions(T, outDir, options)
h = new_fig('Endpoint Distributions', options, [1.05 1.0]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
histogram(T.AUC_mg_h_L, options.hist_bins, 'FaceColor', options.color_auc, 'EdgeColor', 'none', ...
    'Normalization', options.hist_normalization);
xline(options.auc_target(1), '--k', 'AUC target min');
xline(options.auc_target(2), '--k', 'AUC target max');
title('AUC Distribution'); xlabel('AUC (mg h/L)'); ylabel('Probability'); grid on;

nexttile;
histogram(T.Cmax_uM, options.hist_bins, 'FaceColor', options.color_cmax, 'EdgeColor', 'none', ...
    'Normalization', options.hist_normalization);
xline(options.cmax_target(1), '--k', 'Cmax target min');
xline(options.cmax_target(2), '--k', 'Cmax target max');
if options.auto_frame_cmax
    xlim(compute_auto_limits(T.Cmax_uM, options.axis_margin_fraction, true));
end
title('Cmax Distribution'); xlabel('Cmax (uM)'); ylabel('Probability'); grid on;

save_fig(h, outDir, '01_endpoint_distributions', options);
end

function h = plot_endpoint_ecdf(T, outDir, options)
h = new_fig('Endpoint ECDF', options, [1.05 1.0]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
ecdf(T.AUC_mg_h_L); hold on;
xline(options.auc_target(1), '--k'); xline(options.auc_target(2), '--k');
title('AUC CDF'); xlabel('AUC (mg h/L)'); ylabel('F(x)'); grid on;

nexttile;
ecdf(T.Cmax_uM); hold on;
xline(options.cmax_target(1), '--k'); xline(options.cmax_target(2), '--k');
if options.auto_frame_cmax
    xlim(compute_auto_limits(T.Cmax_uM, options.axis_margin_fraction, true));
end
title('Cmax CDF'); xlabel('Cmax (uM)'); ylabel('F(x)'); grid on;

save_fig(h, outDir, '02_endpoint_ecdf', options);
end

function h = plot_violin_or_fallback(T, outDir, options)
h = new_fig('Violin Plots', options, [0.95 1.0]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

hasViolin = isfield(options, 'enable_violin') && options.enable_violin && exist('violinplot', 'file') == 2;

nexttile;
if hasViolin
    ok = try_violinplot(T.AUC_mg_h_L, options.color_auc);
    if ~ok
        boxchart(ones(height(T),1), T.AUC_mg_h_L, 'BoxFaceColor', options.color_auc);
    end
else
    boxchart(ones(height(T),1), T.AUC_mg_h_L, 'BoxFaceColor', options.color_auc);
end
title('AUC'); ylabel('AUC (mg h/L)'); xticks([]); grid on;

nexttile;
if hasViolin
    ok = try_violinplot(T.Cmax_uM, options.color_cmax);
    if ~ok
        boxchart(ones(height(T),1), T.Cmax_uM, 'BoxFaceColor', options.color_cmax);
    end
else
    boxchart(ones(height(T),1), T.Cmax_uM, 'BoxFaceColor', options.color_cmax);
end
if options.auto_frame_cmax
    ylim(compute_auto_limits(T.Cmax_uM, options.axis_margin_fraction, true));
end
title('Cmax'); ylabel('Cmax (uM)'); xticks([]); grid on;

save_fig(h, outDir, '03_violin_or_box', options);
end

function h = plot_auc_vs_cmax_scatter(T, outDir, options)
h = new_fig('AUC vs Cmax', options, [1.0 1.0]);
scatter(T.AUC_mg_h_L, T.Cmax_uM, 8, 'MarkerFaceColor', options.color_joint, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.25);
hold on;

xline(options.auc_target(1), '--', 'Color', [0.4 0.4 0.4]);
xline(options.auc_target(2), '--', 'Color', [0.4 0.4 0.4]);
yline(options.cmax_target(1), '--', 'Color', [0.4 0.4 0.4]);
yline(options.cmax_target(2), '--', 'Color', [0.4 0.4 0.4]);
if options.auto_frame_cmax
    ylim(compute_auto_limits(T.Cmax_uM, options.axis_margin_fraction, true));
end

xlabel('AUC (mg h/L)');
ylabel('Cmax (uM)');
title('Joint Exposure Space (AUC vs Cmax)');
grid on;

save_fig(h, outDir, '04_auc_vs_cmax_scatter', options);
end

function h = plot_target_attainment_bars(targetTbl, outDir, options)
h = new_fig('Target Attainment', options, [0.8 0.95]);
b = bar(categorical(targetTbl.metric), targetTbl.percent_attained, 0.6);
b.FaceColor = [0.30 0.60 0.30];
ylabel('Attainment (%)');
title('Probability of Target Attainment');
grid on;
ylim([0 100]);

for i = 1:height(targetTbl)
    text(i, targetTbl.percent_attained(i) + 1.5, sprintf('%.1f%%', targetTbl.percent_attained(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

save_fig(h, outDir, '05_target_attainment_bars', options);
end

function h = plot_tornado(sensAUC, sensCmax, outDir, options)
n = options.top_n_tornado;
sa = sensAUC(1:min(n,height(sensAUC)), :);
sc = sensCmax(1:min(n,height(sensCmax)), :);

h = new_fig('Tornado Sensitivity', options, [1.2 1.1]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
barh(categorical(sa.parameter), sa.abs_rho, 'FaceColor', [0.8 0.35 0.35]);
set(gca, 'YDir', 'reverse');
xlabel('|Spearman rho|'); ylabel('Parameter'); title('Top AUC Drivers'); grid on;

nexttile;
barh(categorical(sc.parameter), sc.abs_rho, 'FaceColor', [0.25 0.45 0.80]);
set(gca, 'YDir', 'reverse');
xlabel('|Spearman rho|'); ylabel('Parameter'); title('Top Cmax Drivers'); grid on;

save_fig(h, outDir, '06_tornado_sensitivity', options);
end

function h = plot_top_scatter_grid(T, sensAUC, sensCmax, outDir, options)
topN = options.top_n_scatter;
topSet = unique([string(sensAUC.parameter(1:min(topN,height(sensAUC)))); ...
                 string(sensCmax.parameter(1:min(topN,height(sensCmax))))]);

n = numel(topSet);
if n == 0
    h = [];
    return;
end

h = new_fig('Top Parameter Scatter Grid', options, [max(1.2, 0.5*n) 1.35]);
tiledlayout(2, n, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:n
    p = char(topSet(i));
    x = to_numeric(T.(p));

    nexttile;
    scatter(x, T.AUC_mg_h_L, 6, 'filled', 'MarkerFaceAlpha', 0.25, ...
        'MarkerFaceColor', [0.8 0.35 0.35], 'MarkerEdgeColor', 'none');
    lsline;
    xlabel(strrep(p, '_', '\_'), 'Interpreter', 'tex');
    ylabel('AUC');
    title(sprintf('%s vs AUC', p), 'Interpreter', 'none');
    grid on;

    nexttile;
    scatter(x, T.Cmax_uM, 6, 'filled', 'MarkerFaceAlpha', 0.25, ...
        'MarkerFaceColor', [0.25 0.45 0.80], 'MarkerEdgeColor', 'none');
    lsline;
    xlabel(strrep(p, '_', '\_'), 'Interpreter', 'tex');
    ylabel('Cmax');
    title(sprintf('%s vs Cmax', p), 'Interpreter', 'none');
    grid on;
    if options.auto_frame_cmax
        ylim(compute_auto_limits(T.Cmax_uM, options.axis_margin_fraction, true));
    end
end

save_fig(h, outDir, '07_top_scatter_grid', options);
end

function h = plot_correlation_heatmap(T, paramCols, outDir, options)
if numel(paramCols) > options.max_heatmap_vars
    paramCols = paramCols(1:options.max_heatmap_vars);
end
vars = [paramCols(:)', {'AUC_mg_h_L', 'Cmax_uM'}];
X = nan(height(T), numel(vars));
for i = 1:numel(vars)
    X(:,i) = to_numeric(T.(vars{i}));
end
R = corr(X, 'Type', 'Spearman', 'Rows', 'pairwise');

h = new_fig('Spearman Correlation Heatmap', options, [1.25 1.25]);
imagesc(R);
axis square;
colormap(parula);
colorbar;
caxis([-1 1]);

tickLbls = cellstr(strrep(shorten_labels(vars, options.tick_label_max_chars), '_', '\_'));
set(gca, 'XTick', 1:numel(vars), 'YTick', 1:numel(vars));
set(gca, 'XTickLabel', tickLbls, 'YTickLabel', tickLbls);
set(gca, 'TickLabelInterpreter', 'tex');
try
    xtickangle(45);
catch
    % Older MATLAB releases may not support xtickangle.
end
title('Spearman Correlation Matrix');

save_fig(h, outDir, '08_correlation_heatmap', options);
end

function h = plot_task_level_qc(T, outDir, options)
validTask = isfinite(T.chunk_task_id);
X = T(validTask,:);
if isempty(X)
    h = [];
    return;
end

[G, taskID] = findgroups(X.chunk_task_id);
nRows = splitapply(@numel, X.chunk_task_id, G);
aucMean = splitapply(@mean, X.AUC_mg_h_L, G);
cmaxMean = splitapply(@mean, X.Cmax_uM, G);

h = new_fig('Task-Level QC', options, [1.0 1.3]);
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile;
bar(taskID, nRows, 'FaceColor', [0.40 0.40 0.40]);
ylabel('Rows'); title('Rows per chunk_task_id'); grid on;

nexttile;
plot(taskID, aucMean, 'o-', 'Color', options.color_auc, 'LineWidth', 1.3);
ylabel('Mean AUC'); title('Task-level Mean AUC'); grid on;

nexttile;
plot(taskID, cmaxMean, 'o-', 'Color', options.color_cmax, 'LineWidth', 1.3);
ylabel('Mean Cmax'); xlabel('chunk_task_id'); title('Task-level Mean Cmax'); grid on;

save_fig(h, outDir, '09_task_level_qc', options);
end

function h = new_fig(name, options, sizeScale)
if nargin < 3 || isempty(sizeScale)
    sizeScale = [1 1];
end
figSize = options.fig_size_cm .* sizeScale;
h = figure('Name', name, 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 figSize(1) figSize(2)]);

% Keep on-screen and export canvas consistent.
set(h, 'PaperUnits', 'centimeters');
set(h, 'PaperSize', figSize);
set(h, 'PaperPosition', [0 0 figSize(1) figSize(2)]);
set(findall(h, '-property', 'FontSize'), 'FontSize', options.base_font_size);
set(findall(h, '-property', 'LineWidth'), 'LineWidth', options.axis_line_width);
end

function save_fig(h, outDir, baseName, options)
if isempty(h) || ~isgraphics(h)
    return;
end
pngPath = fullfile(outDir, [baseName '.png']);
figPath = fullfile(outDir, [baseName '.fig']);

% Re-apply publication size in case figure was resized interactively.
set(h, 'Units', 'centimeters');
pos = get(h, 'Position');
figSize = pos(3:4);
set(h, 'Position', pos);
set(h, 'PaperUnits', 'centimeters');
set(h, 'PaperSize', figSize);
set(h, 'PaperPosition', [0 0 figSize(1) figSize(2)]);

style_axes_publication(h, options);

% Disable axes toolbars/interactions before export to avoid toolbar artifacts.
disable_axes_toolbars(h);

drawnow;
if isfield(options, 'use_exportgraphics') && options.use_exportgraphics
    try
        exportgraphics(h, pngPath, 'Resolution', options.export_png_dpi);
    catch ME
        warning('analyze_MC10k_comprehensive:ExportGraphicsFailed', ...
            'exportgraphics failed for %s (%s). Falling back to print().', baseName, ME.message);
        try
            print(h, pngPath, '-dpng', sprintf('-r%d', options.export_png_dpi));
        catch ME2
            warning('analyze_MC10k_comprehensive:PrintFallbackFailed', ...
                'PNG export fallback also failed for %s: %s', baseName, ME2.message);
        end
    end
else
    try
        print(h, pngPath, '-dpng', sprintf('-r%d', options.export_png_dpi));
    catch ME
        warning('analyze_MC10k_comprehensive:PrintPrimaryFailed', ...
            'print() failed for %s (%s). Falling back to exportgraphics.', baseName, ME.message);
        try
            exportgraphics(h, pngPath, 'Resolution', options.export_png_dpi);
        catch ME2
            warning('analyze_MC10k_comprehensive:ExportGraphicsFallbackFailed', ...
                'PNG export fallback also failed for %s: %s', baseName, ME2.message);
        end
    end
end

try
    savefig(h, figPath);
catch ME
    warning('analyze_MC10k_comprehensive:SaveFigFailed', ...
        'savefig failed for %s: %s', baseName, ME.message);
end

if isfield(options, 'export_pdf') && options.export_pdf
    pdfPath = fullfile(outDir, [baseName '.pdf']);
    if isfield(options, 'use_exportgraphics') && options.use_exportgraphics
        ct = 'image';
        if isfield(options, 'pdf_content_type') && ~isempty(options.pdf_content_type)
            ct = char(lower(string(options.pdf_content_type)));
            if ~strcmp(ct, 'image') && ~strcmp(ct, 'vector')
                ct = 'image';
            end
        end
        try
            exportgraphics(h, pdfPath, 'ContentType', ct);
        catch
            try
                print(h, pdfPath, '-dpdf', '-painters');
            catch
            end
        end
    else
        try
            print(h, pdfPath, '-dpdf', '-painters');
        catch
            try
                exportgraphics(h, pdfPath, 'ContentType', 'image');
            catch
            end
        end
    end
end
end

function ok = try_violinplot(y, colorVal)
ok = false;
y = y(isfinite(y));
if isempty(y)
    return;
end

% Different violinplot implementations accept different signatures.
try
    violinplot(y);
    style_last_violin(colorVal);
    ok = true;
    return;
catch
end

try
    violinplot(y, ones(size(y)));
    style_last_violin(colorVal);
    ok = true;
    return;
catch
end

try
    violinplot(y, ones(size(y)), 'ViolinColor', colorVal);
    ok = true;
catch
    ok = false;
end
end

function style_last_violin(colorVal)
try
    ax = gca;
    h = findobj(ax, 'Type', 'Patch');
    if ~isempty(h)
        set(h(1), 'FaceColor', colorVal, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
catch
    % Styling is best-effort only.
end
end

function disable_axes_toolbars(figHandle)
axs = findall(figHandle, 'Type', 'axes');
for i = 1:numel(axs)
    try
        if isprop(axs(i), 'Toolbar') && ~isempty(axs(i).Toolbar)
            axs(i).Toolbar.Visible = 'off';
        end
    catch
        % Keep export resilient across MATLAB versions.
    end
end
end

function lim = compute_auto_limits(x, marginFrac, clampToZero) %#ok<DEFNU>
if nargin < 2 || isempty(marginFrac)
    marginFrac = 0.08;
end
if nargin < 3
    clampToZero = false;
end

x = to_numeric(x);
x = x(isfinite(x));
if isempty(x)
    lim = [0 1];
    return;
end

xMin = min(x);
xMax = max(x);
if xMax <= xMin
    span = max(abs(xMax), 1);
    pad = 0.1 * span;
else
    span = xMax - xMin;
    pad = max(span * marginFrac, eps(max(abs([xMin, xMax]))));
end

lo = xMin - pad;
hi = xMax + pad;
if clampToZero
    lo = max(0, lo);
end
if hi <= lo
    hi = lo + 1;
end

lim = [lo hi];
end

function out = shorten_labels(in, maxChars)
out = string(in);
for i = 1:numel(out)
    if strlength(out(i)) > maxChars
        out(i) = extractBefore(out(i), maxChars - 2) + "..";
    end
end
end

function style_axes_publication(figHandle, options)
axs = findall(figHandle, 'Type', 'axes');
for i = 1:numel(axs)
    try
        set(axs(i), 'FontSize', options.base_font_size, 'LineWidth', options.axis_line_width, 'Box', 'off');
        if ~isempty(axs(i).XLabel)
            axs(i).XLabel.FontSize = options.label_font_size;
        end
        if ~isempty(axs(i).YLabel)
            axs(i).YLabel.FontSize = options.label_font_size;
        end
        if ~isempty(axs(i).Title)
            axs(i).Title.FontSize = options.title_font_size;
            axs(i).Title.FontWeight = 'bold';
        end
    catch
    end
end
lg = findall(figHandle, 'Type', 'Legend');
for i = 1:numel(lg)
    try
        lg(i).FontSize = options.legend_font_size;
        lg(i).Box = 'off';
    catch
    end
end
end

function save_qc_report_txt(p, T, summaryTbl, targetTbl, results, options)
fid = fopen(p, 'w');
if fid < 0
    warning('Could not write report file: %s', p);
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'MC10k Comprehensive Analysis Report\n');
fprintf(fid, 'Generated: %s\n\n', char(datetime('now')));

fprintf(fid, 'Input file: %s\n', results.input_file);
fprintf(fid, 'Output dir: %s\n\n', results.output_dir);

fprintf(fid, 'Data completeness\n');
fprintf(fid, '  Rows raw:   %d\n', results.n_rows_raw);
fprintf(fid, '  Rows clean: %d\n', results.n_rows_clean);
fprintf(fid, '  Dropped:    %d\n\n', results.n_rows_dropped);

fprintf(fid, 'Targets\n');
fprintf(fid, '  AUC target range:  [%.3g, %.3g]\n', options.auc_target(1), options.auc_target(2));
fprintf(fid, '  Cmax target range: [%.3g, %.3g]\n\n', options.cmax_target(1), options.cmax_target(2));

fprintf(fid, 'Summary statistics\n');
for i = 1:height(summaryTbl)
    fprintf(fid, '  %s: mean=%.4g, median=%.4g, IQR=[%.4g, %.4g], p05=%.4g, p95=%.4g\n', ...
        char(summaryTbl.metric(i)), summaryTbl.mean(i), summaryTbl.median(i), summaryTbl.q25(i), summaryTbl.q75(i), ...
        summaryTbl.p05(i), summaryTbl.p95(i));
end
fprintf(fid, '\n');

fprintf(fid, 'Target attainment\n');
for i = 1:height(targetTbl)
    fprintf(fid, '  %s: %.2f%%\n', char(targetTbl.metric(i)), targetTbl.percent_attained(i));
end
fprintf(fid, '\n');

if ismember('chunk_task_id', T.Properties.VariableNames)
    x = to_numeric(T.chunk_task_id);
    x = x(isfinite(x));
    fprintf(fid, 'Task metadata\n');
    fprintf(fid, '  Distinct chunk_task_id: %d\n', numel(unique(x)));
end

end
