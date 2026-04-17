function results = compare_MC10k_regimens(referenceCsv, infusionCsv, outDir, options)
% compare_MC10k_regimens
% Compare two merged MC regimen outputs (e.g., baseline vs infusion-dominant)
% and generate publication-ready comparative figures and statistics.
%
% Outputs:
%   - Comparative endpoint distribution figures
%   - Comparative violin/box figures
%   - Target attainment comparison bars
%   - AUC/Cmax scatter comparison
%   - Sensitivity delta tornado (Spearman rho differences)
%   - Percentile shift plots (emergent trend view)
%   - Correlation difference heatmap
%   - CSV tables and text report
%
% Usage:
%   compare_MC10k_regimens
%   compare_MC10k_regimens('MC_results/MC10k_chunks/MC_10k_merged_outputs.csv', ...
%                          'MC_results/MC10k_infusion/MC_10k_merged_outputs.csv')
%
%   options = struct();
%   options.reference_name = 'Bolus-dominant';
%   options.infusion_name = 'Infusion-dominant';
%   compare_MC10k_regimens(refCsv, infCsv, 'comparison_pack', options)

if nargin < 1 || isempty(referenceCsv)
    referenceCsv = fullfile('MC_results', 'MC10k_chunks', 'MC_10k_merged_outputs.csv');
end
if nargin < 2 || isempty(infusionCsv)
    infusionCsv = fullfile('MC_results', 'MC10k_infusion_chunks', 'MC_10k_merged_outputs.csv');
end
if nargin < 3 || isempty(outDir)
    stamp = char(string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    outDir = fullfile(pwd, ['comparison_pack_' stamp]);
end
if nargin < 4 || isempty(options)
    options = struct();
end

options = apply_defaults(options);
if ~isfile(referenceCsv)
    error('Reference CSV not found: %s', referenceCsv);
end
if ~isfile(infusionCsv)
    error('Infusion CSV not found: %s', infusionCsv);
end
if ~isfolder(outDir)
    mkdir(outDir);
end

fprintf('\n=== Regimen Comparison Analysis ===\n');
fprintf('Reference: %s\n', referenceCsv);
fprintf('Infusion : %s\n', infusionCsv);
fprintf('Output   : %s\n\n', outDir);

R = readtable(referenceCsv, 'VariableNamingRule', 'preserve');
I = readtable(infusionCsv, 'VariableNamingRule', 'preserve');

required = {'AUC_mg_h_L','Cmax_uM'};
for k = 1:numel(required)
    if ~ismember(required{k}, R.Properties.VariableNames)
        error('Reference file missing required column: %s', required{k});
    end
    if ~ismember(required{k}, I.Properties.VariableNames)
        error('Infusion file missing required column: %s', required{k});
    end
end

R = clean_table(R);
I = clean_table(I);

results = struct();
results.reference_file = referenceCsv;
results.infusion_file = infusionCsv;
results.output_dir = outDir;
results.reference_n = height(R);
results.infusion_n = height(I);

fprintf('Reference rows (clean): %d\n', results.reference_n);
fprintf('Infusion rows  (clean): %d\n\n', results.infusion_n);

% Common numeric parameter columns
rMask = varfun(@isnumeric, R, 'OutputFormat', 'uniform');
iMask = varfun(@isnumeric, I, 'OutputFormat', 'uniform');
rNum = R.Properties.VariableNames(rMask);
iNum = I.Properties.VariableNames(iMask);

excludeCols = {'AUC_mg_h_L','Cmax_uM','sample_local_index','chunk_task_id'};
commonParams = intersect(setdiff(rNum, excludeCols, 'stable'), setdiff(iNum, excludeCols, 'stable'), 'stable');
results.common_parameters = commonParams;

% Summary and target attainment tables
summaryTbl = build_regimen_summary(R, I, options);
writetable(summaryTbl, fullfile(outDir, 'regimen_summary_statistics.csv'));

targetTbl = build_target_attainment_comparison(R, I, options);
writetable(targetTbl, fullfile(outDir, 'target_attainment_comparison.csv'));

% Sensitivity comparisons
if ~isempty(commonParams)
    [sensAucRef, sensAucInf] = spearman_sensitivity_pair(R, I, commonParams, 'AUC_mg_h_L');
    [sensCmaxRef, sensCmaxInf] = spearman_sensitivity_pair(R, I, commonParams, 'Cmax_uM');

    sensAucDelta = join_sensitivity(sensAucRef, sensAucInf, options.reference_name, options.infusion_name);
    sensCmaxDelta = join_sensitivity(sensCmaxRef, sensCmaxInf, options.reference_name, options.infusion_name);

    writetable(sensAucDelta, fullfile(outDir, 'sensitivity_delta_auc.csv'));
    writetable(sensCmaxDelta, fullfile(outDir, 'sensitivity_delta_cmax.csv'));
else
    sensAucDelta = table();
    sensCmaxDelta = table();
end

% Emergent percentile shift tables
pctTblA = build_percentile_shift_table(R.AUC_mg_h_L, I.AUC_mg_h_L, 'AUC_mg_h_L');
pctTblC = build_percentile_shift_table(R.Cmax_uM, I.Cmax_uM, 'Cmax_uM');
writetable(pctTblA, fullfile(outDir, 'percentile_shift_auc.csv'));
writetable(pctTblC, fullfile(outDir, 'percentile_shift_cmax.csv'));

% Figure pack
fig_endpoint_distributions(R, I, outDir, options);
fig_violin_or_box(R, I, outDir, options);
fig_auc_cmax_scatter(R, I, outDir, options);
fig_target_attainment(targetTbl, outDir, options);
fig_percentile_shift(pctTblA, pctTblC, outDir, options);

if ~isempty(commonParams)
    fig_tornado_delta(sensAucDelta, sensCmaxDelta, outDir, options);
    fig_correlation_delta_heatmap(R, I, commonParams, outDir, options);
end

write_report(fullfile(outDir, 'comparison_report.txt'), targetTbl, pctTblA, pctTblC, results, options);
save(fullfile(outDir, 'comparison_workspace.mat'), 'results', 'options', 'summaryTbl', 'targetTbl', ...
    'sensAucDelta', 'sensCmaxDelta', 'pctTblA', 'pctTblC');

fprintf('Comparison analysis complete.\nOutput folder: %s\n', outDir);

end

function options = apply_defaults(options)
def = struct();
def.reference_name = 'Reference regimen';
def.infusion_name = 'Infusion-dominant regimen';
def.auc_target = [20 30];
def.cmax_target = [300 500];
def.alpha = 0.05;
def.top_n_tornado = 12;
def.percentiles = (1:1:99)';
def.export_png_dpi = 450;

% Publication sizing, A4-aware.
def.a4_width_cm = 21.0;
def.figure_width_fraction_a4 = 0.8;
def.figure_height_to_width_ratio = 0.62;
def.figure_min_height_cm = 7.4;
def.fig_size_cm = [];

% Style
ndef = fieldnames(def);
for i = 1:numel(ndef)
    if ~isfield(options, ndef{i}) || isempty(options.(ndef{i}))
        options.(ndef{i}) = def.(ndef{i});
    end
end

if ~isfield(options, 'fig_size_cm') || isempty(options.fig_size_cm)
    w = options.a4_width_cm * options.figure_width_fraction_a4;
    h = max(options.figure_min_height_cm, w * options.figure_height_to_width_ratio);
    options.fig_size_cm = [w h];
end

options.color_ref = [0.2 0.45 0.8];
options.color_inf = [0.82 0.30 0.30];
options.color_delta = [0.1 0.1 0.1];
end

function T = clean_table(T)
T.AUC_mg_h_L = to_num(T.AUC_mg_h_L);
T.Cmax_uM = to_num(T.Cmax_uM);
keep = isfinite(T.AUC_mg_h_L) & isfinite(T.Cmax_uM);
T = T(keep, :);

vars = T.Properties.VariableNames;
for i = 1:numel(vars)
    if ~isnumeric(T.(vars{i}))
        x = to_num(T.(vars{i}));
        if any(isfinite(x))
            T.(vars{i}) = x;
        end
    end
end
end

function x = to_num(v)
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

function summaryTbl = build_regimen_summary(R, I, options)
metrics = {'AUC_mg_h_L','Cmax_uM'};
name = {};
regimen = {};
n = [];
meanv = [];
sd = [];
medianv = [];
q25 = [];
q75 = [];
p05 = [];
p95 = [];

for m = 1:numel(metrics)
    for g = 1:2
        if g == 1
            x = R.(metrics{m});
            lab = options.reference_name;
        else
            x = I.(metrics{m});
            lab = options.infusion_name;
        end
        x = x(isfinite(x));
        name{end+1,1} = metrics{m}; %#ok<AGROW>
        regimen{end+1,1} = lab; %#ok<AGROW>
        n(end+1,1) = numel(x); %#ok<AGROW>
        meanv(end+1,1) = mean(x); %#ok<AGROW>
        sd(end+1,1) = std(x); %#ok<AGROW>
        medianv(end+1,1) = median(x); %#ok<AGROW>
        q25(end+1,1) = prctile(x,25); %#ok<AGROW>
        q75(end+1,1) = prctile(x,75); %#ok<AGROW>
        p05(end+1,1) = prctile(x,5); %#ok<AGROW>
        p95(end+1,1) = prctile(x,95); %#ok<AGROW>
    end
end

summaryTbl = table(string(name), string(regimen), n, meanv, sd, medianv, q25, q75, p05, p95, ...
    'VariableNames', {'metric','regimen','n','mean','sd','median','q25','q75','p05','p95'});
end

function targetTbl = build_target_attainment_comparison(R, I, options)
aucR = R.AUC_mg_h_L; cmaxR = R.Cmax_uM;
aucI = I.AUC_mg_h_L; cmaxI = I.Cmax_uM;

AUC_R = mean(aucR >= options.auc_target(1) & aucR <= options.auc_target(2)) * 100;
AUC_I = mean(aucI >= options.auc_target(1) & aucI <= options.auc_target(2)) * 100;
CMX_R = mean(cmaxR >= options.cmax_target(1) & cmaxR <= options.cmax_target(2)) * 100;
CMX_I = mean(cmaxI >= options.cmax_target(1) & cmaxI <= options.cmax_target(2)) * 100;
BTH_R = mean((aucR >= options.auc_target(1) & aucR <= options.auc_target(2)) & ...
             (cmaxR >= options.cmax_target(1) & cmaxR <= options.cmax_target(2))) * 100;
BTH_I = mean((aucI >= options.auc_target(1) & aucI <= options.auc_target(2)) & ...
             (cmaxI >= options.cmax_target(1) & cmaxI <= options.cmax_target(2))) * 100;

targetTbl = table();
targetTbl.metric = ["AUC"; "Cmax"; "Both"];
targetTbl.reference_percent = [AUC_R; CMX_R; BTH_R];
targetTbl.infusion_percent = [AUC_I; CMX_I; BTH_I];
targetTbl.delta_percent = targetTbl.infusion_percent - targetTbl.reference_percent;
end

function [Sref, Sinf] = spearman_sensitivity_pair(R, I, params, outcome)
Sref = table();
Sinf = table();
Sref.parameter = string(params(:));
Sinf.parameter = string(params(:));
Sref.rho = nan(numel(params),1);
Sinf.rho = nan(numel(params),1);
Sref.p_value = nan(numel(params),1);
Sinf.p_value = nan(numel(params),1);

yR = R.(outcome);
yI = I.(outcome);

for i = 1:numel(params)
    xR = to_num(R.(params{i}));
    xI = to_num(I.(params{i}));
    [Sref.rho(i), Sref.p_value(i)] = corr(xR, yR, 'Type', 'Spearman', 'Rows', 'complete');
    [Sinf.rho(i), Sinf.p_value(i)] = corr(xI, yI, 'Type', 'Spearman', 'Rows', 'complete');
end
end

function Tout = join_sensitivity(Sref, Sinf, refName, infName)
Tout = table();
Tout.parameter = Sref.parameter;
Tout.([sanitize_name(refName) '_rho']) = Sref.rho;
Tout.([sanitize_name(infName) '_rho']) = Sinf.rho;
Tout.delta_rho = Sinf.rho - Sref.rho;
Tout.delta_abs_rho = abs(Sinf.rho) - abs(Sref.rho);
Tout = sortrows(Tout, 'delta_abs_rho', 'descend');
end

function s = sanitize_name(x)
s = lower(char(string(x)));
s = regexprep(s, '[^a-z0-9]+', '_');
s = regexprep(s, '^_+|_+$', '');
if isempty(s)
    s = 'regimen';
end
end

function T = build_percentile_shift_table(xRef, xInf, metricName)
p = (1:1:99)';
qRef = prctile(xRef, p);
qInf = prctile(xInf, p);
T = table();
T.metric = repmat(string(metricName), numel(p), 1);
T.percentile = p;
T.reference_value = qRef(:);
T.infusion_value = qInf(:);
T.absolute_shift = T.infusion_value - T.reference_value;
T.relative_shift_percent = 100 * T.absolute_shift ./ max(abs(T.reference_value), eps);
end

function h = fig_endpoint_distributions(R, I, outDir, options)
h = new_pub_fig('Endpoint Comparison Distributions', options);
tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
histogram(R.AUC_mg_h_L, 40, 'FaceColor', options.color_ref, 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
histogram(I.AUC_mg_h_L, 40, 'FaceColor', options.color_inf, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
xline(options.auc_target(1), '--k'); xline(options.auc_target(2), '--k');
legend(options.reference_name, options.infusion_name, 'Location', 'best');
title('AUC Histogram Overlay'); xlabel('AUC (mg h/L)'); ylabel('Count'); grid on;

nexttile;
ecdf(R.AUC_mg_h_L); hold on; ecdf(I.AUC_mg_h_L);
legend(options.reference_name, options.infusion_name, 'Location', 'best');
title('AUC CDF Comparison'); xlabel('AUC (mg h/L)'); ylabel('F(x)'); grid on;

nexttile;
histogram(R.Cmax_uM, 40, 'FaceColor', options.color_ref, 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
histogram(I.Cmax_uM, 40, 'FaceColor', options.color_inf, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
xline(options.cmax_target(1), '--k'); xline(options.cmax_target(2), '--k');
legend(options.reference_name, options.infusion_name, 'Location', 'best');
title('Cmax Histogram Overlay'); xlabel('Cmax (uM)'); ylabel('Count'); grid on;

nexttile;
ecdf(R.Cmax_uM); hold on; ecdf(I.Cmax_uM);
legend(options.reference_name, options.infusion_name, 'Location', 'best');
title('Cmax CDF Comparison'); xlabel('Cmax (uM)'); ylabel('F(x)'); grid on;

save_pub_fig(h, outDir, '01_endpoint_distribution_comparison', options);
end

function h = fig_violin_or_box(R, I, outDir, options)
h = new_pub_fig('Violin Comparison', options);
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot_violin_or_box_pair(R.AUC_mg_h_L, I.AUC_mg_h_L, options);
title('AUC'); ylabel('AUC (mg h/L)'); grid on;

nexttile;
plot_violin_or_box_pair(R.Cmax_uM, I.Cmax_uM, options);
title('Cmax'); ylabel('Cmax (uM)'); grid on;

save_pub_fig(h, outDir, '02_violin_box_comparison', options);
end

function plot_violin_or_box_pair(xRef, xInf, options)
hasViolin = exist('violinplot', 'file') == 2;
if hasViolin
    ok1 = try_violin_single(xRef, 1, options.color_ref);
    ok2 = try_violin_single(xInf, 2, options.color_inf);
    if ~(ok1 && ok2)
        fallback_boxplot(xRef, xInf, options);
        return;
    end
    xlim([0.5 2.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {options.reference_name, options.infusion_name});
else
    fallback_boxplot(xRef, xInf, options);
end
end

function ok = try_violin_single(y, groupPos, colorVal)
ok = false;
y = y(isfinite(y));
if isempty(y)
    return;
end

try
    violinplot(y, groupPos);
    style_last_violin(colorVal);
    ok = true;
catch
    ok = false;
end
end

function style_last_violin(colorVal)
try
    ax = gca;
    hPatch = findobj(ax, 'Type', 'Patch');
    if ~isempty(hPatch)
        set(hPatch(1), 'FaceColor', colorVal, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
catch
end
end

function fallback_boxplot(xRef, xInf, options)
xRef = xRef(isfinite(xRef));
xInf = xInf(isfinite(xInf));
try
    g = [ones(numel(xRef),1); 2*ones(numel(xInf),1)];
    y = [xRef; xInf];
    boxchart(g, y, 'BoxFaceColor', options.color_ref);
catch
    boxplot([xRef(:), xInf(:)], 'Labels', {options.reference_name, options.infusion_name});
end
set(gca, 'XTick', [1 2], 'XTickLabel', {options.reference_name, options.infusion_name});
end

function h = fig_auc_cmax_scatter(R, I, outDir, options)
h = new_pub_fig('AUC Cmax Scatter Comparison', options);
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
scatter(R.AUC_mg_h_L, R.Cmax_uM, 7, 'MarkerFaceColor', options.color_ref, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.18); hold on;
xline(options.auc_target(1), '--k'); xline(options.auc_target(2), '--k');
yline(options.cmax_target(1), '--k'); yline(options.cmax_target(2), '--k');
title(options.reference_name, 'Interpreter', 'none');
xlabel('AUC (mg h/L)'); ylabel('Cmax (uM)'); grid on;

nexttile;
scatter(I.AUC_mg_h_L, I.Cmax_uM, 7, 'MarkerFaceColor', options.color_inf, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.18); hold on;
xline(options.auc_target(1), '--k'); xline(options.auc_target(2), '--k');
yline(options.cmax_target(1), '--k'); yline(options.cmax_target(2), '--k');
title(options.infusion_name, 'Interpreter', 'none');
xlabel('AUC (mg h/L)'); ylabel('Cmax (uM)'); grid on;

save_pub_fig(h, outDir, '03_auc_cmax_scatter_comparison', options);
end

function h = fig_target_attainment(T, outDir, options)
h = new_pub_fig('Target Attainment Comparison', options);
X = [T.reference_percent, T.infusion_percent];
bar(categorical(T.metric), X, 0.72);
colormap([options.color_ref; options.color_inf]);
legend(options.reference_name, options.infusion_name, 'Location', 'best', 'Interpreter', 'none');
ylabel('Attainment (%)'); ylim([0 100]); grid on;
title('Target Attainment Comparison');

for i = 1:height(T)
    text(i, max(X(i,:)) + 1.2, sprintf('Delta %+0.1f%%', T.delta_percent(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

save_pub_fig(h, outDir, '04_target_attainment_comparison', options);
end

function h = fig_percentile_shift(pA, pC, outDir, options)
h = new_pub_fig('Percentile Shift Emergent Trends', options);
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(pA.percentile, pA.absolute_shift, '-', 'Color', options.color_delta, 'LineWidth', 1.8); hold on;
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('Percentile'); ylabel('Infusion - Reference');
title('AUC Percentile Shift'); grid on;

nexttile;
plot(pC.percentile, pC.absolute_shift, '-', 'Color', options.color_delta, 'LineWidth', 1.8); hold on;
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('Percentile'); ylabel('Infusion - Reference');
title('Cmax Percentile Shift'); grid on;

save_pub_fig(h, outDir, '05_percentile_shift_emergent_trends', options);
end

function h = fig_tornado_delta(sA, sC, outDir, options)
n = min(options.top_n_tornado, height(sA));
a = sA(1:n,:);
c = sC(1:min(options.top_n_tornado, height(sC)), :);

h = new_pub_fig('Sensitivity Delta Tornado', options);
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
barh(categorical(a.parameter), a.delta_abs_rho, 'FaceColor', [0.45 0.45 0.45]);
set(gca, 'YDir', 'reverse');
xlabel('Delta |rho| (Infusion - Reference)');
title('AUC Driver Shift'); grid on;

nexttile;
barh(categorical(c.parameter), c.delta_abs_rho, 'FaceColor', [0.25 0.25 0.25]);
set(gca, 'YDir', 'reverse');
xlabel('Delta |rho| (Infusion - Reference)');
title('Cmax Driver Shift'); grid on;

save_pub_fig(h, outDir, '06_sensitivity_delta_tornado', options);
end

function h = fig_correlation_delta_heatmap(R, I, params, outDir, options)
vars = [params(:)', {'AUC_mg_h_L','Cmax_uM'}];
XR = nan(height(R), numel(vars));
XI = nan(height(I), numel(vars));
for i = 1:numel(vars)
    XR(:,i) = to_num(R.(vars{i}));
    XI(:,i) = to_num(I.(vars{i}));
end
RR = corr(XR, 'Type', 'Spearman', 'Rows', 'pairwise');
RI = corr(XI, 'Type', 'Spearman', 'Rows', 'pairwise');
D = RI - RR;

h = new_pub_fig('Correlation Delta Heatmap', options);
imagesc(D);
axis square;
colormap(parula);
colorbar;
caxis([-0.5 0.5]);

lbl = cellstr(strrep(vars, '_', '\_'));
set(gca, 'XTick', 1:numel(vars), 'YTick', 1:numel(vars), 'XTickLabel', lbl, 'YTickLabel', lbl, ...
    'TickLabelInterpreter', 'tex');
try
    xtickangle(45);
catch
end
title('Spearman Correlation Shift (Infusion - Reference)');

save_pub_fig(h, outDir, '07_correlation_delta_heatmap', options);
end

function h = new_pub_fig(figName, options)
h = figure('Name', figName, 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 options.fig_size_cm(1) options.fig_size_cm(2)]);
set(h, 'PaperUnits', 'centimeters');
set(h, 'PaperSize', options.fig_size_cm);
set(h, 'PaperPosition', [0 0 options.fig_size_cm(1) options.fig_size_cm(2)]);
end

function save_pub_fig(h, outDir, baseName, options)
if isempty(h) || ~isgraphics(h)
    return;
end
pngPath = fullfile(outDir, [baseName '.png']);
figPath = fullfile(outDir, [baseName '.fig']);

set(h, 'Units', 'centimeters');
pos = get(h, 'Position');
pos(3:4) = options.fig_size_cm;
set(h, 'Position', pos);
set(h, 'PaperUnits', 'centimeters');
set(h, 'PaperSize', options.fig_size_cm);
set(h, 'PaperPosition', [0 0 options.fig_size_cm(1) options.fig_size_cm(2)]);

disable_axes_toolbars(h);
exportgraphics(h, pngPath, 'Resolution', options.export_png_dpi);
savefig(h, figPath);
end

function disable_axes_toolbars(figHandle)
axs = findall(figHandle, 'Type', 'axes');
for i = 1:numel(axs)
    try
        if isprop(axs(i), 'Toolbar') && ~isempty(axs(i).Toolbar)
            axs(i).Toolbar.Visible = 'off';
        end
    catch
    end
end
end

function write_report(outPath, targetTbl, pctA, pctC, results, options)
fid = fopen(outPath, 'w');
if fid < 0
    warning('Could not write report: %s', outPath);
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Regimen Comparison Report\n');
fprintf(fid, 'Generated: %s\n\n', char(datetime('now')));

fprintf(fid, 'Reference file: %s\n', results.reference_file);
fprintf(fid, 'Infusion file : %s\n', results.infusion_file);
fprintf(fid, 'Output dir    : %s\n\n', results.output_dir);

fprintf(fid, 'Sample counts\n');
fprintf(fid, '  %s: %d\n', options.reference_name, results.reference_n);
fprintf(fid, '  %s: %d\n\n', options.infusion_name, results.infusion_n);

fprintf(fid, 'Target attainment deltas\n');
for i = 1:height(targetTbl)
    fprintf(fid, '  %s delta: %+0.2f%%\n', char(targetTbl.metric(i)), targetTbl.delta_percent(i));
end
fprintf(fid, '\n');

fprintf(fid, 'Most shifted AUC percentiles (absolute)\n');
[~, ia] = maxk(abs(pctA.absolute_shift), 5);
for k = 1:numel(ia)
    i = ia(k);
    fprintf(fid, '  P%02d: shift %+0.4g\n', pctA.percentile(i), pctA.absolute_shift(i));
end
fprintf(fid, '\n');

fprintf(fid, 'Most shifted Cmax percentiles (absolute)\n');
[~, ic] = maxk(abs(pctC.absolute_shift), 5);
for k = 1:numel(ic)
    i = ic(k);
    fprintf(fid, '  P%02d: shift %+0.4g\n', pctC.percentile(i), pctC.absolute_shift(i));
end
fprintf(fid, '\n');

fprintf(fid, 'Summary statistics table exported to regimen_summary_statistics.csv\n');
fprintf(fid, 'Target table exported to target_attainment_comparison.csv\n');
end
