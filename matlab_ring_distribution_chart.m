% Stacked bar chart of aggregated radial rings using the tuple from
% crc_tumor_caf_model.py (RING_DATA), without loading cell location CSV.

outputDir = "d:\PhysiCell";
analysisDate = datestr(now, 'yyyy-mm-dd');

% Group adjacent 50 um rings into larger bins for clearer interpretation.
ringsPerGroup = 5;  % 5 x 50um = 250um wide bins

if ~isfolder(outputDir)
    mkdir(outputDir);
end

% RING_DATA columns from crc_tumor_caf_model.py:
% [ring_num, CSC, CP, CD, CD8, M2, CAF]
ringNum = (1:30)';
countsByRing = [ ...
      0,    0,    2,   0,   0,   0;
      4,    0,   17,   0,   1,   0;
     21,    0,   56,   0,   4,   0;
     53,    0,  125,   0,   9,   0;
     88,    0,  231,   0,  16,   0;
    103,    0,  374,   0,  28,   0;
     87,    2,  552,   0,  43,   0;
     54,    6,  756,   0,  61,   0;
     25,   15,  971,   0,  84,   0;
      9,   34, 1179,   0, 110,   0;
      2,   70, 1362,   0, 139,   0;
      0,  130, 1503,   0, 170,   0;
      0,  223, 1589,   1, 201,   0;
      0,  356, 1613,   4, 230,   1;
      0,  527, 1575,  12, 258,   3;
      1,  727, 1482,  31, 281,   7;
     26,  934, 1345,  68, 299,  15;
    224, 1118, 1178, 135, 311,  30;
    681, 1249,  997, 235, 317,  55;
    756, 1304,  816, 365, 315,  95;
    307, 1272,  646, 503, 308, 154;
     45, 1160,  495, 619, 294, 231;
      2,  989,  367, 678, 275, 325;
      0,  789,  264, 662, 253, 428;
      0,  589,  184, 576, 228, 527;
      0,  412,  124, 447, 201, 608;
      0,  270,   81, 309, 175, 657;
      0,  165,   51, 191, 149, 665;
      0,   95,   31, 105, 124, 630;
      0,   51,   19,  51, 102, 560];

typeCats = ["Cancer_stem_cell", "cancer_proliferating", "cancer_differentiated", ...
            "CD8_T_cell", "M2_macrophage", "CAF"];

% Professional, colorblind-friendly palette (one fixed color per cell type)
typeColors = [
    0.16 0.37 0.59;  % Cancer_stem_cell (deep blue)
    0.85 0.37 0.01;  % cancer_proliferating (orange)
    0.00 0.62 0.45;  % cancer_differentiated (teal-green)
    0.57 0.44 0.86;  % CD8_T_cell (purple)
    0.90 0.60 0.00;  % M2_macrophage (amber)
    0.20 0.20 0.20;  % CAF (charcoal)
];

nRings = numel(ringNum);
if mod(nRings, ringsPerGroup) ~= 0
    error('ringsPerGroup (%d) must divide total rings (%d).', ringsPerGroup, nRings);
end

nGroups = nRings / ringsPerGroup;
countMat = zeros(nGroups, numel(typeCats));
ringLabels = strings(nGroups, 1);

for g = 1:nGroups
    idx0 = (g-1) * ringsPerGroup + 1;
    idx1 = g * ringsPerGroup;

    countMat(g, :) = sum(countsByRing(idx0:idx1, :), 1);

    innerUm = (idx0 - 1) * 50;
    outerUm = idx1 * 50;
    ringLabels(g) = sprintf('%d-%d µm', innerUm, outerUm);
end

% Plot stacked bar chart
fig = figure('Color', 'w', 'Position', [120 120 1150 650]);
ringCat = categorical(ringLabels, ringLabels, 'Ordinal', true);
b = bar(ringCat, countMat, 'stacked', 'EdgeColor', 'none');
for c = 1:numel(typeCats)
    b(c).FaceColor = typeColors(c, :);
end

xlabel('Aggregated radial ring from tumor center');
ylabel('Cell count');
title(sprintf('Cell Counts by Larger Rings (%d rings per bin)', ringsPerGroup));
grid on;
legend(typeCats, 'Interpreter', 'none', 'Location', 'eastoutside');

% Add total count annotations above each bar
totals = sum(countMat, 2);
yTop = max(totals);
for g = 1:nGroups
    text(g, totals(g) + max(1, 0.015 * yTop), sprintf('%d', round(totals(g))), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end

set_publication_figure_style(fig, 'font_size', 11, 'font_name', 'Arial');

totalCells = sum(countsByRing, 'all');
add_figure_metadata('crc_tumor_caf_model.py:RING_DATA', totalCells, analysisDate);

quickPng = fullfile(outputDir, 'ring_counts_stacked.png');
exportgraphics(fig, quickPng, 'Resolution', 300);
fprintf('Saved PNG only: %s\n', quickPng);

close(fig);

% Line + scatter plot: per-ring counts (50 um rings) for each cell type
fig2 = figure('Color', 'w', 'Position', [140 140 1150 650]);
hold on;

ringCenterUm = (ringNum - 0.5) * 50;
lineHandles = gobjects(numel(typeCats), 1);
for c = 1:numel(typeCats)
    y = countsByRing(:, c);
    lineHandles(c) = plot(ringCenterUm, y, '-', 'LineWidth', 2.0, ...
        'Color', typeColors(c, :), 'DisplayName', typeCats(c));
    scatter(ringCenterUm, y, 30, ...
        'MarkerFaceColor', typeColors(c, :), ...
        'MarkerEdgeColor', 'white', ...
        'LineWidth', 0.5, ...
        'MarkerFaceAlpha', 0.9, ...
        'HandleVisibility', 'off');
end

hold off;
xlabel('Ring center radius (\mum)');
ylabel('Cell count per 50 \mum ring');
title('Per-Ring Cell Counts by Cell Type (Line + Scatter)');
grid on;
legend(lineHandles, typeCats, 'Interpreter', 'none', 'Location', 'eastoutside');
set_publication_figure_style(fig2, 'font_size', 11, 'font_name', 'Arial');

% Apply metadata and save PNG only (per request)
add_figure_metadata('crc_tumor_caf_model.py:RING_DATA', totalCells, analysisDate);

lineScatterPng = fullfile(outputDir, 'ring_counts_line_scatter.png');
exportgraphics(fig2, lineScatterPng, 'Resolution', 300);
fprintf('Saved PNG only: %s\n', lineScatterPng);

close(fig2);

% 100% stacked composition plot across larger radial bins
fig3 = figure('Color', 'w', 'Position', [160 160 1150 650]);

rowTotals = sum(countMat, 2);
pctMat = zeros(size(countMat));
validRows = rowTotals > 0;
pctMat(validRows, :) = 100 * countMat(validRows, :) ./ rowTotals(validRows);

ringCatPct = categorical(ringLabels, ringLabels, 'Ordinal', true);
b3 = bar(ringCatPct, pctMat, 'stacked', 'EdgeColor', 'none');
for c = 1:numel(typeCats)
    b3(c).FaceColor = typeColors(c, :);
end

xlabel('Aggregated radial ring from tumor center');
ylabel('Composition (%)');
title(sprintf('Relative Cell-Type Composition by Larger Rings (%d rings per bin)', ringsPerGroup));
ylim([0 100]);
grid on;
legend(typeCats, 'Interpreter', 'none', 'Location', 'eastoutside');
set_publication_figure_style(fig3, 'font_size', 11, 'font_name', 'Arial');
add_figure_metadata('crc_tumor_caf_model.py:RING_DATA', totalCells, analysisDate);

compositionPng = fullfile(outputDir, 'ring_counts_stacked_percent.png');
exportgraphics(fig3, compositionPng, 'Resolution', 300);
fprintf('Saved PNG only: %s\n', compositionPng);

close(fig3);

% Cumulative radial burden by cell type (outward from center)
fig4 = figure('Color', 'w', 'Position', [180 180 1150 650]);
hold on;

cumulativeCounts = cumsum(countsByRing, 1);
lineHandlesCum = gobjects(numel(typeCats), 1);
for c = 1:numel(typeCats)
    yCum = cumulativeCounts(:, c);
    lineHandlesCum(c) = plot(ringCenterUm, yCum, '-', 'LineWidth', 2.2, ...
        'Color', typeColors(c, :), 'DisplayName', typeCats(c));
    scatter(ringCenterUm, yCum, 22, ...
        'MarkerFaceColor', typeColors(c, :), ...
        'MarkerEdgeColor', 'white', ...
        'LineWidth', 0.4, ...
        'MarkerFaceAlpha', 0.85, ...
        'HandleVisibility', 'off');
end

hold off;
xlabel('Ring center radius (\mum)');
ylabel('Cumulative cell count (0 \rightarrow r)');
title('Cumulative Radial Cell Burden by Type');
grid on;
legend(lineHandlesCum, typeCats, 'Interpreter', 'none', 'Location', 'eastoutside');
set_publication_figure_style(fig4, 'font_size', 11, 'font_name', 'Arial');
add_figure_metadata('crc_tumor_caf_model.py:RING_DATA', totalCells, analysisDate);

cumulativePng = fullfile(outputDir, 'ring_counts_cumulative_by_type.png');
exportgraphics(fig4, cumulativePng, 'Resolution', 300);
fprintf('Saved PNG only: %s\n', cumulativePng);

close(fig4);

% Immune pressure vs tumor burden by radius:
% ratio = CD8 / (CSC + CP + CD), computed per 50 um ring
fig5 = figure('Color', 'w', 'Position', [200 200 1150 650]);

idxCSC = 1;
idxCP = 2;
idxCD = 3;
idxCD8 = 4;

tumorBurden = countsByRing(:, idxCSC) + countsByRing(:, idxCP) + countsByRing(:, idxCD);
immuneCD8 = countsByRing(:, idxCD8);

ratioImmuneToTumor = zeros(size(tumorBurden));
valid = tumorBurden > 0;
ratioImmuneToTumor(valid) = immuneCD8(valid) ./ tumorBurden(valid);

yyaxis left;
plot(ringCenterUm, tumorBurden, '-', 'LineWidth', 2.4, 'Color', [0.20 0.20 0.20]);
hold on;
scatter(ringCenterUm, tumorBurden, 30, 'filled', ...
    'MarkerFaceColor', [0.20 0.20 0.20], 'MarkerEdgeColor', 'white', ...
    'LineWidth', 0.5, 'HandleVisibility', 'off');
ylabel('Tumor burden per ring (CSC+CP+CD)');

yyaxis right;
plot(ringCenterUm, ratioImmuneToTumor, '-', 'LineWidth', 2.4, 'Color', typeColors(idxCD8, :));
scatter(ringCenterUm, ratioImmuneToTumor, 30, 'filled', ...
    'MarkerFaceColor', typeColors(idxCD8, :), 'MarkerEdgeColor', 'white', ...
    'LineWidth', 0.5, 'HandleVisibility', 'off');
ylabel('Immune pressure ratio: CD8 / (CSC+CP+CD)');

xlabel('Ring center radius (\mum)');
title('Immune Pressure vs Tumor Burden by Radius');
grid on;
legend({'Tumor burden (CSC+CP+CD)', 'Immune pressure ratio (CD8/Tumor)'}, ...
    'Interpreter', 'none', 'Location', 'northwest');

set_publication_figure_style(fig5, 'font_size', 11, 'font_name', 'Arial');
add_figure_metadata('crc_tumor_caf_model.py:RING_DATA', totalCells, analysisDate);

immuneRatioPng = fullfile(outputDir, 'ring_immune_pressure_vs_tumor_burden.png');
exportgraphics(fig5, immuneRatioPng, 'Resolution', 300);
fprintf('Saved PNG only: %s\n', immuneRatioPng);

close(fig5);
