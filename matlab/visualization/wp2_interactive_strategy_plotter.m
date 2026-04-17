function wp2_interactive_strategy_plotter(wp2Root)
% WP2 interactive plotting and table explorer.
% Scans strategy folders (WP2S1_CSV, WP2S2_CSV, ... or legacy WP2/1_Hourly,
% WP2/2_Square, ...), labels them as S1, S2, ..., and provides interactive
% plotting for:
%   1) Strategy base CSV (e.g., summary_by_time.csv)
%   2) curve_fitter_highlights.csv (if present)
%   3) curve_fitter_table.csv (if present)
%
% Usage:
%   wp2_interactive_strategy_plotter
%   wp2_interactive_strategy_plotter('C:\path\to\repo')
%   wp2_interactive_strategy_plotter('C:\path\to\repo\WP2')

if nargin < 1 || isempty(wp2Root)
    wp2Root = pwd;
end

if ~isfolder(wp2Root)
    error('Root folder not found: %s', wp2Root);
end

strategies = discover_strategies(wp2Root);
if isempty(strategies)
    error('No strategy folders found under: %s', wp2Root);
end

styles = build_style_map(numel(strategies));

state = struct();
state.wp2Root = wp2Root;
state.strategies = strategies;
state.styles = styles;
state.datasetKind = 'strategy_csv';

screen = get(groot, 'ScreenSize');
figW = min(1680, round(0.94 * screen(3)));
figH = min(920, round(0.90 * screen(4)));
figX = max(20, round((screen(3) - figW) / 2));
figY = max(20, round((screen(4) - figH) / 2));
fig = uifigure('Name', 'WP2 Strategy Interactive Explorer', 'Position', [figX figY figW figH]);
mainGrid = uigridlayout(fig, [1 2]);
mainGrid.ColumnWidth = {300, '1x'};
mainGrid.RowHeight = {'1x'};
mainGrid.Padding = [10 10 10 10];
mainGrid.ColumnSpacing = 10;

ctrl = uigridlayout(mainGrid, [12 1]);
ctrl.Layout.Row = 1;
ctrl.Layout.Column = 1;
ctrl.RowHeight = {22, 28, 22, 28, 22, 28, 22, 120, 28, 28, '1x', 30};
ctrl.Padding = [6 6 6 6];
ctrl.RowSpacing = 6;

lblData = uilabel(ctrl, 'Text', 'Data source');
lblData.Layout.Row = 1;
datasetDrop = uidropdown(ctrl, ...
    'Items', {'strategy_csv', 'curve_fitter_highlights', 'curve_fitter_table'}, ...
    'Value', 'strategy_csv');
datasetDrop.Layout.Row = 2;

lblX = uilabel(ctrl, 'Text', 'X variable');
lblX.Layout.Row = 3;
xDrop = uidropdown(ctrl, 'Items', {'(loading...)'}, 'Value', '(loading...)');
xDrop.Layout.Row = 4;

lblY = uilabel(ctrl, 'Text', 'Y variable');
lblY.Layout.Row = 5;
yDrop = uidropdown(ctrl, 'Items', {'(loading...)'}, 'Value', '(loading...)');
yDrop.Layout.Row = 6;

lblStrategies = uilabel(ctrl, 'Text', 'Strategies to plot');
lblStrategies.Layout.Row = 7;

strategyItems = strategy_item_labels(strategies);
strategyList = uilistbox(ctrl, ...
    'Items', strategyItems, ...
    'Value', strategyItems, ...
    'Multiselect', 'on');
strategyList.Layout.Row = 8;

plotBtn = uibutton(ctrl, 'Text', 'Plot Selected', 'ButtonPushedFcn', @on_plot);
plotBtn.Layout.Row = 9;

refreshBtn = uibutton(ctrl, 'Text', 'Refresh Columns', 'ButtonPushedFcn', @on_refresh_columns);
refreshBtn.Layout.Row = 10;

statusArea = uitextarea(ctrl, ...
    'Editable', 'off', ...
    'Value', {'Ready.'});
statusArea.Layout.Row = 11;

closeBtn = uibutton(ctrl, 'Text', 'Close', 'ButtonPushedFcn', @(~,~) close(fig));
closeBtn.Layout.Row = 12;

rightGrid = uigridlayout(mainGrid, [2 1]);
rightGrid.Layout.Row = 1;
rightGrid.Layout.Column = 2;
rightGrid.RowHeight = {58, '1x'};
rightGrid.RowSpacing = 8;

plotCtrlGrid = uigridlayout(rightGrid, [2 6]);
plotCtrlGrid.Layout.Row = 1;
plotCtrlGrid.Layout.Column = 1;
plotCtrlGrid.RowHeight = {22, 28};
plotCtrlGrid.ColumnWidth = {40, '1x', 40, '1x', 95, 95};
plotCtrlGrid.ColumnSpacing = 6;
plotCtrlGrid.Padding = [4 4 4 4];

titleLbl = uilabel(plotCtrlGrid, 'Text', 'Title');
titleLbl.Layout.Row = 1;
titleLbl.Layout.Column = 1;
titleEdit = uieditfield(plotCtrlGrid, 'text', 'Value', 'WP2 Strategy Comparison');
titleEdit.Layout.Row = 1;
titleEdit.Layout.Column = [2 5];
applyLabelsBtn = uibutton(plotCtrlGrid, 'Text', 'Apply', 'ButtonPushedFcn', @on_apply_labels);
applyLabelsBtn.Layout.Row = 1;
applyLabelsBtn.Layout.Column = 6;

xLbl = uilabel(plotCtrlGrid, 'Text', 'X');
xLbl.Layout.Row = 2;
xLbl.Layout.Column = 1;
xLabelEdit = uieditfield(plotCtrlGrid, 'text', 'Value', '');
xLabelEdit.Layout.Row = 2;
xLabelEdit.Layout.Column = 2;

yLbl = uilabel(plotCtrlGrid, 'Text', 'Y');
yLbl.Layout.Row = 2;
yLbl.Layout.Column = 3;
yLabelEdit = uieditfield(plotCtrlGrid, 'text', 'Value', '');
yLabelEdit.Layout.Row = 2;
yLabelEdit.Layout.Column = 4;

exportPngBtn = uibutton(plotCtrlGrid, 'Text', 'Export PNG', 'ButtonPushedFcn', @on_export_png);
exportPngBtn.Layout.Row = 2;
exportPngBtn.Layout.Column = 5;

exportFigBtn = uibutton(plotCtrlGrid, 'Text', 'Export FIG', 'ButtonPushedFcn', @on_export_fig);
exportFigBtn.Layout.Row = 2;
exportFigBtn.Layout.Column = 6;

ax = uiaxes(rightGrid);
ax.Layout.Row = 2;
ax.Layout.Column = 1;
grid(ax, 'on');
title(ax, 'WP2 Strategy Comparison');
xlabel(ax, 'X');
ylabel(ax, 'Y');

datasetDrop.ValueChangedFcn = @on_dataset_changed;
strategyList.ValueChangedFcn = @on_strategy_selection_changed;
xDrop.ValueChangedFcn = @on_axis_changed;
yDrop.ValueChangedFcn = @on_axis_changed;
titleEdit.ValueChangedFcn = @on_apply_labels;
xLabelEdit.ValueChangedFcn = @on_apply_labels;
yLabelEdit.ValueChangedFcn = @on_apply_labels;

on_refresh_columns();
on_plot();

    function on_dataset_changed(~, ~)
        state.datasetKind = datasetDrop.Value;
        on_refresh_columns();
        on_plot();
    end

    function on_strategy_selection_changed(~, ~)
        on_refresh_columns();
        on_plot();
    end

    function on_axis_changed(~, ~)
        on_plot();
    end

    function on_apply_labels(~, ~)
        on_plot();
    end

    function on_refresh_columns(~, ~)
        selected = selected_strategy_indices(strategyList, strategies);
        if isempty(selected)
            xDrop.Items = {'(none)'};
            xDrop.Value = '(none)';
            yDrop.Items = {'(none)'};
            yDrop.Value = '(none)';
            statusArea.Value = {'No strategies selected.'};
            return;
        end

        commonCols = {};
        unionCols = {};
        missingFiles = strings(0,1);

        for i = 1:numel(selected)
            s = strategies(selected(i));
            [tbl, srcPath, ok, msg] = load_dataset_table(s, state.datasetKind);
            if ~ok
                missingFiles(end+1,1) = string(sprintf('%s: %s', s.code, msg)); %#ok<AGROW>
                continue;
            end

            numCols = numeric_like_columns(tbl);
            unionCols = union(unionCols, numCols, 'stable');
            if isempty(commonCols)
                commonCols = numCols;
            else
                commonCols = intersect(commonCols, numCols, 'stable');
            end

            if isempty(srcPath)
                missingFiles(end+1,1) = string(sprintf('%s: no file path resolved', s.code)); %#ok<AGROW>
            end
        end

        if isempty(commonCols)
            commonCols = unionCols;
        end

        if isempty(commonCols)
            commonCols = {'(none)'};
        end

        prevX = xDrop.Value;
        prevY = yDrop.Value;
        xDrop.Items = commonCols;
        yDrop.Items = commonCols;

        if ismember(prevX, commonCols)
            xDrop.Value = prevX;
        else
            xDrop.Value = commonCols{1};
        end

        if ismember(prevY, commonCols)
            yDrop.Value = prevY;
        elseif numel(commonCols) >= 2
            yDrop.Value = commonCols{2};
        else
            yDrop.Value = commonCols{1};
        end

        lines = {
            sprintf('Dataset: %s', state.datasetKind), ...
            sprintf('Selected strategies: %d', numel(selected)), ...
            sprintf('Numeric columns available: %d', numel(commonCols))
            };

        if ~isempty(missingFiles)
            lines{end+1} = 'Missing/unavailable files:';
            for k = 1:numel(missingFiles)
                lines{end+1} = char(missingFiles(k)); %#ok<AGROW>
            end
        end

        statusArea.Value = lines;
    end

    function on_plot(~, ~)
        cla(ax);
        hold(ax, 'on');

        selected = selected_strategy_indices(strategyList, strategies);
        xName = xDrop.Value;
        yName = yDrop.Value;

        if isempty(selected)
            title(ax, 'No strategies selected');
            hold(ax, 'off');
            return;
        end

        if strcmp(xName, '(none)') || strcmp(yName, '(none)')
            title(ax, 'No plottable numeric columns found for selected data source');
            hold(ax, 'off');
            return;
        end

        plotCount = 0;
        skipped = strings(0,1);
        xAll = [];
        yAll = [];

        for i = 1:numel(selected)
            idx = selected(i);
            s = strategies(idx);
            [tbl, ~, ok, msg] = load_dataset_table(s, state.datasetKind);
            if ~ok
                skipped(end+1,1) = string(sprintf('%s skipped: %s', s.code, msg)); %#ok<AGROW>
                continue;
            end

            if ~ismember(xName, tbl.Properties.VariableNames) || ~ismember(yName, tbl.Properties.VariableNames)
                skipped(end+1,1) = string(sprintf('%s skipped: x/y not found', s.code)); %#ok<AGROW>
                continue;
            end

            x = to_numeric(tbl.(xName));
            y = to_numeric(tbl.(yName));
            valid = isfinite(x) & isfinite(y);
            if nnz(valid) < 2
                skipped(end+1,1) = string(sprintf('%s skipped: <2 valid points', s.code)); %#ok<AGROW>
                continue;
            end

            x = x(valid);
            y = y(valid);
            [x, order] = sort(x);
            y = y(order);

            st = state.styles(idx);
            plot(ax, x, y, ...
                'DisplayName', s.code, ...
                'Color', st.color, ...
                'LineStyle', st.lineStyle, ...
                'LineWidth', 2.0);
            plotCount = plotCount + 1;
            xAll = [xAll; x(:)]; %#ok<AGROW>
            yAll = [yAll; y(:)]; %#ok<AGROW>
        end

        customX = strtrim(xLabelEdit.Value);
        customY = strtrim(yLabelEdit.Value);
        customTitle = strtrim(titleEdit.Value);

        if isempty(customX)
            customX = xName;
        end
        if isempty(customY)
            customY = yName;
        end
        if isempty(customTitle)
            customTitle = sprintf('Dataset: %s', state.datasetKind);
        end

        [plotX, xInterp] = resolve_label_text(customX, xName);
        [plotY, yInterp] = resolve_label_text(customY, yName);
        [plotTitle, tInterp] = resolve_label_text(customTitle, sprintf('Dataset: %s', state.datasetKind));

        xlabel(ax, plotX, 'Interpreter', xInterp);
        ylabel(ax, plotY, 'Interpreter', yInterp);
        title(ax, plotTitle, 'Interpreter', tInterp);

        if plotCount > 0
            legend(ax, 'Location', 'best');
            apply_tight_limits(ax, xAll, yAll);
        end

        hold(ax, 'off');

        lines = {
            sprintf('Plotted %d strategy lines.', plotCount), ...
            sprintf('X: %s | Y: %s', xName, yName)
            };
        for k = 1:numel(skipped)
            lines{end+1} = char(skipped(k)); %#ok<AGROW>
        end
        statusArea.Value = lines;
    end

    function on_export_png(~, ~)
        [file, path] = uiputfile({'*.png'; '*.pdf'}, 'Export Plot', fullfile(state.wp2Root, 'wp2_strategy_plot.png'));
        if isequal(file, 0)
            return;
        end

        outPath = fullfile(path, file);
        try
            figOut = build_publication_figure();
            exportgraphics(figOut, outPath, 'Resolution', 600);
            close(figOut);
            statusArea.Value = {sprintf('Exported plot to: %s', outPath)};
        catch ME
            statusArea.Value = {sprintf('Export failed: %s', ME.message)};
        end
    end

    function on_export_fig(~, ~)
        [file, path] = uiputfile({'*.fig'}, 'Export MATLAB FIG', fullfile(state.wp2Root, 'wp2_strategy_plot.fig'));
        if isequal(file, 0)
            return;
        end

        outPath = fullfile(path, file);
        try
            figOut = build_publication_figure();
        catch ME
            statusArea.Value = {sprintf('FIG export failed: %s', ME.message)};
            return;
        end

        try
            savefig(figOut, outPath);
            statusArea.Value = {sprintf('Exported FIG to: %s', outPath)};
        catch ME
            statusArea.Value = {sprintf('FIG export failed: %s', ME.message)};
        end
        close(figOut);
    end

    function [figOut, axOut] = build_publication_figure()
        selected = selected_strategy_indices(strategyList, strategies);
        xName = xDrop.Value;
        yName = yDrop.Value;

        if isempty(selected) || strcmp(xName, '(none)') || strcmp(yName, '(none)')
            error('Nothing to export: select valid strategies and x/y columns first.');
        end

        % Fixed publication export size.
        figWidthIn = 14 / 2.54;
        figHeightIn = 7.4 / 2.54;

        figOut = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 figWidthIn figHeightIn]);
        figOut.PaperUnits = 'inches';
        figOut.PaperSize = [figWidthIn figHeightIn];
        figOut.PaperPosition = [0 0 figWidthIn figHeightIn];

        axOut = axes(figOut, 'Position', [0.10 0.12 0.86 0.83]);
        hold(axOut, 'on');

        plotted = 0;
        xAll = [];
        yAll = [];
        for i = 1:numel(selected)
            idx = selected(i);
            s = strategies(idx);
            [tbl, ~, ok, ~] = load_dataset_table(s, state.datasetKind);
            if ~ok
                continue;
            end
            if ~ismember(xName, tbl.Properties.VariableNames) || ~ismember(yName, tbl.Properties.VariableNames)
                continue;
            end

            x = to_numeric(tbl.(xName));
            y = to_numeric(tbl.(yName));
            valid = isfinite(x) & isfinite(y);
            if nnz(valid) < 2
                continue;
            end

            x = x(valid);
            y = y(valid);
            [x, order] = sort(x);
            y = y(order);

            st = state.styles(idx);
            plot(axOut, x, y, ...
                'DisplayName', s.code, ...
                'Color', st.color, ...
                'LineStyle', st.lineStyle, ...
                'LineWidth', 2.0);

            plotted = plotted + 1;
            xAll = [xAll; x(:)]; %#ok<AGROW>
            yAll = [yAll; y(:)]; %#ok<AGROW>
        end

        customX = strtrim(xLabelEdit.Value);
        customY = strtrim(yLabelEdit.Value);
        customTitle = strtrim(titleEdit.Value);
        if isempty(customX)
            customX = xName;
        end
        if isempty(customY)
            customY = yName;
        end
        if isempty(customTitle)
            customTitle = sprintf('Dataset: %s', state.datasetKind);
        end

        [plotX, xInterp] = resolve_label_text(customX, xName);
        [plotY, yInterp] = resolve_label_text(customY, yName);
        [plotTitle, tInterp] = resolve_label_text(customTitle, sprintf('Dataset: %s', state.datasetKind));

        grid(axOut, 'on');
        xlabel(axOut, plotX, 'Interpreter', xInterp);
        ylabel(axOut, plotY, 'Interpreter', yInterp);
        title(axOut, plotTitle, 'Interpreter', tInterp);
        if plotted > 0
            legend(axOut, 'Location', 'best');
            apply_tight_limits(axOut, xAll, yAll);
        end

        hold(axOut, 'off');
        fit_axes_inside_figure(axOut);
    end

    function fit_axes_inside_figure(axHandle)
        % Minimize whitespace while guaranteeing labels/titles/ticks stay inside figure.
        set(axHandle, 'Units', 'normalized');
        drawnow;

        ti = get(axHandle, 'TightInset');
        if isempty(ti) || numel(ti) ~= 4
            return;
        end

        % Small safety padding so export rasterization doesn't clip glyphs.
        pad = 0.012;

        left = max(ti(1) + pad, 0.04);
        bottom = max(ti(2) + pad, 0.04);
        right = max(ti(3) + pad, 0.04);
        top = max(ti(4) + pad, 0.04);

        width = max(1 - left - right, 0.1);
        height = max(1 - bottom - top, 0.1);
        set(axHandle, 'Position', [left, bottom, width, height]);

        % Re-evaluate once after repositioning; keep the tighter safe fit.
        drawnow;
        ti2 = get(axHandle, 'TightInset');
        left2 = max(ti2(1) + pad, left);
        bottom2 = max(ti2(2) + pad, bottom);
        right2 = max(ti2(3) + pad, right);
        top2 = max(ti2(4) + pad, top);
        width2 = max(1 - left2 - right2, 0.1);
        height2 = max(1 - bottom2 - top2, 0.1);
        set(axHandle, 'Position', [left2, bottom2, width2, height2]);
    end

end

function strategies = discover_strategies(wp2Root)
d = strategy_dirs_from_root(wp2Root);

strategies = struct('code', {}, 'name', {}, 'folder', {}, 'strategyCsv', {}, 'highlightsCsv', {}, 'curveTableCsv', {});
for i = 1:numel(d)
    folderName = d(i).name;

    tokS = regexp(folderName, '^WP2S(\d+)_CSV$', 'tokens', 'once');
    tokLegacy = regexp(folderName, '^(\d+)_', 'tokens', 'once');
    if ~isempty(tokS)
        n = str2double(tokS{1});
    elseif ~isempty(tokLegacy)
        n = str2double(tokLegacy{1});
    else
        continue;
    end

    code = sprintf('S%d', n);
    folderPath = fullfile(d(i).folder, folderName);

    if ~isempty(tokS)
        preferred = {
            'summary_by_time.csv', ...
            'time_aligned_long.csv', ...
            'metadata.csv', ...
            'cell_type_counts_wide.csv'
            };
        strategyCsv = '';
        for p = 1:numel(preferred)
            c = fullfile(folderPath, preferred{p});
            if isfile(c)
                strategyCsv = c;
                break;
            end
        end
        if isempty(strategyCsv)
            fallback = dir(fullfile(folderPath, '*.csv'));
            if ~isempty(fallback)
                strategyCsv = fullfile(fallback(1).folder, fallback(1).name);
            end
        end
    else
        strategyCsv = fullfile(folderPath, sprintf('%s.csv', folderName));
        if ~isfile(strategyCsv)
            fallback = dir(fullfile(folderPath, '*.csv'));
            if ~isempty(fallback)
                strategyCsv = fullfile(fallback(1).folder, fallback(1).name);
            else
                strategyCsv = '';
            end
        end
    end

    if isempty(strategyCsv)
        fallback = dir(fullfile(folderPath, '*.csv'));
        if ~isempty(fallback)
            strategyCsv = fullfile(fallback(1).folder, fallback(1).name);
        else
            strategyCsv = '';
        end
    end

    highlightsCsv = first_existing(fullfile(folderPath, 'curve_fitter_highlights.csv'), ...
                                   fullfile(folderPath, '*curve*fitter*highlights*.csv'));
    curveTableCsv = first_existing(fullfile(folderPath, 'curve_fitter_table.csv'), ...
                                   fullfile(folderPath, '*curve*fitter*table*.csv'));

    strategies(end+1) = struct( ...
        'code', code, ...
        'name', folderName, ...
        'folder', folderPath, ...
        'strategyCsv', strategyCsv, ...
        'highlightsCsv', highlightsCsv, ...
        'curveTableCsv', curveTableCsv); %#ok<AGROW>
end

if ~isempty(strategies)
    nums = cellfun(@(c) str2double(c(2:end)), {strategies.code});
    [~, order] = sort(nums);
    strategies = strategies(order);
end
end

function out = first_existing(exactPath, wildcardPath)
if isfile(exactPath)
    out = exactPath;
    return;
end

matches = dir(wildcardPath);
if isempty(matches)
    out = '';
else
    out = fullfile(matches(1).folder, matches(1).name);
end
end

function labels = strategy_item_labels(strategies)
labels = arrayfun(@(s) sprintf('%s - %s', s.code, s.name), strategies, 'UniformOutput', false);
end

function idx = selected_strategy_indices(selectorComponent, strategies)
labels = strategy_item_labels(strategies);
selectedLabels = selectorComponent.Value;
if ischar(selectedLabels) || isstring(selectedLabels)
    selectedLabels = cellstr(selectedLabels);
end

idx = [];
for i = 1:numel(selectedLabels)
    m = find(strcmp(labels, selectedLabels{i}), 1, 'first');
    if ~isempty(m)
        idx(end+1) = m; %#ok<AGROW>
    end
end
end

function styles = build_style_map(n)
% Fixed strategy style order requested by user:
% 1) purple, 2) yellow, 3) pink, 4) orange, 5) blue, 6) green,
% 7) dashed purple, 8) dashed yellow.
baseColors = [
    0.4940 0.1840 0.5560  % purple
    0.9290 0.6940 0.1250  % yellow
    0.8500 0.3250 0.6500  % pink
    0.8500 0.3250 0.0980  % orange
    0.0000 0.4470 0.7410  % blue
    0.4660 0.6740 0.1880  % green
    0.4940 0.1840 0.5560  % dashed purple
    0.9290 0.6940 0.1250  % dashed yellow
    ];
lineStyles = {'-', '-', '-', '-', '-', '-', '--', '--'};

styles = struct('color', {}, 'lineStyle', {});
for i = 1:n
    styleIdx = mod(i - 1, size(baseColors, 1)) + 1;
    c = baseColors(styleIdx, :);
    ls = lineStyles{styleIdx};
    styles(end+1) = struct('color', c, 'lineStyle', ls); %#ok<AGROW>
end
end

function [tbl, sourcePath, ok, msg] = load_dataset_table(strategy, datasetKind)
tbl = table();
sourcePath = '';
ok = false;
msg = '';

switch datasetKind
    case 'strategy_csv'
        sourcePath = strategy.strategyCsv;
    case 'curve_fitter_highlights'
        sourcePath = strategy.highlightsCsv;
    case 'curve_fitter_table'
        sourcePath = strategy.curveTableCsv;
    otherwise
        msg = sprintf('Unknown dataset kind: %s', datasetKind);
        return;
end

if isempty(sourcePath)
    msg = sprintf('No %s file in %s', datasetKind, strategy.name);
    return;
end

if ~isfile(sourcePath)
    msg = sprintf('File missing: %s', sourcePath);
    return;
end

try
    tbl = readtable(sourcePath, 'VariableNamingRule', 'preserve');
    ok = true;
catch ME
    msg = sprintf('Failed to read %s (%s)', sourcePath, ME.message);
end
end

function cols = numeric_like_columns(tbl)
cols = {};
if isempty(tbl) || width(tbl) == 0
    return;
end

vars = tbl.Properties.VariableNames;
for i = 1:numel(vars)
    v = tbl.(vars{i});
    x = to_numeric(v);
    if nnz(isfinite(x)) >= 2
        cols{end+1} = vars{i}; %#ok<AGROW>
    end
end
end

function x = to_numeric(v)
if isnumeric(v)
    x = double(v);
    return;
end

if islogical(v)
    x = double(v);
    return;
end

if isdatetime(v)
    x = datenum(v);
    return;
end

if isduration(v)
    x = seconds(v);
    return;
end

if iscellstr(v) || isstring(v) || ischar(v) || iscategorical(v)
    try
        x = str2double(string(v));
    catch
        x = nan(size(v));
    end
    x = double(x);
    return;
end

if iscell(v)
    try
        x = str2double(string(v));
    catch
        x = nan(size(v));
    end
    x = double(x);
    return;
end

x = nan(size(v));
x = double(x);
end

function [outText, interp] = resolve_label_text(customText, fallbackText)
customText = char(string(customText));
fallbackText = char(string(fallbackText));

if isempty(strtrim(customText))
    outText = fallbackText;
    interp = 'none';
    return;
end

outText = customText;

% Allow TeX-like scientific labels, e.g. "Tumor Envelope Volume (\mu m^3)".
if contains(customText, '\\') || contains(customText, '^') || contains(customText, '{') || contains(customText, '}')
    interp = 'tex';
else
    interp = 'none';
end
end

function apply_tight_limits(axHandle, xAll, yAll)
if isempty(xAll) || isempty(yAll)
    return;
end

xAll = xAll(isfinite(xAll));
yAll = yAll(isfinite(yAll));
if isempty(xAll) || isempty(yAll)
    return;
end

xMin = min(xAll);
xMax = max(xAll);
yMin = min(yAll);
yMax = max(yAll);

if xMax == xMin
    xPad = max(1, abs(xMin) * 0.02);
else
    xPad = 0;
end

if yMax == yMin
    yPad = max(1, abs(yMin) * 0.02);
else
    yPad = 0;
end

xlim(axHandle, [xMin - xPad, xMax + xPad]);
ylim(axHandle, [yMin - yPad, yMax + yPad]);
end

function d = strategy_dirs_from_root(wp2Root)
% Prefer root-level WP2S*_CSV folders when present.
d = dir(fullfile(wp2Root, 'WP2S*_CSV'));
d = d([d.isdir]);
if ~isempty(d)
    return;
end

% If repo root is passed, support legacy layout in ./WP2.
legacyRoot = fullfile(wp2Root, 'WP2');
if isfolder(legacyRoot)
    d = dir(legacyRoot);
    d = d([d.isdir]);
    d = d(~ismember({d.name}, {'.', '..'}));
    return;
end

% Otherwise treat provided folder itself as the strategy container.
d = dir(wp2Root);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.', '..'}));
end