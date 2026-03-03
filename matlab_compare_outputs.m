% Compare two PhysiCell export directories (e.g., untreated vs treated)
% and generate side-by-side trajectory + final-state comparison figures.
%
% Usage options:
% 1) Edit DIR_A / DIR_B below, then run this script.
% 2) Or call from MATLAB command window:
%      compare_two_physicell_dirs("physicell_csvs_untreated", "physicell_csvs_treated")

DIR_A = "physicell_csvs_7";
DIR_B = "physicell_csvs_8";
LABEL_A = "Treated with Chronomodulated";
LABEL_B = "Treated With de Gramont";
OUT_SUBDIR = "comparison_figures_png2";
HEADLESS_PNG_BATCH_MODE = false;
DOCK_ALL_FIGURES = true;

if ~isempty(getenv('PHYSICELL_HEADLESS_PNG'))
    envFlag = string(getenv('PHYSICELL_HEADLESS_PNG'));
    HEADLESS_PNG_BATCH_MODE = any(strcmpi(envFlag, ["1","true","yes","on"]));
end

compare_two_physicell_dirs(DIR_A, DIR_B, LABEL_A, LABEL_B, OUT_SUBDIR, HEADLESS_PNG_BATCH_MODE, DOCK_ALL_FIGURES);


function compare_two_physicell_dirs(dirA, dirB, labelA, labelB, outSubdir, headlessMode, dockFigures)
    if nargin < 3 || strlength(string(labelA)) == 0
        labelA = "Condition_A";
    end
    if nargin < 4 || strlength(string(labelB)) == 0
        labelB = "Condition_B";
    end
    if nargin < 5 || strlength(string(outSubdir)) == 0
        outSubdir = "comparison_figures_png";
    end
    if nargin < 6
        headlessMode = false;
    end
    if nargin < 7
        dockFigures = true;
    end

    dirA = string(dirA);
    dirB = string(dirB);

    if ~isfolder(dirA)
        error("Directory A not found: %s", dirA);
    end
    if ~isfolder(dirB)
        error("Directory B not found: %s", dirB);
    end

    defaultFigureVisibleOriginal = get(groot, 'defaultFigureVisible');
    defaultFigureWindowStyleOriginal = get(groot, 'defaultFigureWindowStyle');
    if headlessMode
        set(groot, 'defaultFigureVisible', 'off');
    end
    if dockFigures && ~headlessMode
        set(groot, 'defaultFigureWindowStyle', 'docked');
    end

    cleanupObj = onCleanup(@() restore_figure_defaults(defaultFigureVisibleOriginal, defaultFigureWindowStyleOriginal));

    outDir = fullfile(dirB, outSubdir);
    if ~exist(outDir, "dir")
        mkdir(outDir);
    end

    A = load_condition_tables(dirA);
    B = load_condition_tables(dirB);

    totalFigures = 8;
    figureCount = 0;
    if headlessMode
        hwb = false;
    else
        hwb = waitbar(0, sprintf('Generating comparison figures: %d/%d', figureCount, totalFigures), ...
            'Name', 'PhysiCell Comparison Progress');
    end

    tA = A.summary.time_min;
    tB = B.summary.time_min;

    % 1) Total cells over time
    fig = figure("Name", "Total Cells Comparison");
    plot(tA, A.summary.total_cells, "-", "LineWidth", 2.0); hold on;
    plot(tB, B.summary.total_cells, "-", "LineWidth", 2.0); hold off;
    xlabel("Time (min)"); ylabel("Total cells");
    title("Total Cells Over Time");
    legend({labelA, labelB}, "Location", "best", "Interpreter", "none");
    grid on;
    exportPngPaper(fig, fullfile(outDir, "01_total_cells_comparison.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Figure %d/%d', figureCount, totalFigures), totalFigures);

    % 2) Live/dead over time (two-panel)
    fig = figure("Name", "Live Dead Comparison");
    tiledlayout(1,2, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    [liveA, deadA] = get_live_dead(A.summary);
    [liveB, deadB] = get_live_dead(B.summary);
    plot(tA, liveA, "-", "LineWidth", 1.8); hold on;
    plot(tB, liveB, "-", "LineWidth", 1.8); hold off;
    xlabel("Time (min)"); ylabel("Live cells"); title("Live Cells"); grid on;
    legend({labelA, labelB}, "Location", "best", "Interpreter", "none");

    nexttile;
    plot(tA, deadA, "-", "LineWidth", 1.8); hold on;
    plot(tB, deadB, "-", "LineWidth", 1.8); hold off;
    xlabel("Time (min)"); ylabel("Dead cells"); title("Dead Cells"); grid on;
    legend({labelA, labelB}, "Location", "best", "Interpreter", "none");

    exportPngPaper(fig, fullfile(outDir, "02_live_dead_comparison.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Figure %d/%d', figureCount, totalFigures), totalFigures);

    % 3) Net growth rate dN/dt
    fig = figure("Name", "Growth Rate Comparison");
    [tgA, gA] = finite_diff_rate(tA, A.summary.total_cells);
    [tgB, gB] = finite_diff_rate(tB, B.summary.total_cells);
    plot(tgA, gA, "-", "LineWidth", 1.8); hold on;
    plot(tgB, gB, "-", "LineWidth", 1.8); hold off;
    yline(0, "--", "LineWidth", 1.0);
    xlabel("Time (min)"); ylabel("dN/dt (cells/min)");
    title("Net Growth Rate Over Time");
    legend({labelA, labelB}, "Location", "best", "Interpreter", "none");
    grid on;
    exportPngPaper(fig, fullfile(outDir, "03_growth_rate_comparison.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Figure %d/%d', figureCount, totalFigures), totalFigures);

    % 4) Cell-type composition (%) over time for both conditions
    fig = figure("Name", "Cell Type Composition Percent Comparison");
    tiledlayout(1,2, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    [tTypeA, pctMatA, typeNamesA] = type_percent_matrix(A.typeLong);
    if ~isempty(tTypeA)
        area(tTypeA, pctMatA, "LineStyle", "none");
        ylim([0 100]);
        xlabel("Time (min)"); ylabel("Composition (%)");
        title(labelA + " - Cell Type Composition");
        grid on;
        legend(typeNamesA, "Interpreter", "none", "Location", "eastoutside");
    else
        text(0.5, 0.5, "No type data", "HorizontalAlignment", "center"); axis off;
    end

    nexttile;
    [tTypeB, pctMatB, typeNamesB] = type_percent_matrix(B.typeLong);
    if ~isempty(tTypeB)
        area(tTypeB, pctMatB, "LineStyle", "none");
        ylim([0 100]);
        xlabel("Time (min)"); ylabel("Composition (%)");
        title(labelB + " - Cell Type Composition");
        grid on;
        legend(typeNamesB, "Interpreter", "none", "Location", "eastoutside");
    else
        text(0.5, 0.5, "No type data", "HorizontalAlignment", "center"); axis off;
    end

    exportPngPaper(fig, fullfile(outDir, "04_cell_type_composition_percent_comparison.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Figure %d/%d', figureCount, totalFigures), totalFigures);

    % 5) Final-time cell-type counts side-by-side
    fig = figure("Name", "Final Cell Type Counts Comparison");
    [cats, finalA, finalB] = final_type_counts_pair(A.typeLong, B.typeLong);
    X = [finalA(:), finalB(:)];
    bar(categorical(cats, cats, 'Ordinal', true), X, 'grouped', 'EdgeColor', 'none');
    xlabel("Cell type"); ylabel("Cell count");
    title("Final Cell Type Counts");
    legend({labelA, labelB}, "Location", "best", "Interpreter", "none");
    grid on;
    exportPngPaper(fig, fullfile(outDir, "05_final_cell_type_counts_comparison.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Figure %d/%d', figureCount, totalFigures), totalFigures);

    % 6) Microenvironment means comparison (first common substrate)
    fig = figure("Name", "Microenvironment Mean Comparison");
    [subName, subA, subB] = first_common_substrate(A.micro, B.micro);
    if strlength(subName) > 0
        plot(subA.time_min, subA.mean, "-", "LineWidth", 1.8); hold on;
        plot(subB.time_min, subB.mean, "-", "LineWidth", 1.8); hold off;
        xlabel("Time (min)"); ylabel("Mean concentration");
        title("Microenvironment Mean: " + subName);
        legend({labelA, labelB}, "Location", "best", "Interpreter", "none");
        grid on;
    else
        text(0.5, 0.5, "No common substrate data", "HorizontalAlignment", "center");
        axis off;
    end
    exportPngPaper(fig, fullfile(outDir, "06_microenvironment_mean_comparison.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Figure %d/%d', figureCount, totalFigures), totalFigures);

    % 7) Ring totals over time
    fig = figure("Name", "Ring Totals Comparison");
    tiledlayout(1,2, "Padding", "compact", "TileSpacing", "compact");
    nexttile;
    plot_ring_totals(A.region, labelA);
    nexttile;
    plot_ring_totals(B.region, labelB);
    exportPngPaper(fig, fullfile(outDir, "07_ring_totals_comparison.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Figure %d/%d', figureCount, totalFigures), totalFigures);

    % 8) Final summary table export (as CSV + bar chart)
    finalTbl = build_final_summary_table(A.summary, B.summary, labelA, labelB);
    writetable(finalTbl, fullfile(outDir, "08_final_summary_table.csv"));

    fig = figure("Name", "Final Summary Delta");
    metrics = string(finalTbl.metric);
    deltaPct = finalTbl.delta_percent_B_vs_A;
    bar(categorical(metrics, metrics, 'Ordinal', true), deltaPct, 'EdgeColor', 'none');
    yline(0, "--", "LineWidth", 1.0);
    ylabel("Delta % (B vs A)"); xlabel("Metric");
    title("Final Metrics Percent Change (" + labelB + " vs " + labelA + ")");
    grid on;
    exportPngPaper(fig, fullfile(outDir, "08_final_summary_delta_percent.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Figure %d/%d', figureCount, totalFigures), totalFigures);

    if ~isequal(hwb, false) && ~isempty(hwb) && isgraphics(hwb)
        close(hwb);
    end

    clear cleanupObj;

    disp("Comparison complete.");
    disp("Figures written to: " + string(outDir));
end


function S = load_condition_tables(rootDir)
    S.summary = must_read_table(rootDir, "summary_by_time.csv");
    S.typeLong = must_read_table(rootDir, "cell_type_counts_long.csv");
    S.phaseLong = must_read_table(rootDir, "cell_phase_counts_long.csv");
    S.micro = must_read_table(rootDir, "microenvironment_stats.csv");
    S.region = must_read_table(rootDir, "region_counts.csv");
end


function T = must_read_table(rootDir, filename)
    f = fullfile(rootDir, filename);
    if ~isfile(f)
        error("Required file missing: %s", f);
    end
    T = readtable(f);
end


function code = parseCode(value)
    if isnumeric(value)
        if isempty(value)
            code = NaN;
        else
            code = round(double(value(1)));
        end
        return;
    end

    if islogical(value)
        code = double(value);
        return;
    end

    if iscell(value)
        if isempty(value)
            code = NaN;
        else
            code = parseCode(value{1});
        end
        return;
    end

    s = string(value);
    if strlength(s) == 0 || strcmpi(s, "nan") || strcmpi(s, "missing")
        code = NaN;
        return;
    end

    s = regexprep(s, '^(var|x)', '', 'ignorecase');
    numVal = str2double(s);
    if isnan(numVal)
        code = NaN;
    else
        code = round(numVal);
    end
end


function name = typeName(typeNum)
    % Cell types from initial.xml <cell_types>
    % 0 = Cancer stem, 1 = cancer proliferating, 2 = cancer differentiated,
    % 3 = CAF, 4 = CD8 T cell, 5 = M2 macrophage.
    typeMap = containers.Map(...
        {0, 1, 2, 3, 4, 5}, ...
        {"Cancer stem cell", "Cancer proliferating cell", "Cancer differentiated cell", "CAF", "CD8 T cell", "M2 macrophage"});

    key = parseCode(typeNum);
    if isnan(key)
        name = string(typeNum);
        return;
    end

    if isKey(typeMap, key)
        name = typeMap(key);
    else
        name = "Type_" + string(key);
    end
end


function labels = normalizeTypeLabels(col)
    s = string(col);
    labels = strings(size(s));
    for i = 1:numel(s)
        labels(i) = typeName(s(i));
    end
end


function [liveVals, deadVals] = get_live_dead(summaryTbl)
    if ismember("live_cells", summaryTbl.Properties.VariableNames)
        liveVals = summaryTbl.live_cells;
    elseif ismember("dead_cells", summaryTbl.Properties.VariableNames)
        liveVals = summaryTbl.total_cells - summaryTbl.dead_cells;
    else
        liveVals = summaryTbl.total_cells;
    end

    if ismember("dead_cells", summaryTbl.Properties.VariableNames)
        deadVals = summaryTbl.dead_cells;
    else
        deadVals = summaryTbl.total_cells - liveVals;
    end
    deadVals(deadVals < 0) = 0;
end


function [tMid, rate] = finite_diff_rate(t, y)
    t = double(t(:));
    y = double(y(:));
    if numel(t) < 2 || numel(y) < 2
        tMid = t;
        rate = zeros(size(t));
        return;
    end

    dt = diff(t);
    dy = diff(y);
    dt(dt == 0) = NaN;
    rate = dy ./ dt;
    tMid = t(1:end-1) + dt / 2;

    valid = isfinite(rate) & isfinite(tMid);
    tMid = tMid(valid);
    rate = rate(valid);
end


function [timeVals, pctMat, labels] = type_percent_matrix(typeLong)
    if isempty(typeLong) || height(typeLong) == 0
        timeVals = [];
        pctMat = [];
        labels = {};
        return;
    end

    timeVals = unique(typeLong.time_min, "sorted");
    rawTypes = normalizeTypeLabels(typeLong.cell_type);
    typeIDs = unique(rawTypes, "stable");

    counts = zeros(numel(timeVals), numel(typeIDs));
    for i = 1:numel(typeIDs)
        m = rawTypes == typeIDs(i);
        tLocal = typeLong.time_min(m);
        cLocal = typeLong.count(m);
        [ok, idx] = ismember(tLocal, timeVals);
        if any(ok)
            counts(:, i) = counts(:, i) + accumarray(idx(ok), cLocal(ok), [numel(timeVals), 1], @sum, 0);
        end
    end

    totals = sum(counts, 2);
    pctMat = zeros(size(counts));
    nz = totals > 0;
    pctMat(nz, :) = 100 * counts(nz, :) ./ totals(nz);

    labels = cellstr(typeIDs);
end


function [cats, finalA, finalB] = final_type_counts_pair(typeA, typeB)
    [catsA, valsA] = final_type_counts(typeA);
    [catsB, valsB] = final_type_counts(typeB);

    cats = unique([catsA; catsB], "stable");
    finalA = zeros(numel(cats), 1);
    finalB = zeros(numel(cats), 1);

    for i = 1:numel(cats)
        ia = find(catsA == cats(i), 1, 'first');
        ib = find(catsB == cats(i), 1, 'first');
        if ~isempty(ia), finalA(i) = valsA(ia); end
        if ~isempty(ib), finalB(i) = valsB(ib); end
    end
end


function [cats, vals] = final_type_counts(typeLong)
    if isempty(typeLong) || height(typeLong) == 0
        cats = strings(0,1);
        vals = zeros(0,1);
        return;
    end

    tfinal = max(typeLong.time_min);
    T = typeLong(typeLong.time_min == tfinal, :);
    raw = normalizeTypeLabels(T.cell_type);
    cats = unique(raw, "stable");
    vals = zeros(numel(cats), 1);
    for i = 1:numel(cats)
        vals(i) = sum(T.count(raw == cats(i)));
    end
end


function [name, subA, subB] = first_common_substrate(microA, microB)
    name = "";
    subA = table();
    subB = table();

    if isempty(microA) || isempty(microB) || height(microA) == 0 || height(microB) == 0
        return;
    end

    namesA = unique(string(microA.substrate), "stable");
    namesB = unique(string(microB.substrate), "stable");
    common = intersect(namesA, namesB, "stable");
    if isempty(common)
        return;
    end

    name = common(1);
    subA = microA(string(microA.substrate) == name, :);
    subB = microB(string(microB.substrate) == name, :);
end


function plot_ring_totals(regionTbl, labelText)
    if isempty(regionTbl) || height(regionTbl) == 0
        text(0.5, 0.5, "No region data", "HorizontalAlignment", "center");
        axis off;
        return;
    end

    tvals = unique(regionTbl.time_min, "sorted");
    regs = unique(string(regionTbl.region), "stable");
    M = zeros(numel(tvals), numel(regs));

    for i = 1:numel(regs)
        m = string(regionTbl.region) == regs(i);
        tLocal = regionTbl.time_min(m);
        cLocal = regionTbl.count(m);
        [ok, idx] = ismember(tLocal, tvals);
        if any(ok)
            M(:, i) = M(:, i) + accumarray(idx(ok), cLocal(ok), [numel(tvals), 1], @sum, 0);
        end
    end

    area(tvals, M, "LineStyle", "none");
    title(labelText + " - Ring Counts");
    xlabel("Time (min)");
    ylabel("Cell count");
    grid on;
    legend(cellstr(regs), "Interpreter", "none", "Location", "eastoutside");
end


function finalTbl = build_final_summary_table(sumA, sumB, labelA, labelB)
    idxA = find(sumA.time_min == max(sumA.time_min), 1, 'first');
    idxB = find(sumB.time_min == max(sumB.time_min), 1, 'first');

    [liveA, deadA] = get_live_dead(sumA);
    [liveB, deadB] = get_live_dead(sumB);

    metrics = ["total_cells"; "live_cells"; "dead_cells"];
    valsA = [sumA.total_cells(idxA); liveA(idxA); deadA(idxA)];
    valsB = [sumB.total_cells(idxB); liveB(idxB); deadB(idxB)];

    deltaAbs = valsB - valsA;
    deltaPct = zeros(size(deltaAbs));
    nz = valsA ~= 0;
    deltaPct(nz) = 100 * deltaAbs(nz) ./ valsA(nz);

    varA = make_valid_var_name(labelA, "condition_A");
    varB = make_valid_var_name(labelB, "condition_B");

    varNames = ["metric", string(varA), string(varB), "delta_abs_B_minus_A", "delta_percent_B_vs_A"];
    varNames = matlab.lang.makeUniqueStrings(varNames);

    finalTbl = table(metrics, valsA, valsB, deltaAbs, deltaPct, ...
        'VariableNames', cellstr(varNames));
end


function outName = make_valid_var_name(labelText, fallback)
    s = strtrim(char(string(labelText)));
    if isempty(s)
        s = fallback;
    end
    outName = matlab.lang.makeValidName(s, 'ReplacementStyle', 'delete');
    if isempty(outName)
        outName = fallback;
    end
end


function restore_figure_defaults(defaultVisible, defaultWindowStyle)
    try
        set(groot, 'defaultFigureVisible', defaultVisible);
    catch
    end
    try
        set(groot, 'defaultFigureWindowStyle', defaultWindowStyle);
    catch
    end
end


function h = safeWaitbarUpdate(progressValue, h, messageText, totalCount)
    progressValue = max(0, min(1, progressValue));
    if nargin < 4 || isempty(totalCount)
        totalCount = 0;
    end

    if isequal(h, false)
        return;
    end

    if ~isempty(h) && isgraphics(h)
        waitbar(progressValue, h, char(messageText));
    else
        if totalCount > 0
            h = waitbar(progressValue, char(messageText), 'Name', sprintf('PhysiCell Comparison Progress (%d figs)', totalCount));
        else
            h = waitbar(progressValue, char(messageText), 'Name', 'PhysiCell Comparison Progress');
        end
    end
end


function exportPngPaper(figHandle, outPath)
    figWidthCm = 17.5;
    figHeightCm = 11.0;
    printDpi = 600;
    baseFontSize = 9;

    try
        set(figHandle, 'Color', 'w');
        set(figHandle, 'Units', 'centimeters');
        p = get(figHandle, 'Position');
        set(figHandle, 'Position', [p(1), p(2), figWidthCm, figHeightCm]);
    catch
    end

    axAll = findall(figHandle, 'Type', 'axes');
    for i = 1:numel(axAll)
        try
            set(axAll(i), 'FontSize', baseFontSize, 'LineWidth', 1.0, 'Box', 'off');
            if ~isempty(axAll(i).Toolbar)
                axAll(i).Toolbar.Visible = 'off';
            end
        catch
            try
                axtoolbar(axAll(i), 'none');
            catch
            end
        end
    end

    lgAll = findall(figHandle, 'Type', 'legend');
    for i = 1:numel(lgAll)
        try
            set(lgAll(i), 'FontSize', baseFontSize);
        catch
        end
    end

    exportgraphics(figHandle, outPath, "Resolution", printDpi);

    try
        [figFolder, figName, ~] = fileparts(char(outPath));
        savefig(figHandle, fullfile(figFolder, figName + ".fig"));
    catch
    end

    try
        if strcmpi(string(get(figHandle, 'Visible')), "off")
            close(figHandle);
        end
    catch
    end
end
