% MATLAB analysis helpers for PhysiCell CSV exports
% Update outDir to the folder created by export_outputs_to_excel.py

outDir = "Sim_results_8_csv"; %CHANGE ME TODO run as arg

% ===== OUTPUT / DISPLAY MODES =====
% When true: save all PNGs without showing figure windows (batch export mode).
HEADLESS_PNG_BATCH_MODE = false;

% Optional override from environment variable PHYSICELL_HEADLESS_PNG = 1/true/yes/on
headlessEnv = string(getenv('PHYSICELL_HEADLESS_PNG'));
if strlength(headlessEnv) > 0
    HEADLESS_PNG_BATCH_MODE = any(strcmpi(headlessEnv, ["1","true","yes","on"]));
end

% Preserve and optionally override default figure visibility.
defaultFigureVisibleOriginal = get(groot, 'defaultFigureVisible');
if HEADLESS_PNG_BATCH_MODE
    set(groot, 'defaultFigureVisible', 'off');
end

% ===== PROGRESS TRACKING =====
figureCount = 0;
totalFigures = 9;  % Total expected figures
if HEADLESS_PNG_BATCH_MODE
    hwb = false; % sentinel to disable waitbar creation/updates in headless mode
else
    hwb = waitbar(0, sprintf('Generating figures: 0/%d', totalFigures), 'Name', 'PhysiCell Analysis Progress');
end

% Helper function to extract display name from column name (types or phases)
function displayName = extractDisplayName(colName, isType)
    numVal = parseCode(colName);
    if ~isnan(numVal)
        if isType
            displayName = typeName(numVal);
        else
            displayName = phaseName(numVal);
        end
    else
        displayName = string(colName);
    end
end

% Robustly parse encoded IDs such as: 0, 0.0, "0", "var0", "x0"
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

% Convert table columns (numeric/string/cell/categorical) into numeric IDs
function codes = columnToCodes(col)
    if isnumeric(col)
        codes = round(double(col));
        return;
    end

    if islogical(col)
        codes = double(col);
        return;
    end

    if iscell(col)
        codes = NaN(size(col));
        for i = 1:numel(col)
            codes(i) = parseCode(col{i});
        end
        return;
    end

    s = string(col);
    codes = NaN(size(s));
    for i = 1:numel(s)
        codes(i) = parseCode(s(i));
    end
end

% Helper function to map phase numbers to names
function name = phaseName(phaseNum)
    phaseMap = containers.Map(...
        {0, 1, 2, 3, 4, 5, 14, 100, 101, 102, 103}, ...
        {"G0/G1", "S", "G2", "M", "Apoptotic", "Necrotic_Swelling", "Necrotic_Lysed", "Ki67_G0", "Ki67_G1", "Ki67_S", "Ki67_G2"});

    key = parseCode(phaseNum);
    if isnan(key)
        name = "Unknown_Phase";
        return;
    end

    if isKey(phaseMap, key)
        name = phaseMap(key);
    else
        name = "Phase_" + string(key);
    end
end

% Helper function to map type numbers to names (customize based on your model)
function name = typeName(typeNum)
    % Cell types from initial.xml <cell_types>
    typeMap = containers.Map(...
        {0, 1, 2, 3, 4, 5}, ...
        {"Cancer_stem", "cancer_proliferating", "cancer_differentiated", "CAF", "CD8_T_cell", "M2_macrophage"});

    key = parseCode(typeNum);
    if isnan(key)
        name = "Unknown_Type";
        return;
    end

    if isKey(typeMap, key)
        name = typeMap(key);
    else
        name = "Type_" + string(key);
    end
end

% Normalize cell_type column values into readable type labels
function labels = normalizeTypeLabels(col)
    if isnumeric(col)
        labels = strings(size(col));
        for i = 1:numel(col)
            labels(i) = typeName(col(i));
        end
        return;
    end

    s = string(col);
    labels = strings(size(s));
    for i = 1:numel(s)
        c = parseCode(s(i));
        if ~isnan(c)
            labels(i) = typeName(c);
        else
            labels(i) = s(i);
        end
    end
end

% Classification helpers for derived tumor/stromal summaries
function tf = isCancerLike(typeLabel)
    t = lower(string(typeLabel));
    tf = contains(t, "cancer") | contains(t, "stem") | contains(t, "prolif") | contains(t, "differ");
end

function tf = isStromaLike(typeLabel)
    t = lower(string(typeLabel));
    tf = contains(t, "caf") | contains(t, "fibro") | contains(t, "m2") | contains(t, "macro");
end

% Stable radial ordering for raw region labels (when labels are qualitative).
function orderedRaw = orderRawRegions(rawVals)
    raw = unique(string(rawVals), "stable");
    if isempty(raw)
        orderedRaw = strings(0,1);
        return;
    end

    preferred = ["core", "inner", "outer", "rim"];
    orderedRaw = strings(0,1);
    rawLower = lower(raw);
    for i = 1:numel(preferred)
        hit = find(rawLower == preferred(i), 1, 'first');
        if ~isempty(hit)
            orderedRaw(end+1,1) = raw(hit); %#ok<AGROW>
        end
    end
    for i = 1:numel(raw)
        if ~any(orderedRaw == raw(i))
            orderedRaw(end+1,1) = raw(i); %#ok<AGROW>
        end
    end
end

% Build fixed 250 um bins from 0 to 2000 um.
% Example: 0-250, 250-500, ..., 1750-2000 um.
function binLabels = buildDimRingLabels()
    stepUm = 250;
    maxUm = 2000;
    edges = (0:stepUm:maxUm)';
    nBins = numel(edges) - 1;

    binLabels = strings(nBins,1);
    for i = 1:nBins
        binLabels(i) = sprintf("%.0f-%.0f um", edges(i), edges(i+1));
    end
end

% Map raw region labels onto dimensioned bins using radial order.
function outLabels = mapRegionsToDimLabels(regionCol, ~)
    raw = string(regionCol);
    outLabels = raw;
    if isempty(raw)
        return;
    end

    orderedRaw = orderRawRegions(raw);
    binLabels = buildDimRingLabels(); %#ok<NASGU>
    nBins = numel(binLabels);

    for i = 1:numel(raw)
        idxRaw = find(orderedRaw == raw(i), 1, 'first');
        if isempty(idxRaw)
            continue;
        end
        idxBin = min(idxRaw, nBins);
        outLabels(i) = binLabels(idxBin);
    end
end

% Export helper to avoid axes toolbar artifacts in output images
function exportPngNoToolbar(figHandle, outPath)
    % Publication-oriented defaults (A4-compatible figure placement)
    % Full-width panel target for A4 manuscripts: ~175 mm wide.
    figWidthCm = 17.5;
    figHeightCm = 11.0;
    printDpi = 600;
    baseFontSize = 9;
    axisLineWidth = 1.0;

    try
        set(figHandle, 'Color', 'w');
        set(figHandle, 'Units', 'centimeters');
        pos = get(figHandle, 'Position');
        set(figHandle, 'Position', [pos(1), pos(2), figWidthCm, figHeightCm]);
    catch
    end

    % Harmonize axis/text styling for print readability.
    axAll = findall(figHandle, 'Type', 'axes');
    for i = 1:numel(axAll)
        try
            set(axAll(i), 'FontSize', baseFontSize, 'LineWidth', axisLineWidth, 'Box', 'off');
        catch
        end
    end

    lgAll = findall(figHandle, 'Type', 'legend');
    for i = 1:numel(lgAll)
        try
            set(lgAll(i), 'FontSize', baseFontSize);
        catch
        end
    end

    % Axis labels and titles slightly larger than base axis tick text.
    txtAll = findall(figHandle, 'Type', 'text');
    for i = 1:numel(txtAll)
        try
            currentSize = get(txtAll(i), 'FontSize');
            if isempty(currentSize) || ~isnumeric(currentSize)
                set(txtAll(i), 'FontSize', baseFontSize + 1);
            else
                set(txtAll(i), 'FontSize', max(currentSize, baseFontSize + 1));
            end
        catch
        end
    end

    axAll = findall(figHandle, 'Type', 'axes');
    for i = 1:numel(axAll)
        try
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
    exportgraphics(figHandle, outPath, "Resolution", printDpi);

    try
        [figFolder, figName, ~] = fileparts(char(outPath));
        savefig(figHandle, fullfile(figFolder, figName + ".fig"));
    catch
    end

    % In headless mode, figures are invisible; close them after export.
    try
        if strcmpi(string(get(figHandle, 'Visible')), "off")
            close(figHandle);
        end
    catch
    end
end

% Robust waitbar updater: reuses handle when valid, recreates if needed
function h = safeWaitbarUpdate(progressValue, h, messageText, totalCount)
    if isequal(h, false)
        return;
    end

    progressValue = max(0, min(1, progressValue));
    if nargin < 4 || isempty(totalCount)
        totalCount = 0;
    end

    if ~isempty(h) && isgraphics(h)
        waitbar(progressValue, h, char(messageText));
    else
        if totalCount > 0
            h = waitbar(progressValue, char(messageText), 'Name', sprintf('PhysiCell Analysis Progress (%d figs)', totalCount));
        else
            h = waitbar(progressValue, char(messageText), 'Name', 'PhysiCell Analysis Progress');
        end
    end
end

summary = readtable(fullfile(outDir, "summary_by_time.csv"));
cellTypeLong = readtable(fullfile(outDir, "cell_type_counts_long.csv"));
cellPhaseLong = readtable(fullfile(outDir, "cell_phase_counts_long.csv"));
cellTypeWide = readtable(fullfile(outDir, "cell_type_counts_wide.csv"));
cellPhaseWide = readtable(fullfile(outDir, "cell_phase_counts_wide.csv"));
cellTypePhase = readtable(fullfile(outDir, "cell_type_phase_counts.csv"));
spatialStats = readtable(fullfile(outDir, "spatial_stats_by_type.csv"));
regionCounts = readtable(fullfile(outDir, "region_counts.csv"));
attrStats = readtable(fullfile(outDir, "cell_attribute_stats.csv"));
microStats = readtable(fullfile(outDir, "microenvironment_stats.csv"));
metadata = readtable(fullfile(outDir, "metadata.csv"));

rMaxForRegionLabels = NaN;
if ~isempty(spatialStats) && ismember("radius_max", spatialStats.Properties.VariableNames)
    rMaxForRegionLabels = max(spatialStats.radius_max, [], "omitnan");
end

% base (9) + additional outDir-derived spatial proxy plots (4)
totalFigures = totalFigures + 4;
hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d', figureCount, totalFigures), totalFigures);

% Figure export setup
figOutDir = fullfile(outDir, "analysis_figures_png");
if ~exist(figOutDir, "dir")
    mkdir(figOutDir);
end

% Total cells over time
figure("Name", "Total Cells");
plot(summary.time_min, summary.total_cells, "-o", "LineWidth", 1.2);
xlabel("Time (min)");
ylabel("Total cells");
title("Total Cells Over Time");
grid on;
exportPngNoToolbar(gcf, fullfile(figOutDir, "01_total_cells.png"));

% Update progress
figureCount = figureCount + 1;
hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Total Cells', figureCount, totalFigures), totalFigures);

% Line plot of cell types over time (build from long table to avoid wide-header name issues)
if ~isempty(cellTypeLong) && height(cellTypeLong) > 0
    figure("Name", "Cell Types Over Time");

    timeVals = unique(cellTypeLong.time_min, "sorted");

    if isnumeric(cellTypeLong.cell_type)
        observedTypeIDs = unique(cellTypeLong.cell_type(~isnan(cellTypeLong.cell_type)), "stable");
    else
        observedTypeIDs = unique(str2double(string(cellTypeLong.cell_type)), "stable");
        observedTypeIDs = observedTypeIDs(~isnan(observedTypeIDs));
    end

    % Ensure key tumor types are always present in the plot
    % 0 = Cancer_stem, 1 = cancer_proliferating
    typeIDs = unique([0; 1; observedTypeIDs(:)], "stable");

    countsMat = zeros(numel(timeVals), numel(typeIDs));
    for iType = 1:numel(typeIDs)
        id = typeIDs(iType);
        if isnumeric(cellTypeLong.cell_type)
            mask = (cellTypeLong.cell_type == id);
        else
            mask = (str2double(string(cellTypeLong.cell_type)) == id);
        end

        if any(mask)
            tLocal = cellTypeLong.time_min(mask);
            cLocal = cellTypeLong.count(mask);
            [isMember, tIdx] = ismember(tLocal, timeVals);
            validRows = tIdx(isMember);
            countsMat(validRows, iType) = countsMat(validRows, iType) + cLocal(isMember);
        end
    end

    hold on;
    for iType = 1:numel(typeIDs)
        lw = 1.5;
        if typeIDs(iType) == 0 || typeIDs(iType) == 1
            lw = 2.2; % emphasize stem/proliferating
        end
        plot(timeVals, countsMat(:, iType), "LineWidth", lw);
    end
    hold off;

    displayNames = arrayfun(@(id) typeName(id), typeIDs, 'UniformOutput', false);

    xlabel("Time (min)");
    ylabel("Cell count");
    title("Cell Types Over Time");
    legend(displayNames, "Interpreter", "none", "Location", "bestoutside");
    grid on;
    exportPngNoToolbar(gcf, fullfile(figOutDir, "02_cell_types_over_time.png"));

    % Update progress
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Cell Types', figureCount, totalFigures), totalFigures);
end

% Line plot of cell phases over time
if ~isempty(cellPhaseLong) && height(cellPhaseLong) > 0
    figure("Name", "Cell Phases Over Time");

    timeVals = unique(cellPhaseLong.time_min, "sorted");
    phaseIDsRaw = columnToCodes(cellPhaseLong.phase);
    phaseIDs = unique(phaseIDsRaw(~isnan(phaseIDsRaw)), "stable");

    countsMat = zeros(numel(timeVals), numel(phaseIDs));
    for iPhase = 1:numel(phaseIDs)
        id = phaseIDs(iPhase);
        mask = (phaseIDsRaw == id);

        if any(mask)
            tLocal = cellPhaseLong.time_min(mask);
            cLocal = cellPhaseLong.count(mask);
            [isMember, tIdx] = ismember(tLocal, timeVals);
            validRows = tIdx(isMember);
            countsMat(validRows, iPhase) = countsMat(validRows, iPhase) + cLocal(isMember);
        end
    end

    hold on;
    for iPhase = 1:numel(phaseIDs)
        plot(timeVals, countsMat(:, iPhase), "LineWidth", 1.5);
    end
    hold off;

    displayNames = arrayfun(@(id) phaseName(id), phaseIDs, 'UniformOutput', false);

    xlabel("Time (min)");
    ylabel("Cell count");
    title("Cell Phases Over Time");
    legend(displayNames, "Interpreter", "none", "Location", "bestoutside");
    grid on;
    exportPngNoToolbar(gcf, fullfile(figOutDir, "03_cell_phases_over_time.png"));
    
    % Update progress
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Cell Phases', figureCount, totalFigures), totalFigures);
end

% Cell-type composition over time (stacked %)
if ~isempty(cellTypeLong) && height(cellTypeLong) > 0
    timeVals = unique(cellTypeLong.time_min, "sorted");
    typeIDsRaw = columnToCodes(cellTypeLong.cell_type);
    valid = ~isnan(typeIDsRaw);
    typeIDs = unique(typeIDsRaw(valid), "stable");

    countsMat = zeros(numel(timeVals), numel(typeIDs));
    for iType = 1:numel(typeIDs)
        id = typeIDs(iType);
        mask = (typeIDsRaw == id);
        if any(mask)
            tLocal = cellTypeLong.time_min(mask);
            cLocal = cellTypeLong.count(mask);
            [isMember, tIdx] = ismember(tLocal, timeVals);
            validRows = tIdx(isMember);
            countsMat(validRows, iType) = countsMat(validRows, iType) + cLocal(isMember);
        end
    end

    totals = sum(countsMat, 2);
    pctMat = zeros(size(countsMat));
    nz = totals > 0;
    pctMat(nz, :) = 100 * countsMat(nz, :) ./ totals(nz);

    figure("Name", "Cell Type Composition (%)");
    area(timeVals, pctMat, "LineStyle", "none");
    xlabel("Time (min)");
    ylabel("Composition (%)");
    title("Cell Type Composition Over Time");
    ylim([0 100]);
    grid on;
    typeLabels = arrayfun(@(id) typeName(id), typeIDs, 'UniformOutput', false);
    legend(typeLabels, "Interpreter", "none", "Location", "eastoutside");
    exportPngNoToolbar(gcf, fullfile(figOutDir, "04_cell_type_composition_percent.png"));

    % Update progress
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Composition', figureCount, totalFigures), totalFigures);
end

% Cell-type composition over time (stacked absolute counts)
if ~isempty(cellTypeLong) && height(cellTypeLong) > 0
    timeVals = unique(cellTypeLong.time_min, "sorted");
    typeIDsRaw = columnToCodes(cellTypeLong.cell_type);
    valid = ~isnan(typeIDsRaw);
    typeIDs = unique(typeIDsRaw(valid), "stable");

    countsMat = zeros(numel(timeVals), numel(typeIDs));
    for iType = 1:numel(typeIDs)
        id = typeIDs(iType);
        mask = (typeIDsRaw == id);
        if any(mask)
            tLocal = cellTypeLong.time_min(mask);
            cLocal = cellTypeLong.count(mask);
            [isMember, tIdx] = ismember(tLocal, timeVals);
            validRows = tIdx(isMember);
            countsMat(validRows, iType) = countsMat(validRows, iType) + cLocal(isMember);
        end
    end

    figure("Name", "Cell Type Composition (Counts)");
    area(timeVals, countsMat, "LineStyle", "none");
    xlabel("Time (min)");
    ylabel("Cell count");
    title("Cell Type Composition Over Time (Absolute Counts)");
    grid on;
    typeLabels = arrayfun(@(id) typeName(id), typeIDs, 'UniformOutput', false);
    legend(typeLabels, "Interpreter", "none", "Location", "eastoutside");
    exportPngNoToolbar(gcf, fullfile(figOutDir, "09_cell_type_composition_counts.png"));

    % Update progress
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Composition Counts', figureCount, totalFigures), totalFigures);
end

% Microenvironment substrate trends (mean over time)
if ~isempty(microStats) && height(microStats) > 0
    substrates = unique(microStats.substrate, "stable");
    figure("Name", "Microenvironment Means");
    hold on;
    
    % Vectorized approach: group by substrate
    for i = 1:numel(substrates)
        if iscell(substrates)
            mask = strcmp(string(microStats.substrate), string(substrates{i}));
        else
            mask = strcmp(microStats.substrate, substrates{i});
        end
        sub_data = microStats(mask, :);
        if height(sub_data) > 0
            plot(sub_data.time_min, sub_data.mean, "LineWidth", 1.1);
        end
    end
    hold off;
    xlabel("Time (min)");
    ylabel("Mean concentration");
    title("Microenvironment Means Over Time");
    legend(substrates, "Interpreter", "none", "Location", "best");
    grid on;
    exportPngNoToolbar(gcf, fullfile(figOutDir, "05_microenvironment_means.png"));
    
    % Update progress
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Microenvironment', figureCount, totalFigures), totalFigures);
end

% Spatial radius trends for a selected cell type
if ~isempty(spatialStats) && height(spatialStats) > 0
    selectedType = spatialStats.cell_type(1);
    sub = spatialStats(spatialStats.cell_type == selectedType, :);
    if ismember("radius_mean", sub.Properties.VariableNames)
        figure("Name", "Radius Mean Over Time");
        plot(sub.time_min, sub.radius_mean, "-o", "LineWidth", 1.2);
        xlabel("Time (min)");
        ylabel("Mean radius");
        title("Mean Radius Over Time");
        grid on;
        exportPngNoToolbar(gcf, fullfile(figOutDir, "06_radius_mean_over_time.png"));
        
        % Update progress
        figureCount = figureCount + 1;
        hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Spatial Radius', figureCount, totalFigures), totalFigures);
    end
end

% Viability burden over time (live vs dead)
if ~isempty(summary) && height(summary) > 0
    if ismember("live_cells", summary.Properties.VariableNames)
        liveVals = summary.live_cells;
    elseif ismember("dead_cells", summary.Properties.VariableNames)
        liveVals = summary.total_cells - summary.dead_cells;
    else
        liveVals = summary.total_cells;
    end

    if ismember("dead_cells", summary.Properties.VariableNames)
        deadVals = summary.dead_cells;
    else
        deadVals = summary.total_cells - liveVals;
    end

    deadVals(deadVals < 0) = 0;

    figure("Name", "Live vs Dead Cells Over Time");
    area(summary.time_min, [liveVals deadVals], "LineStyle", "none");
    xlabel("Time (min)");
    ylabel("Cell count");
    title("Tumor Burden and Viability Over Time");
    legend({"Live cells", "Dead cells"}, "Location", "best", "Interpreter", "none");
    grid on;
    exportPngNoToolbar(gcf, fullfile(figOutDir, "07_live_dead_burden_over_time.png"));

    % Update progress
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Viability', figureCount, totalFigures), totalFigures);
end

% Region counts over time for a selected type
if ~isempty(regionCounts) && height(regionCounts) > 0
    selectedType = regionCounts.cell_type(1);
    sub = regionCounts(regionCounts.cell_type == selectedType, :);
    rawRegionsOrdered = orderRawRegions(sub.region);
    regionDisplay = buildDimRingLabels();
    nBins = numel(regionDisplay);

    tVals = unique(sub.time_min, "sorted");
    ringMat = zeros(numel(tVals), nBins);
    rawSub = string(sub.region);
    [~, rawIdx] = ismember(rawSub, rawRegionsOrdered);
    ringIdx = min(max(rawIdx, 1), nBins);
    [~, tIdx] = ismember(sub.time_min, tVals);
    v = tIdx > 0 & ringIdx > 0;
    if any(v)
        linIdx = sub2ind(size(ringMat), tIdx(v), ringIdx(v));
        ringMat = ringMat + reshape(accumarray(linIdx, sub.count(v), [numel(ringMat), 1], @sum, 0), size(ringMat));
    end

    figure("Name", "Ring Counts");
    hold on;
    for i = 1:nBins
        plot(tVals, ringMat(:, i), "LineWidth", 1.1);
    end
    hold off;
    xlabel("Time (min)");
    ylabel("Cell count");
    title("Ring Counts Over Time");
    legend(cellstr(regionDisplay), "Interpreter", "none", "Location", "best");
    grid on;
    exportPngNoToolbar(gcf, fullfile(figOutDir, "08_region_counts_over_time.png"));
    
    % Update progress
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Ring Counts', figureCount, totalFigures), totalFigures);
end

% ====================== Additional outDir-only spatial proxy plots ======================
if ~isempty(regionCounts) && height(regionCounts) > 0
    RC = regionCounts;
    RC.region_raw = string(RC.region);
    RC.type_label = normalizeTypeLabels(RC.cell_type);

    rawRegionsOrdered = orderRawRegions(RC.region_raw);
    ringLabels = buildDimRingLabels();
    nRings = numel(ringLabels);

    [~, rawIdxAll] = ismember(RC.region_raw, rawRegionsOrdered);
    RC.ring_idx = min(max(rawIdxAll, 1), nRings);
    RC.region_label = ringLabels(RC.ring_idx);

    ringTotalsAll = accumarray(RC.ring_idx, RC.count, [nRings, 1], @sum, 0);
    activeRingIdx = find(ringTotalsAll > 0);
    if isempty(activeRingIdx)
        activeRingIdx = (1:nRings)';
    end

    % Figure 10: Region composition over time (all cells)
    timeVals = unique(RC.time_min, "sorted");
    regions = ringLabels(activeRingIdx);
    regionMat = zeros(numel(timeVals), numel(regions));
    for iR = 1:numel(activeRingIdx)
        m = RC.ring_idx == activeRingIdx(iR);
        tLocal = RC.time_min(m);
        cLocal = RC.count(m);
        [isMember, tIdx] = ismember(tLocal, timeVals);
        vRows = tIdx(isMember);
        cValid = cLocal(isMember);
        if ~isempty(vRows)
            regionMat(:, iR) = regionMat(:, iR) + accumarray(vRows, cValid, [numel(timeVals), 1], @sum, 0);
        end
    end

    figure("Name", "Ring Composition Over Time");
    area(timeVals, regionMat, "LineStyle", "none");
    xlabel("Time (min)");
    ylabel("Cell count");
    title("Ring Composition Over Time (All Cells)");
    legend(cellstr(regions), "Interpreter", "none", "Location", "eastoutside");
    grid on;
    exportPngNoToolbar(gcf, fullfile(figOutDir, "10_region_composition_over_time.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Ring Composition', figureCount, totalFigures), totalFigures);

    % Figure 11: Final-time stacked bar by region (cell-type amounts)
    tFinal = max(RC.time_min);
    RF = RC(RC.time_min == tFinal, :);
    typeList = unique(RF.type_label, "stable");
    regList = ringLabels(activeRingIdx);

    H = zeros(numel(regList), numel(typeList));
    for iT = 1:numel(typeList)
        for iR = 1:numel(activeRingIdx)
            m = (RF.type_label == typeList(iT)) & (RF.ring_idx == activeRingIdx(iR));
            H(iR, iT) = sum(RF.count(m));
        end
    end

    figure("Name", "Type by Ring Stacked Bar (Final Time)");
    bar(categorical(regList, regList, 'Ordinal', true), H, 'stacked', 'EdgeColor', 'none');
    xlabel("Ring");
    ylabel("Cell count");
    title(sprintf("Cell Type Amounts by Ring at t = %.1f min", tFinal));
    legend(cellstr(typeList), "Interpreter", "none", "Location", "eastoutside");
    grid on;
    exportPngNoToolbar(gcf, fullfile(figOutDir, "11_type_region_stackedbar_final_time.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Type-Ring Stacked Bar', figureCount, totalFigures), totalFigures);

    % Figure 12: Tumor vs stromal burden by region over time
    isTumor = arrayfun(@(x) isCancerLike(x), RC.type_label);
    isStroma = arrayfun(@(x) isStromaLike(x), RC.type_label);
    regions2 = ringLabels(activeRingIdx);

    figure("Name", "Tumor vs Stroma by Ring");
    nReg = numel(regions2);
    nCols = min(2, nReg);
    nRows = ceil(nReg / nCols);
    tiledlayout(nRows, nCols, "Padding", "compact", "TileSpacing", "compact");
    for iR = 1:nReg
        nexttile;
        rr = regions2(iR);
        mR = (RC.ring_idx == activeRingIdx(iR));

        tvals = unique(RC.time_min(mR), "sorted");
        tumorSeries = zeros(numel(tvals),1);
        stromaSeries = zeros(numel(tvals),1);

        for it = 1:numel(tvals)
            mt = mR & (RC.time_min == tvals(it));
            tumorSeries(it) = sum(RC.count(mt & isTumor));
            stromaSeries(it) = sum(RC.count(mt & isStroma));
        end

        plot(tvals, tumorSeries, "-", "LineWidth", 1.8); hold on;
        plot(tvals, stromaSeries, "-", "LineWidth", 1.8); hold off;
        title(sprintf("Ring: %s", rr), "Interpreter", "none");
        xlabel("Time (min)"); ylabel("Cell count");
        grid on;
        if iR == 1
            legend({"Tumor-like", "Stroma-like"}, "Location", "best", "Interpreter", "none");
        end
    end
    exportPngNoToolbar(gcf, fullfile(figOutDir, "12_tumor_vs_stroma_by_region.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Tumor vs Stroma by Ring', figureCount, totalFigures), totalFigures);

    % Figure 13: Cancer fraction by ring over time
    figure("Name", "Cancer Fraction by Ring Over Time");
    hold on;
    for iR = 1:numel(regions2)
        rr = regions2(iR);
        mR = (RC.ring_idx == activeRingIdx(iR));
        tvals = unique(RC.time_min(mR), "sorted");
        frac = zeros(numel(tvals),1);
        for it = 1:numel(tvals)
            mt = mR & (RC.time_min == tvals(it));
            tumorN = sum(RC.count(mt & isTumor));
            stromaN = sum(RC.count(mt & isStroma));
            denom = tumorN + stromaN;
            if denom > 0
                frac(it) = tumorN / denom;
            else
                frac(it) = 0;
            end
        end
        plot(tvals, frac, "LineWidth", 1.8, "DisplayName", rr);
    end
    hold off;
    ylim([0 1]);
    xlabel("Time (min)");
    ylabel("Cancer fraction = Tumor/(Tumor+Stroma)");
    title("Cancer Fraction by Ring Over Time");
    legend("Interpreter", "none", "Location", "bestoutside");
    grid on;
    exportPngNoToolbar(gcf, fullfile(figOutDir, "13_cancer_fraction_by_region_over_time.png"));
    figureCount = figureCount + 1;
    hwb = safeWaitbarUpdate(figureCount/totalFigures, hwb, sprintf('Generating figures: %d/%d - Cancer Fraction by Ring', figureCount, totalFigures), totalFigures);
end

% Close progress bar
if ~isempty(hwb) && isgraphics(hwb)
    close(hwb);
end

% Restore user/default figure visibility preference.
set(groot, 'defaultFigureVisible', defaultFigureVisibleOriginal);

disp('Analysis complete! All figures generated.');

% Show metadata
disp(" ");
disp("========== METADATA ==========");
disp(metadata);
