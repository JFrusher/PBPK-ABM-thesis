function [outputPath, doseSummary] = plotDosingBarcode(csvFilename, outputPath, exportPhysiCellShifted)
%PLOTDOSINGBARCODE Create standardized barcode plot for dosing windows.
%   outputPath = plotDosingBarcode(csvFilename)
%   outputPath = plotDosingBarcode(csvFilename, outputPath)
%   outputPath = plotDosingBarcode(csvFilename, outputPath, exportPhysiCellShifted)
%   [outputPath, doseSummary] = plotDosingBarcode(...)
%
% Inputs:
%   csvFilename : Path to dosing CSV with at least:
%                 start_time_min,end_time_min,dosing_type
%   outputPath  : Optional PNG path. Default: <csvName>_barcode.png
%   exportPhysiCellShifted : Optional logical flag (default false). If true,
%                 writes a separate CSV with start/end times shifted by -1000
%                 minutes for Physicell/cellular run-in workflows.
%
% Behavior:
%   - Uses a one-line barcode timeline
%   - Encodes continuous-dose intensity as relative bar height
%   - Keeps bolus doses as event ticks
%   - Uses publication-oriented panel sizing for A4 4x2 figure grids
%   - Computes total dose and dose/day from CSV timing + dosing fields

    if nargin < 1 || (~ischar(csvFilename) && ~isstring(csvFilename))
        error('csvFilename must be a character vector or string.');
    end

    csvFilename = char(csvFilename);
    if ~isfile(csvFilename)
        error('CSV file not found: %s', csvFilename);
    end

    if nargin < 2 || isempty(outputPath)
        [~, baseName, ~] = fileparts(csvFilename);
        scheduleLabel = makeScheduleLabel(baseName);
        outputPath = sprintf('%s_barcode.png', scheduleLabel);
    end
    [outFolder, outName, outExt] = fileparts(outputPath);
    if isempty(outExt)
        outputPath = fullfile(outFolder, [outName '.png']);
    end

    if nargin < 3 || isempty(exportPhysiCellShifted)
        exportPhysiCellShifted = false;
    end
    exportPhysiCellShifted = logical(exportPhysiCellShifted);

    dosing = readtable(csvFilename);

    if exportPhysiCellShifted
        shiftedCsvPath = exportPhysiCellVersion(dosing, csvFilename, 1000);
        fprintf('Saved Physicell-shifted CSV (-1000 min): %s\n', shiftedCsvPath);
    end

    requiredColumns = {'start_time_min', 'end_time_min', 'dosing_type'};
    for i = 1:numel(requiredColumns)
        if ~ismember(requiredColumns{i}, dosing.Properties.VariableNames)
            error('Missing required column: %s', requiredColumns{i});
        end
    end

    startMin = double(dosing.start_time_min(:));
    endMin = double(dosing.end_time_min(:));
    typeRaw = dosing.dosing_type;

    % Normalize dosing type to lowercase strings
    if iscell(typeRaw)
        type = lower(strtrim(string(typeRaw)));
    else
        type = lower(strtrim(string(typeRaw(:))));
    end

    % Guard against inverted ranges
    swapMask = endMin < startMin;
    if any(swapMask)
        tmp = startMin(swapMask);
        startMin(swapMask) = endMin(swapMask);
        endMin(swapMask) = tmp;
    end

    durationMin = endMin - startMin;

    % Interpret instantaneous or near-instantaneous rows as bolus events
    bolusMask = type == "bolus" | durationMin <= 1e-9;

    % Continuous-like windows (chronotherapy half-hour rows become a single
    % contiguous visual block because adjacent rectangles touch)
    continuousMask = (type == "continuous" | type == "constant") & ~bolusMask;

    % Other windowed dosing types (e.g., sinusoidal/custom)
    otherWindowMask = ~bolusMask & ~continuousMask;

    bolusTimes = startMin(bolusMask);

    bolusValues = ones(size(startMin));
    if ismember('dose_amount', dosing.Properties.VariableNames)
        bolusValues = double(dosing.dose_amount(:));
    end
    bolusValues(~isfinite(bolusValues) | bolusValues < 0) = 0;

    bolusDoseValues = bolusValues(bolusMask);
    if ~isempty(bolusTimes)
        [bolusTimes, bolusOrder] = sort(bolusTimes);
        bolusDoseValues = bolusDoseValues(bolusOrder);
    end

    windowValue = zeros(size(startMin));
    if ismember('infusion_rate', dosing.Properties.VariableNames)
        windowValue = double(dosing.infusion_rate(:));
    elseif ismember('mean_rate', dosing.Properties.VariableNames)
        windowValue = double(dosing.mean_rate(:));
    elseif ismember('dose_amount', dosing.Properties.VariableNames)
        validDur = durationMin > 0;
        windowValue(validDur) = double(dosing.dose_amount(validDur)) ./ durationMin(validDur);
    else
        windowValue(:) = 1;
    end
    windowValue(~isfinite(windowValue) | windowValue < 0) = 0;

    continuousRate = windowValue(continuousMask);
    otherWindowValue = windowValue(otherWindowMask);

    % Dose summary: rate-based rows use mg/min * minutes, bolus uses dose_amount (mg)
    bolusDoseMg = 0;
    if ismember('dose_amount', dosing.Properties.VariableNames)
        bolusDoseVals = double(dosing.dose_amount(bolusMask));
        bolusDoseVals(~isfinite(bolusDoseVals) | bolusDoseVals < 0) = 0;
        bolusDoseMg = sum(bolusDoseVals);
    end

    rateDoseMg = 0;
    if ismember('infusion_rate', dosing.Properties.VariableNames)
        rowRate = double(dosing.infusion_rate(:));
        rowRate(~isfinite(rowRate) | rowRate < 0) = 0;
        rateDoseMg = sum(rowRate(~bolusMask) .* max(durationMin(~bolusMask), 0));
    end

    totalDoseMg = bolusDoseMg + rateDoseMg;

    contStart = startMin(continuousMask);
    contEnd = endMin(continuousMask);
    contVal = continuousRate(:);
    nCont = min([numel(contStart), numel(contEnd), numel(contVal)]);
    if nCont > 0
        contSegments = [contStart(1:nCont), contEnd(1:nCont), contVal(1:nCont)];
        contSegments = sortrows(contSegments, 1);
    else
        contSegments = zeros(0, 3);
    end

    otherStart = startMin(otherWindowMask);
    otherEnd = endMin(otherWindowMask);
    otherVal = otherWindowValue(:);
    nOther = min([numel(otherStart), numel(otherEnd), numel(otherVal)]);
    if nOther > 0
        otherSegments = [otherStart(1:nOther), otherEnd(1:nOther), otherVal(1:nOther)];
        otherSegments = sortrows(otherSegments, 1);
    else
        otherSegments = zeros(0, 3);
    end

    hasWindowDose = ~isempty(contSegments) || ~isempty(otherSegments);
    hasBolusDose = ~isempty(bolusTimes);

    maxWindowDose = max([contVal; otherVal; 0]);
    if maxWindowDose <= 0
        maxWindowDose = 1;
    end

    maxBolusDose = max([bolusDoseValues; 0]);
    if maxBolusDose <= 0
        maxBolusDose = 1;
    end

    allTimes = [startMin; endMin];
    tMin = min(allTimes);
    tMax = max(allTimes);

    if isempty(allTimes) || ~isfinite(tMin) || ~isfinite(tMax)
        error('No valid timing data found in CSV.');
    end

    if tMax <= tMin
        tMax = tMin + 60;
    end

    spanMinutes = max(tMax - tMin, 1);
    spanDays = spanMinutes / 1440;
    dosePerDayMg = totalDoseMg / spanDays;

    doseSummary = struct();
    doseSummary.totalDose_mg = totalDoseMg;
    doseSummary.totalDosePerDay_mg_day = dosePerDayMg;
    doseSummary.bolusDose_mg = bolusDoseMg;
    doseSummary.rateDose_mg = rateDoseMg;
    doseSummary.span_days = spanDays;

    % --- Publication-oriented panel sizing (one slot in an A4 4x2 layout) ---
    panelSizeCm = getA4GridPanelSizeCm(4, 2);
    fig = figure('Color', 'w', ...
                 'Name', 'Dosing Barcode', 'NumberTitle', 'off', ...
                 'Visible', 'on', ...
                 'Units', 'centimeters', ...
                 'Position', [2 2 panelSizeCm(1) panelSizeCm(2)], ...
                 'PaperUnits', 'centimeters', ...
                 'PaperPositionMode', 'auto');
    ax = axes(fig, 'Position', [0.10 0.18 0.86 0.72]);
    hold(ax, 'on');

    baseFont = 8.5;
    showLegend = strcmpi(getenv('PBPK_SHOW_BARCODE_LEGEND'), '1');
    set(ax, 'FontName', 'Helvetica', 'FontSize', baseFont, 'LineWidth', 0.9, ...
            'TickDir', 'out', 'Box', 'off', 'Layer', 'top');

    % Colorblind-safe, print-friendly palette.
    cContinuousFace = [0.20 0.48 0.72];
    cContinuousEdge = [0.11 0.30 0.47];
    cWindowFace = [0.93 0.68 0.24];
    cWindowEdge = [0.72 0.49 0.12];
    cBolus = [0.72 0.20 0.20];
    cBaseline = [0.18 0.18 0.18];

    hasWindowMagnitude = hasWindowDose && any([contVal; otherVal] > 0);
    hasBolusMagnitude = hasBolusDose && any(bolusDoseValues > 0);
    bothDoseTypes = hasWindowMagnitude && hasBolusMagnitude;

    % Baseline timeline
    % If there is a dose at tMin (e.g., 0), pad the left axis a bit so it's visible
    padLeft = 0;
    if any(abs([startMin; endMin] - tMin) < 1e-6)
        padLeft = 0.5; % hours
    end
    line(ax, [tMin-padLeft tMax] / 60, [0 0], 'Color', cBaseline, 'LineWidth', 1.0);

    % Draw low-intensity windows for non-continuous, non-bolus entries
    hOther = gobjects(0);
    for i = 1:size(otherSegments,1)
        x0 = otherSegments(i,1) / 60;
        x1 = otherSegments(i,2) / 60;
        w = max(x1 - x0, 1e-6);
        hRect = rectangle(ax, 'Position', [x0, 0, w, otherSegments(i,3)], ...
                  'FaceColor', cWindowFace, ...
                  'EdgeColor', cWindowEdge, ...
                  'LineWidth', 0.8, ...
                  'Curvature', [0.02 0.02]);
        if i == 1
            hOther = hRect;
        end
    end

    % Draw continuous windows with relative height (key visual signal)
    hCont = gobjects(0);
    if ~isempty(contSegments)
        for i = 1:size(contSegments,1)
            x0 = contSegments(i,1) / 60;
            x1 = contSegments(i,2) / 60;
            w = max(x1 - x0, 1e-6);
            hRect = rectangle(ax, 'Position', [x0, 0, w, contSegments(i,3)], ...
                      'FaceColor', cContinuousFace, ...
                      'EdgeColor', cContinuousEdge, ...
                      'LineWidth', 0.6, ...
                      'Curvature', [0.02 0.02]);
            if i == 1
                hCont = hRect;
            end
        end
    end

    hBolus = gobjects(0);
    if bothDoseTypes
        yyaxis(ax, 'right');
    end

    if ~isempty(bolusTimes)
        xHours = bolusTimes / 60;
        yBolus = bolusDoseValues;
        for i = 1:numel(xHours)
            hLine = line(ax, [xHours(i) xHours(i)], [0 yBolus(i)], ...
                'Color', cBolus, 'LineWidth', 1.8, 'LineStyle', '-');
            if i == 1
                hBolus = hLine;
            end
        end
        scatter(ax, xHours, yBolus, 46, ...
                'MarkerFaceColor', cBolus, ...
                'MarkerEdgeColor', 'none');
    end

    xlim(ax, [tMin tMax] / 60);
    if bothDoseTypes
        yyaxis(ax, 'left');
        ylim(ax, [0 1.08 * maxWindowDose]);
        setDoseAxisTicks(ax, maxWindowDose);
        ylabel(ax, 'Infusion/window value', 'FontSize', baseFont + 1, 'FontWeight', 'bold');

        yyaxis(ax, 'right');
        ylim(ax, [0 1.08 * maxBolusDose]);
        setDoseAxisTicks(ax, maxBolusDose);
        ylabel(ax, 'Bolus dose', 'FontSize', baseFont + 1, 'FontWeight', 'bold');
    elseif hasWindowDose
        ylim(ax, [0 1.08 * maxWindowDose]);
        setDoseAxisTicks(ax, maxWindowDose);
        ylabel(ax, 'Infusion/window value', 'FontSize', baseFont + 1, 'FontWeight', 'bold');
    else
        ylim(ax, [0 1.08 * maxBolusDose]);
        setDoseAxisTicks(ax, maxBolusDose);
        ylabel(ax, 'Bolus dose', 'FontSize', baseFont + 1, 'FontWeight', 'bold');
    end

    if ~bothDoseTypes
        try
            if numel(ax.YAxis) > 1
                ax.YAxis(2).Visible = 'off';
            end
        catch
        end
    end

    applyPublicationTimeTicks(ax, tMin, tMax);

    xlabel(ax, 'Time', 'FontSize', baseFont + 1, 'FontWeight', 'bold');

    [~, titleName, ~] = fileparts(csvFilename);
    scheduleLabel = makeScheduleLabel(titleName);
    title(ax, sprintf('Dosing Schedule: %s', scheduleLabel), ...
            'FontSize', baseFont + 2, 'FontWeight', 'bold', 'Interpreter', 'none');

    legendHandles = gobjects(0);
    legendLabels = {};
    if ~isempty(hCont)
        legendHandles(end+1) = hCont; %#ok<AGROW>
        legendLabels{end+1} = 'Continuous/constant window'; %#ok<AGROW>
    end
    if ~isempty(hOther)
        legendHandles(end+1) = hOther; %#ok<AGROW>
        legendLabels{end+1} = 'Other dosing window'; %#ok<AGROW>
    end
    if ~isempty(hBolus)
        legendHandles(end+1) = hBolus; %#ok<AGROW>
        legendLabels{end+1} = 'Bolus event'; %#ok<AGROW>
    end
    if showLegend && ~isempty(legendHandles)
        lgd = legend(ax, legendHandles, legendLabels, 'Location', 'northoutside', ...
                     'Orientation', 'horizontal', 'Box', 'off');
        lgd.FontSize = baseFont;
    end

    grid(ax, 'on');
    ax.GridAlpha = 0.18;
    ax.MinorGridAlpha = 0.10;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'on';
    ax.GridColor = [0.2 0.2 0.2];

    if bothDoseTypes
        yyaxis(ax, 'right');
        ax.YColor = cBolus;
        yyaxis(ax, 'left');
    end

    exportgraphics(fig, outputPath, 'Resolution', 600);
    fprintf('Saved barcode plot: %s\n', outputPath);
    fprintf('Total dose: %.3f mg\n', doseSummary.totalDose_mg);
    fprintf('Total dose/day: %.3f mg/day (span %.3f days)\n', doseSummary.totalDosePerDay_mg_day, doseSummary.span_days);
end

function applyPublicationTimeTicks(ax, tMin, tMax)
    spanHours = max((tMax - tMin) / 60, 1e-6);
    if spanHours <= 72
        tickStepHours = 12;
        if spanHours <= 24
            tickStepHours = 6;
        end
        t0 = floor((tMin / 60) / tickStepHours) * tickStepHours;
        t1 = ceil((tMax / 60) / tickStepHours) * tickStepHours;
        ticks = t0:tickStepHours:t1;
        xticks(ax, ticks);
        xtickformat(ax, '%.0f h');
    else
        tickStepDays = 1;
        if spanHours > 14 * 24
            tickStepDays = 2;
        end
        t0d = floor((tMin / 1440) / tickStepDays) * tickStepDays;
        t1d = ceil((tMax / 1440) / tickStepDays) * tickStepDays;
        ticksDays = t0d:tickStepDays:t1d;
        xticks(ax, ticksDays * 24);
        labels = compose('Day %d', round(ticksDays));
        xticklabels(ax, labels);
    end
end

function setDoseAxisTicks(ax, maxDoseValue)
    nTicks = 5;
    ticks = linspace(0, maxDoseValue, nTicks);
    yticks(ax, ticks);

    if maxDoseValue >= 100
        ytickformat(ax, '%.0f');
    elseif maxDoseValue >= 1
        ytickformat(ax, '%.2g');
    else
        ytickformat(ax, '%.2f');
    end
end

function panelSizeCm = getA4GridPanelSizeCm(nRows, nCols)
    % Designed for portrait A4 with practical margins/gaps in dissertations.
    pageW = 21.0;
    pageH = 29.7;
    marginL = 1.2;
    marginR = 1.2;
    marginT = 1.2;
    marginB = 1.2;
    gapW = 0.6;
    gapH = 0.55;

    usableW = pageW - marginL - marginR - gapW * (nCols - 1);
    usableH = pageH - marginT - marginB - gapH * (nRows - 1);
    tileW = usableW / nCols;
    tileH = usableH / nRows;
    panelSizeCm = [tileW, tileH];
end

function scheduleLabel = makeScheduleLabel(baseName)
    baseName = char(string(baseName));
    tok = regexp(baseName, '^\s*(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        scheduleLabel = sprintf('S%d', str2double(tok{1}));
    else
        scheduleLabel = strrep(baseName, '_', ' ');
    end
end

function shiftedCsvPath = exportPhysiCellVersion(dosing, sourceCsvPath, shiftMinutes)
    shifted = dosing;

    if ismember('start_time_min', shifted.Properties.VariableNames)
        shifted.start_time_min = shifted.start_time_min - shiftMinutes;
    else
        error('Missing required column: start_time_min');
    end

    if ismember('end_time_min', shifted.Properties.VariableNames)
        shifted.end_time_min = shifted.end_time_min - shiftMinutes;
    else
        error('Missing required column: end_time_min');
    end

    [folderPath, baseName, ~] = fileparts(sourceCsvPath);
    shiftedCsvPath = fullfile(folderPath, sprintf('%s_physicell_shifted.csv', baseName));
    writetable(shifted, shiftedCsvPath);
end

