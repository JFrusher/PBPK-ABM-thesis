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
%   - Uses fixed standardized figure size for consistent visuals
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
        outputPath = sprintf('%s_barcode.png', baseName);
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

    bolusTimes = sort(startMin(bolusMask));

    continuousRate = zeros(size(startMin));
    if ismember('infusion_rate', dosing.Properties.VariableNames)
        continuousRate(continuousMask) = double(dosing.infusion_rate(continuousMask));
    elseif ismember('mean_rate', dosing.Properties.VariableNames)
        continuousRate(continuousMask) = double(dosing.mean_rate(continuousMask));
    elseif ismember('dose_amount', dosing.Properties.VariableNames)
        validDur = durationMin > 0;
        tempRate = zeros(size(startMin));
        tempRate(validDur) = double(dosing.dose_amount(validDur)) ./ durationMin(validDur);
        continuousRate(continuousMask) = tempRate(continuousMask);
    else
        continuousRate(continuousMask) = 1;
    end
    continuousRate(~isfinite(continuousRate) | continuousRate < 0) = 0;

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

    contSegments = [startMin(continuousMask), endMin(continuousMask), continuousRate(continuousMask)];
    if ~isempty(contSegments)
        contSegments = sortrows(contSegments, 1);
    end

    otherSegments = [startMin(otherWindowMask), endMin(otherWindowMask)];
    if ~isempty(otherSegments)
        otherSegments = sortrows(otherSegments, 1);
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

    % --- Standardized figure layout ---
    fig = figure('Color', 'w', 'Position', [100 100 1400 360], ...
                 'Name', 'Dosing Barcode', 'NumberTitle', 'off', ...
                 'Visible', 'on');
    ax = axes(fig, 'Position', [0.07 0.20 0.78 0.70]);
    hold(ax, 'on');

    % Baseline timeline
    line(ax, [tMin tMax] / 60, [0 0], 'Color', [0.15 0.15 0.15], 'LineWidth', 1.0);

    % Draw low-intensity windows for non-continuous, non-bolus entries
    for i = 1:size(otherSegments,1)
        x0 = otherSegments(i,1) / 60;
        x1 = otherSegments(i,2) / 60;
        w = max(x1 - x0, 1e-6);
        rectangle(ax, 'Position', [x0, 0, w, 0.22], ...
                  'FaceColor', [0.95 0.65 0.20], ...
                  'EdgeColor', [0.75 0.48 0.12], ...
                  'LineWidth', 0.8, ...
                  'Curvature', [0.02 0.02]);
    end

    % Draw continuous windows with relative height (key visual signal)
    if ~isempty(contSegments)
        maxRate = max(contSegments(:,3));
        if maxRate <= 0
            maxRate = 1;
        end

        for i = 1:size(contSegments,1)
            x0 = contSegments(i,1) / 60;
            x1 = contSegments(i,2) / 60;
            w = max(x1 - x0, 1e-6);
            h = contSegments(i,3) / maxRate;
            h = max(0.06, min(1.0, h));

            rectangle(ax, 'Position', [x0, 0, w, h], ...
                      'FaceColor', [0.20 0.60 0.95], ...
                      'EdgeColor', [0.10 0.35 0.65], ...
                      'LineWidth', 0.6, ...
                      'Curvature', [0.02 0.02]);
        end
    else
        maxRate = 0;
    end

    if ~isempty(bolusTimes)
        xHours = bolusTimes / 60;
        for i = 1:numel(xHours)
            line(ax, [xHours(i) xHours(i)], [0 1], ...
                'Color', [0.90 0.25 0.25], 'LineWidth', 2.4);
        end
        scatter(ax, xHours, repmat(1, size(xHours)), 46, ...
                'MarkerFaceColor', [0.90 0.25 0.25], ...
                'MarkerEdgeColor', 'none');
    end

    xlim(ax, [tMin tMax] / 60);
    ylim(ax, [0 1.15]);
    yticks(ax, [0 0.25 0.50 0.75 1.0]);
    yticklabels(ax, {'0', '0.25', '0.50', '0.75', '1.00'});

    xlabel(ax, 'Time (hours)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax, 'Relative dosing intensity', 'FontSize', 12, 'FontWeight', 'bold');

    [~, titleName, ~] = fileparts(csvFilename);
    title(ax, sprintf('Dosing Barcode: %s', titleName), ...
            'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

    grid(ax, 'on');
    ax.GridAlpha = 0.20;
    ax.Layer = 'top';

    exportgraphics(fig, outputPath, 'Resolution', 200);
    fprintf('Saved barcode plot: %s\n', outputPath);
    fprintf('Total dose: %.3f mg\n', doseSummary.totalDose_mg);
    fprintf('Total dose/day: %.3f mg/day (span %.3f days)\n', doseSummary.totalDosePerDay_mg_day, doseSummary.span_days);
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

