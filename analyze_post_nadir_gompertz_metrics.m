%% analyze_post_nadir_gompertz_metrics.m
% Batch analysis of ABM tumor-growth CSVs:
% 1) Locate nadir after dosing starts
% 2) Fit post-nadir regrowth with Gompertz model
% 3) Compute clinical-style metrics
% 4) Create a 2-panel figure per file + detailed diagnostics TXT

%clear; clc;

%% ------------------------------- USER SETTINGS -------------------------------
% Option A: Explicit list
csvFiles = {
'WP1M_CSV/curve_fitter_highlights.csv'
    % 'Sim_results_8_csv/curve_fitter_highlights.csv'
};

% Option B: Pattern search when csvFiles is empty
if isempty(csvFiles)
    pattern = fullfile('**', 'curve_fitter_highlights.csv');
    d = dir(pattern);
    csvFiles = cellfun(@(f,p) fullfile(f,p), {d.folder}, {d.name}, 'UniformOutput', false);
end

% Dosing start time in same units as the time column.
doseStartTime = 0;

% Export outputs
outputDir = 'post_nadir_gompertz_analysis';
saveFigures = true;
figureFormat = 'png';

% Interactive dip picker (CLI)
interactiveDipPicker = true;
maxDipCandidates = 8;
dipProminenceFraction = 0.05;     % fraction of post-dose range
minDipSeparationPoints = 3;       % minimum index spacing between listed dips

if ~isfolder(outputDir)
    mkdir(outputDir);
end

if isempty(csvFiles)
    error('No CSV files found. Populate csvFiles or verify the search pattern.');
end

%% ----------------------------- ANALYSIS LOOP --------------------------------
allRows = [];

for i = 1:numel(csvFiles)
    csvPath = char(csvFiles{i});
    if ~isfile(csvPath)
        warning('Skipping missing file: %s', csvPath);
        continue;
    end

    T = readtable(csvPath);
    [timeVals, totalVals, liveVals, deadVals, prolifVals, seriesMeta] = extractCoreSeries(T);

    if numel(timeVals) < 6
        warning('Skipping %s (not enough time points).', csvPath);
        continue;
    end

    validMask = isfinite(timeVals) & isfinite(totalVals) & totalVals > 0;
    timeVals = timeVals(validMask);
    totalVals = totalVals(validMask);
    liveVals = liveVals(validMask);
    deadVals = deadVals(validMask);
    prolifVals = prolifVals(validMask);

    [timeVals, sortIdx] = sort(timeVals(:));
    totalVals = totalVals(sortIdx);
    liveVals = liveVals(sortIdx);
    deadVals = deadVals(sortIdx);
    prolifVals = prolifVals(sortIdx);

    idx0 = find(timeVals >= doseStartTime, 1, 'first');
    if isempty(idx0)
        warning('Skipping %s (no time point at/after doseStartTime).', csvPath);
        continue;
    end

    dosingMask = timeVals >= doseStartTime;
    idxDosing = find(dosingMask);

    % Analysis nadir (for total-cell metrics): absolute global min after dosing start.
    dosingVals = totalVals(dosingMask);
    nadirVal = min(dosingVals);
    % If the same absolute minimum occurs multiple times, use the last one.
    % This avoids anchoring to an earlier trough when a later point is equally low.
    tieTol = max(1e-12, 1e-10 * max(abs(nadirVal), 1));
    relMinAll = find(abs(dosingVals - nadirVal) <= tieTol);
    nRel = relMinAll(end);
    idxNadir = idxDosing(nRel);
    tNadir = timeVals(idxNadir);

    % Plot/evaluator nadir: use live cells when available, else total cells.
    idxNadirPlot = idxNadir;
    if seriesMeta.hasLiveSeries && any(isfinite(liveVals(dosingMask)) & liveVals(dosingMask) > 0)
        dosingLive = liveVals(dosingMask);
        liveNadirVal = min(dosingLive);
        tieTolLive = max(1e-12, 1e-10 * max(abs(liveNadirVal), 1));
        relLiveMinAll = find(abs(dosingLive - liveNadirVal) <= tieTolLive);
        idxNadirPlot = idxDosing(relLiveMinAll(end));
    end

    % Optional interactive override: choose among major dips for inspection.
    if seriesMeta.hasLiveSeries && any(isfinite(liveVals(dosingMask)) & liveVals(dosingMask) > 0)
        dipSeries = liveVals;
        dipSeriesLabel = 'Live Cells';
    else
        dipSeries = totalVals;
        dipSeriesLabel = 'Total Cells';
    end

    dipCandidates = findMajorDipCandidates(timeVals, dipSeries, idxDosing, maxDipCandidates, dipProminenceFraction, minDipSeparationPoints);
    if interactiveDipPicker
        idxNadirPlot = chooseDipViaCLI(dipCandidates, idxNadirPlot, csvPath, dipSeriesLabel);
    end

    postMask = (1:numel(timeVals))' >= idxNadir;
    tPost = timeVals(postMask) - tNadir;
    vPost = totalVals(postMask);

    [Kfit, alphaFit, vFit, fitOk] = fitGompertzPostNadir(tPost, vPost, nadirVal);

    v0 = totalVals(idx0);
    peakVal = max(totalVals);

    % Primary evaluator: time to recover baseline live-cell count
    if seriesMeta.hasLiveSeries && isfinite(liveVals(idx0)) && liveVals(idx0) > 0
        live0 = liveVals(idx0);
        [timeToLiveRegrowth, liveRegrowthReason] = crossingTimeFromNadir(timeVals, liveVals, idxNadirPlot, live0);
    else
        live0 = NaN;
        timeToLiveRegrowth = NaN;
        if ~seriesMeta.hasLiveSeries
            liveRegrowthReason = 'Live-cell column not found';
        else
            liveRegrowthReason = 'Baseline live-cell value at t0 is non-finite or non-positive';
        end
    end

    % Doubling time from nadir: first crossing of 2*nadir
    doublingTarget = 2 * nadirVal;
    [doublingTime, doublingReason] = crossingTimeFromNadir(timeVals, totalVals, idxNadir, doublingTarget);

    % Fractional reduction
    fractionalReduction = (v0 - nadirVal) / v0;

    % Log-kill spread: log10(peak / nadir)
    logKillSpread = log10(peakVal / nadirVal);

    % Necrotic fraction shift: (dead/total at nadir) - (dead/total at t0)
    if seriesMeta.hasDeadSeries && isfinite(deadVals(idx0)) && isfinite(deadVals(idxNadir))
        necFrac0 = deadVals(idx0) / v0;
        necFracNadir = deadVals(idxNadir) / nadirVal;
        necroticFractionShift = necFracNadir - necFrac0;
    else
        necFrac0 = NaN;
        necFracNadir = NaN;
        necroticFractionShift = NaN;
    end

    % Growth fraction (GF) over recovery phase, then summarize as mean GF
    gfRecovery = prolifVals(postMask) ./ vPost;
    gfRecovery(~isfinite(gfRecovery)) = NaN;
    growthFractionMeanRecovery = mean(gfRecovery, 'omitnan');
    if isempty(gfRecovery) || all(~isfinite(gfRecovery))
        growthFractionMeanRecovery = NaN;
    end

    % Time to progression (TTP): from nadir to 120% of original t0 volume
    ttpTarget = 1.2 * v0;
    [ttp, ttpReason] = crossingTimeFromNadir(timeVals, totalVals, idxNadir, ttpTarget);

    % Regression rate: steepest negative slope during dosing phase up to nadir
    idxReg = idx0:idxNadir;
    if numel(idxReg) >= 2
        dv = diff(totalVals(idxReg));
        dt = diff(timeVals(idxReg));
        slopes = dv ./ dt;
        slopes(~isfinite(slopes)) = NaN;
        regressionRate = min(slopes, [], 'omitnan');
        if isempty(slopes) || all(~isfinite(slopes))
            regressionRate = NaN;
        end
    else
        regressionRate = NaN;
    end

    % 8-metric vector (includes Gompertz alpha as regrowth metric)
    metricNames = {
        'Gompertz alpha', ...
        'Doubling time', ...
        'Fractional reduction', ...
        'Log-kill spread', ...
        'Necrotic frac shift', ...
        'Mean GF recovery', ...
        'TTP to 120% V0', ...
        'Regression rate'
    };
    metricVals = [
        alphaFit, ...
        doublingTime, ...
        fractionalReduction, ...
        logKillSpread, ...
        necroticFractionShift, ...
        growthFractionMeanRecovery, ...
        ttp, ...
        regressionRate
    ];

    metricReasons = buildMetricDiagnostics(metricVals, fitOk, seriesMeta, doublingReason, ttpReason, idx0, idxNadir, gfRecovery, v0, nadirVal);

    missingMask = ~isfinite(metricVals);
    if any(missingMask)
        fprintf('Diagnostics for %s:\n', csvPath);
        for k = 1:numel(metricNames)
            if missingMask(k)
                fprintf('  - %s: %s\n', metricNames{k}, metricReasons{k});
            end
        end
    end
    if ~isfinite(timeToLiveRegrowth)
        fprintf('  - Primary evaluator (time to live-cell regrowth): %s\n', liveRegrowthReason);
    end

    % Build summary row
    [~, fileStem, ~] = fileparts(csvPath);
    rowVarNames = {
        'csv_path','sample','gompertz_fit_ok','V0','V_nadir','V_peak','K_fit','alpha_fit', ...
        'live0','time_to_live_regrowth','time_to_live_regrowth_reason', ...
        'doubling_time_from_nadir','fractional_reduction','log_kill_spread', ...
        'necrotic_fraction_shift','growth_fraction_mean_recovery','TTP_120pct_V0', ...
        'regression_rate','metric_diagnostics'
    };
    row = table( ...
        string(csvPath), string(fileStem), fitOk, v0, nadirVal, peakVal, Kfit, alphaFit, ...
        live0, timeToLiveRegrowth, string(liveRegrowthReason), ...
        doublingTime, fractionalReduction, logKillSpread, necroticFractionShift, ...
        growthFractionMeanRecovery, ttp, regressionRate, string(strjoin(metricReasons, ' | ')), ...
        'VariableNames', rowVarNames);

    if isempty(allRows)
        allRows = row;
    else
        allRows = [allRows; row]; %#ok<AGROW>
    end

    % ------------------------- Visualization (2 panels) ----------------------
    hasLiveForPlot = seriesMeta.hasLiveSeries && any(isfinite(liveVals) & liveVals > 0);
    if hasLiveForPlot
        plotVals = liveVals;
        plotSourceLabel = 'Live Cells';
        plotNadirIdx = idxNadirPlot;
    else
        plotVals = totalVals;
        plotSourceLabel = 'Total Cells (fallback: live missing)';
        plotNadirIdx = idxNadir;
    end

    plotNadirVal = plotVals(plotNadirIdx);
    plotNadirTime = timeVals(plotNadirIdx);
    postMaskPlot = (1:numel(timeVals))' >= plotNadirIdx;
    tPostPlot = timeVals(postMaskPlot) - plotNadirTime;
    plotPostVals = plotVals(postMaskPlot);
    [plotKfit, plotAlphaFit, plotFitVals, plotFitOk] = fitGompertzPostNadir(tPostPlot, plotPostVals, max(plotNadirVal, eps));

    fig = figure('Color','w','Position',[80 80 1250 460]);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Panel 1: Full live-cell trajectory with nadir marker
    ax1 = nexttile(tl, 1); hold(ax1, 'on');
    plot(ax1, timeVals, plotVals, 'k-', 'LineWidth', 1.5);
    xline(ax1, doseStartTime, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.0, 'Label', 'Dose start');
    scatter(ax1, plotNadirTime, plotNadirVal, 60, 'r', 'filled');
    title(ax1, sprintf('%s Curve with Nadir', plotSourceLabel));
    xlabel(ax1, 'Time');
    ylabel(ax1, plotSourceLabel);
    grid(ax1, 'on');

    % Panel 2: Post-nadir regrowth with Gompertz overlay (same series as panel 1)
    ax2 = nexttile(tl, 2); hold(ax2, 'on');
    plot(ax2, tPostPlot, plotPostVals, 'bo', 'MarkerSize', 4, 'DisplayName', sprintf('Observed post-nadir (%s)', plotSourceLabel));
    if plotFitOk
        plot(ax2, tPostPlot, plotFitVals, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Gompertz fit (%s)', plotSourceLabel));
    end
    title(ax2, sprintf('Post-Nadir Gompertz Fit (%s)', plotSourceLabel));
    xlabel(ax2, 'Time since Nadir');
    ylabel(ax2, plotSourceLabel);
    legend(ax2, 'Location', 'best');
    grid(ax2, 'on');

    sgtitle(tl, sprintf('Regrowth Analysis: %s', fileStem), 'Interpreter', 'none');

    figOut = 'not saved';
    if saveFigures
        figOut = fullfile(outputDir, sprintf('%s_regrowth_metrics.%s', fileStem, figureFormat));
        exportgraphics(fig, figOut, 'Resolution', 250, 'BackgroundColor', 'white');
    end

    diagTxtPath = fullfile(outputDir, sprintf('%s_regrowth_diagnostics.txt', fileStem));
    writeDiagnosticReport(diagTxtPath, csvPath, metricNames, metricVals, metricReasons, ...
        fitOk, Kfit, alphaFit, v0, nadirVal, peakVal, tNadir, doseStartTime, ...
        doublingTarget, ttpTarget, necFrac0, necFracNadir, live0, timeToLiveRegrowth, ...
        liveRegrowthReason, seriesMeta, timeVals, idx0, idxNadir, gfRecovery, figOut);
end

if isempty(allRows)
    warning('No analyzable CSVs were processed.');
else
    outTable = fullfile(outputDir, 'post_nadir_gompertz_metrics_summary.csv');
    writetable(allRows, outTable);
    fprintf('Saved summary table: %s\n', outTable);
end

%% ------------------------------- FUNCTIONS ----------------------------------
function [t, v, live, dead, prolif, meta] = extractCoreSeries(T)
    tCol = pickColumn(T, {'time_min','time','t_min','t'});
    vCol = pickColumn(T, {'total_cells','tumor_volume','volume','cell_count'});
    liveCol = pickColumn(T, {'live_cells'});
    deadCol = pickColumn(T, {'dead_cells','death_phase_cells','dead','necrotic_cells'});
    prolifCol = pickColumn(T, {'proliferative_phase_cells','proliferative_cells','cycling_cells'});

    if isempty(tCol) || isempty(vCol)
        error('Required columns not found. Need time and total-cell columns.');
    end

    t = toNumeric(T.(tCol));
    v = toNumeric(T.(vCol));

    meta = struct();
    meta.timeColumn = tCol;
    meta.totalColumn = vCol;
    meta.liveColumn = liveCol;
    meta.deadColumn = deadCol;
    meta.prolifColumn = prolifCol;
    meta.hasLiveSeries = ~isempty(liveCol);
    meta.hasDeadSeries = ~isempty(deadCol);
    meta.hasProlifSeries = ~isempty(prolifCol);

    if isempty(liveCol)
        live = nan(height(T),1);
    else
        live = toNumeric(T.(liveCol));
    end

    if isempty(deadCol)
        dead = nan(height(T),1);
    else
        dead = toNumeric(T.(deadCol));
    end

    if isempty(prolifCol)
        prolif = nan(height(T),1);
    else
        prolif = toNumeric(T.(prolifCol));
    end
end

function [Kfit, alphaFit, vFit, fitOk] = fitGompertzPostNadir(tPost, vPost, vNadir)
    fitOk = false;
    Kfit = NaN;
    alphaFit = NaN;
    vFit = nan(size(vPost));

    if numel(tPost) < 4 || any(~isfinite(tPost)) || any(~isfinite(vPost))
        return;
    end

    t = tPost(:);
    y = vPost(:);
    v0 = max(vNadir, eps);

    K0 = max(max(y)*1.2, v0*1.01);
    alpha0 = 1 / max(range(t), 1);

    p0 = [log(max(K0 - v0, eps)); log(max(alpha0, eps))];

    obj = @(p) sum((gompertzModel(t, v0 + exp(p(1)), exp(p(2)), v0) - y).^2, 'omitnan');
    opts = optimset('Display', 'off', 'MaxIter', 5000, 'MaxFunEvals', 10000);

    try
        p = fminsearch(obj, p0, opts);
        Kfit = v0 + exp(p(1));
        alphaFit = exp(p(2));
        vFit = gompertzModel(t, Kfit, alphaFit, v0);
        fitOk = isfinite(Kfit) && isfinite(alphaFit) && Kfit > v0 && alphaFit > 0;
    catch
        fitOk = false;
    end
end

function y = gompertzModel(t, K, alpha, V0)
    y = K .* exp(log(V0 ./ K) .* exp(-alpha .* t));
end

function [dtFromNadir, reason] = crossingTimeFromNadir(t, v, idxNadir, threshold)
    dtFromNadir = NaN;
    reason = 'Target threshold not crossed after nadir';
    tN = t(idxNadir);

    for k = idxNadir:(numel(t)-1)
        v1 = v(k);
        v2 = v(k+1);
        if ~isfinite(v1) || ~isfinite(v2)
            continue;
        end
        if (v1 < threshold && v2 >= threshold) || (v1 > threshold && v2 <= threshold)
            if v2 == v1
                tCross = t(k);
            else
                frac = (threshold - v1) / (v2 - v1);
                tCross = t(k) + frac * (t(k+1) - t(k));
            end
            dtFromNadir = tCross - tN;
            reason = 'OK';
            return;
        end
    end
end

function metricReasons = buildMetricDiagnostics(metricVals, fitOk, seriesMeta, doublingReason, ttpReason, idx0, idxNadir, gfRecovery, v0, nadirVal)
    metricReasons = repmat({'OK'}, 1, numel(metricVals));

    if ~isfinite(metricVals(1))
        if ~fitOk
            metricReasons{1} = 'Gompertz fit failed or insufficient post-nadir data';
        else
            metricReasons{1} = 'Alpha not finite';
        end
    end

    if ~isfinite(metricVals(2))
        metricReasons{2} = doublingReason;
    end

    if ~isfinite(metricVals(3))
        metricReasons{3} = 'Fractional reduction invalid (check V0 and nadir positivity)';
    end

    if ~isfinite(metricVals(4))
        metricReasons{4} = 'Log-kill invalid (peak/nadir non-positive or non-finite)';
    end

    if ~isfinite(metricVals(5))
        if ~seriesMeta.hasDeadSeries
            metricReasons{5} = 'Dead-cell column not found';
        else
            metricReasons{5} = 'Dead/total ratio unavailable at t0 or nadir';
        end
    end

    if ~isfinite(metricVals(6))
        if ~seriesMeta.hasProlifSeries
            metricReasons{6} = 'Proliferative-cell column not found';
        elseif isempty(gfRecovery)
            metricReasons{6} = 'No recovery points available for GF';
        else
            metricReasons{6} = 'GF contains no finite values in recovery phase';
        end
    end

    if ~isfinite(metricVals(7))
        metricReasons{7} = ttpReason;
    end

    if ~isfinite(metricVals(8))
        if idxNadir <= idx0
            metricReasons{8} = 'No dosing-phase interval before nadir';
        else
            metricReasons{8} = 'Regression slope unavailable (non-finite or insufficient points)';
        end
    end

    if ~isfinite(v0) || v0 <= 0
        metricReasons{3} = 'Baseline V0 is non-positive or non-finite';
        metricReasons{7} = 'Baseline V0 is non-positive or non-finite';
    end
    if ~isfinite(nadirVal) || nadirVal <= 0
        metricReasons{2} = 'Nadir is non-positive or non-finite';
        metricReasons{4} = 'Nadir is non-positive or non-finite';
    end
end

function txt = formatMetricValue(v)
    if ~isfinite(v)
        txt = 'NaN';
    elseif abs(v) >= 1e4 || (abs(v) > 0 && abs(v) < 1e-3)
        txt = sprintf('%.2e', v);
    else
        txt = sprintf('%.4g', v);
    end
end

function writeDiagnosticReport(diagTxtPath, csvPath, metricNames, metricVals, metricReasons, ...
    fitOk, Kfit, alphaFit, v0, nadirVal, peakVal, tNadir, doseStartTime, ...
    doublingTarget, ttpTarget, necFrac0, necFracNadir, live0, timeToLiveRegrowth, ...
    liveRegrowthReason, seriesMeta, timeVals, idx0, idxNadir, gfRecovery, figOut)

    fid = fopen(diagTxtPath, 'w');
    if fid < 0
        warning('Could not write diagnostics file: %s', diagTxtPath);
        return;
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, 'Post-Nadir Regrowth Diagnostics\n');
    fprintf(fid, '===========================================\n');
    fprintf(fid, 'Source CSV: %s\n', csvPath);
    fprintf(fid, 'Generated Figure: %s\n', ternary(~isempty(figOut), figOut, 'not saved'));
    fprintf(fid, 'Generated Diagnostics: %s\n\n', diagTxtPath);

    fprintf(fid, 'Data Integrity\n');
    fprintf(fid, '-------------------------------------------\n');
    fprintf(fid, 'Total valid time points: %d\n', numel(timeVals));
    fprintf(fid, 'Time range: [%.6g, %.6g]\n', timeVals(1), timeVals(end));
    fprintf(fid, 'Dose start time: %.6g\n', doseStartTime);
    fprintf(fid, 'Index at/after dose start: %d\n', idx0);
    fprintf(fid, 'Nadir index: %d\n', idxNadir);
    fprintf(fid, 'Nadir time: %.6g\n\n', tNadir);

    fprintf(fid, 'Column Mapping\n');
    fprintf(fid, '-------------------------------------------\n');
    fprintf(fid, 'time column: %s\n', seriesMeta.timeColumn);
    fprintf(fid, 'total column: %s\n', seriesMeta.totalColumn);
    fprintf(fid, 'live column: %s\n', valueOrMissing(seriesMeta.liveColumn));
    fprintf(fid, 'dead column: %s\n', valueOrMissing(seriesMeta.deadColumn));
    fprintf(fid, 'proliferative column: %s\n\n', valueOrMissing(seriesMeta.prolifColumn));

    fprintf(fid, 'Primary Evaluator\n');
    fprintf(fid, '-------------------------------------------\n');
    fprintf(fid, 'Baseline live cells (t0): %s\n', formatMetricValue(live0));
    if isfinite(timeToLiveRegrowth)
        fprintf(fid, 'Time to live-cell regrowth: %s\n\n', formatMetricValue(timeToLiveRegrowth));
    else
        fprintf(fid, 'Time to live-cell regrowth: NaN\n');
        fprintf(fid, 'reason: %s\n\n', liveRegrowthReason);
    end

    fprintf(fid, 'Core Scalars\n');
    fprintf(fid, '-------------------------------------------\n');
    fprintf(fid, 'V0 (at/after dose start): %s\n', formatMetricValue(v0));
    fprintf(fid, 'V_nadir: %s\n', formatMetricValue(nadirVal));
    fprintf(fid, 'V_peak: %s\n', formatMetricValue(peakVal));
    fprintf(fid, 'Doubling target (2 x V_nadir): %s\n', formatMetricValue(doublingTarget));
    fprintf(fid, 'TTP target (1.2 x V0): %s\n', formatMetricValue(ttpTarget));
    fprintf(fid, 'Necrotic fraction @V0: %s\n', formatMetricValue(necFrac0));
    fprintf(fid, 'Necrotic fraction @Nadir: %s\n\n', formatMetricValue(necFracNadir));

    fprintf(fid, 'Gompertz Fit\n');
    fprintf(fid, '-------------------------------------------\n');
    fprintf(fid, 'Fit status: %s\n', ternary(fitOk, 'OK', 'FAILED'));
    fprintf(fid, 'K_fit: %s\n', formatMetricValue(Kfit));
    fprintf(fid, 'alpha_fit: %s\n\n', formatMetricValue(alphaFit));

    fprintf(fid, 'Metric Results\n');
    fprintf(fid, '-------------------------------------------\n');
    for k = 1:numel(metricNames)
        if isfinite(metricVals(k))
            fprintf(fid, '%s: %s\n', metricNames{k}, formatMetricValue(metricVals(k)));
        else
            fprintf(fid, '%s: NaN\n', metricNames{k});
            fprintf(fid, '  reason: %s\n', metricReasons{k});
        end
    end
    fprintf(fid, '\n');

    fprintf(fid, 'Recovery GF Details\n');
    fprintf(fid, '-------------------------------------------\n');
    if isempty(gfRecovery)
        fprintf(fid, 'GF recovery points: 0\n');
    else
        finiteGf = gfRecovery(isfinite(gfRecovery));
        fprintf(fid, 'GF recovery points (total): %d\n', numel(gfRecovery));
        fprintf(fid, 'GF recovery points (finite): %d\n', numel(finiteGf));
        if ~isempty(finiteGf)
            fprintf(fid, 'GF min/max: [%s, %s]\n', formatMetricValue(min(finiteGf)), formatMetricValue(max(finiteGf)));
        end
    end

    clear cleaner;
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end

function out = valueOrMissing(s)
    if isempty(s)
        out = 'MISSING';
    else
        out = s;
    end
end

function colName = pickColumn(T, candidates)
    colName = '';
    vars = string(T.Properties.VariableNames);
    varsLower = lower(vars);

    for i = 1:numel(candidates)
        c = lower(string(candidates{i}));
        idx = find(varsLower == c, 1, 'first');
        if ~isempty(idx)
            colName = char(vars(idx));
            return;
        end
    end

    for i = 1:numel(candidates)
        c = lower(string(candidates{i}));
        idx = find(contains(varsLower, c), 1, 'first');
        if ~isempty(idx)
            colName = char(vars(idx));
            return;
        end
    end
end

function x = toNumeric(v)
    if isnumeric(v)
        x = double(v(:));
        return;
    end
    if iscell(v)
        x = str2double(string(v(:)));
        return;
    end
    x = str2double(string(v(:)));
end

function candidates = findMajorDipCandidates(timeVals, seriesVals, idxWindow, maxCandidates, prominenceFraction, minSeparationPoints)
    candidates = struct('idx', {}, 'time', {}, 'value', {}, 'prominence', {});

    if numel(idxWindow) < 3
        return;
    end

    wVals = seriesVals(idxWindow);
    wTime = timeVals(idxWindow);
    valid = isfinite(wVals);
    if ~any(valid)
        return;
    end

    valueRange = max(wVals(valid)) - min(wVals(valid));
    promThresh = max(1e-12, prominenceFraction * max(valueRange, 1));

    raw = struct('idx', {}, 'time', {}, 'value', {}, 'prominence', {});
    for k = 2:(numel(idxWindow)-1)
        vPrev = wVals(k-1);
        vCur = wVals(k);
        vNext = wVals(k+1);
        if ~isfinite(vPrev) || ~isfinite(vCur) || ~isfinite(vNext)
            continue;
        end
        if vCur <= vPrev && vCur <= vNext
            leftPeak = max(wVals(1:k));
            rightPeak = max(wVals(k:end));
            prom = min(leftPeak, rightPeak) - vCur;
            if isfinite(prom) && prom >= promThresh
                raw(end+1).idx = idxWindow(k); %#ok<AGROW>
                raw(end).time = wTime(k);
                raw(end).value = vCur;
                raw(end).prominence = prom;
            end
        end
    end

    if isempty(raw)
        return;
    end

    [~, order] = sort([raw.prominence], 'descend');
    raw = raw(order);

    selected = struct('idx', {}, 'time', {}, 'value', {}, 'prominence', {});
    for k = 1:numel(raw)
        if isempty(selected)
            selected(end+1) = raw(k); %#ok<AGROW>
        else
            sepOk = all(abs(raw(k).idx - [selected.idx]) >= minSeparationPoints);
            if sepOk
                selected(end+1) = raw(k); %#ok<AGROW>
            end
        end
        if numel(selected) >= maxCandidates
            break;
        end
    end

    [~, tOrder] = sort([selected.time], 'ascend');
    candidates = selected(tOrder);
end

function chosenIdx = chooseDipViaCLI(candidates, defaultIdx, csvPath, seriesLabel)
    chosenIdx = defaultIdx;

    fprintf('\nDip selector for %s\n', csvPath);
    fprintf('Series: %s\n', seriesLabel);

    if isempty(candidates)
        fprintf('No major dips detected. Using default nadir index %d.\n', defaultIdx);
        return;
    end

    fprintf('Major dip candidates:\n');
    for i = 1:numel(candidates)
        isDefault = candidates(i).idx == defaultIdx;
        marker = '';
        if isDefault
            marker = ' [default]';
        end
        fprintf('  %d) t=%s, value=%s, prominence=%s%s\n', ...
            i, formatMetricValue(candidates(i).time), formatMetricValue(candidates(i).value), ...
            formatMetricValue(candidates(i).prominence), marker);
    end

    in = input(sprintf('Select dip number [1-%d], Enter for default: ', numel(candidates)), 's');
    if isempty(strtrim(in))
        fprintf('Using default dip index %d.\n', defaultIdx);
        return;
    end

    sel = str2double(in);
    if isfinite(sel) && sel >= 1 && sel <= numel(candidates) && abs(sel - round(sel)) < 1e-12
        chosenIdx = candidates(round(sel)).idx;
        fprintf('Selected dip #%d (index %d).\n', round(sel), chosenIdx);
    else
        fprintf('Invalid selection. Using default dip index %d.\n', defaultIdx);
    end
end
