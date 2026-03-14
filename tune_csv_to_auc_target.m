function [optimizedCsvPath, summary] = tune_csv_to_auc_target(inputCsv, varargin)
%TUNE_CSV_TO_AUC_TARGET Iteratively scale dosing until target AUC is reached.
%   [optimizedCsvPath, summary] = tune_csv_to_auc_target(inputCsv)
%   [optimizedCsvPath, summary] = tune_csv_to_auc_target(inputCsv, 'Name', value, ...)
%
% Default behavior:
%   - Target central AUC: 25.0 mg·h/L
%   - Tolerance: ±0.1 mg·h/L
%   - Scales dose_amount and infusion_rate rows proportionally each iteration
%
% Name-value options:
%   'TargetAUC'      : desired central AUC (default 25.0)
%   'Tolerance'      : acceptable absolute error (default 0.1)
%   'MaxIterations'  : maximum optimization iterations (default 20)
%   'OutputPrefix'   : output prefix for temp/final files (default '<input>_AUCtuned')
%   'Damping'        : update damping 0-1 (default 0.75)
%   'MaxStepFraction': max relative scaling change per step (default 0.40)
%
% Example:
%   [csvOut, s] = tune_csv_to_auc_target('WP2/6_Cont.csv');

    p = inputParser;
    p.addRequired('inputCsv', @(x) ischar(x) || isstring(x));
    p.addParameter('TargetAUC', 25.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('Tolerance', 0.1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('MaxIterations', 20, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
    p.addParameter('OutputPrefix', '', @(x) ischar(x) || isstring(x));
    p.addParameter('Damping', 0.75, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 1);
    p.addParameter('MaxStepFraction', 0.40, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
    p.parse(inputCsv, varargin{:});

    targetAUC = p.Results.TargetAUC;
    tolerance = p.Results.Tolerance;
    maxIterations = round(p.Results.MaxIterations);
    damping = p.Results.Damping;
    maxStepFraction = p.Results.MaxStepFraction;

    inputCsv = char(inputCsv);
    if ~isfile(inputCsv)
        error('Input CSV not found: %s', inputCsv);
    end

    if strlength(string(p.Results.OutputPrefix)) == 0
        [inFolder, inBase, ~] = fileparts(inputCsv);
        outputPrefix = fullfile(inFolder, sprintf('%s_AUCtuned', inBase));
    else
        outputPrefix = char(p.Results.OutputPrefix);
    end

    originalData = readtable(inputCsv);
    currentData = originalData;

    fprintf('\n=== AUC Target Tuner ===\n');
    fprintf('Input: %s\n', inputCsv);
    fprintf('Target AUC: %.4f mg·h/L\n', targetAUC);
    fprintf('Tolerance: ±%.4f mg·h/L\n', tolerance);
    fprintf('Max iterations: %d\n\n', maxIterations);

    aucHistory = nan(maxIterations,1);
    errorHistory = nan(maxIterations,1);
    cumulativeScaleHistory = nan(maxIterations,1);
    perIterScaleHistory = nan(maxIterations,1);

    cumulativeScale = 1.0;
    converged = false;
    iterationsRun = 0;

    for iter = 1:maxIterations
        iterationsRun = iter;

        iterCsv = sprintf('%s_iter%02d.csv', outputPrefix, iter);
        iterOutPrefix = sprintf('%s_sim%02d', outputPrefix, iter);
        writetable(currentData, iterCsv);

        try
            results = run5FU_PBPK_Simulation(iterCsv, iterOutPrefix);
            auc = results.metrics.AUC_central_mg_h_L;
        catch ME
            error('Simulation failed at iteration %d: %s', iter, ME.message);
        end

        if ~isfinite(auc) || auc <= 0
            error('Invalid AUC at iteration %d: %.6g', iter, auc);
        end

        aucErr = auc - targetAUC;
        aucHistory(iter) = auc;
        errorHistory(iter) = aucErr;
        cumulativeScaleHistory(iter) = cumulativeScale;

        fprintf('Iter %2d | AUC=%.4f | error=%+.4f | cumulativeScale=%.6f\n', ...
            iter, auc, aucErr, cumulativeScale);

        if abs(aucErr) <= tolerance
            converged = true;
            break;
        end

        rawStepScale = targetAUC / auc;
        dampedStepScale = 1 + damping * (rawStepScale - 1);

        minScale = 1 - maxStepFraction;
        maxScale = 1 + maxStepFraction;
        stepScale = min(max(dampedStepScale, minScale), maxScale);

        currentData = applyScaling(currentData, stepScale);
        cumulativeScale = cumulativeScale * stepScale;
        perIterScaleHistory(iter) = stepScale;
    end

    optimizedCsvPath = sprintf('%s_optimized.csv', outputPrefix);
    writetable(currentData, optimizedCsvPath);

    lastAUC = aucHistory(iterationsRun);
    lastErr = errorHistory(iterationsRun);

    summary = struct();
    summary.inputCsv = inputCsv;
    summary.optimizedCsv = optimizedCsvPath;
    summary.targetAUC = targetAUC;
    summary.tolerance = tolerance;
    summary.converged = converged;
    summary.iterations = iterationsRun;
    summary.finalAUC = lastAUC;
    summary.finalError = lastErr;
    summary.cumulativeScale = cumulativeScale;
    summary.aucHistory = aucHistory(1:iterationsRun);
    summary.errorHistory = errorHistory(1:iterationsRun);
    summary.perIterationScale = perIterScaleHistory(1:iterationsRun);

    if converged
        fprintf('\nConverged: final AUC %.4f is within ±%.4f of %.4f\n', lastAUC, tolerance, targetAUC);
    else
        fprintf('\nStopped at max iterations. Final AUC %.4f (error %+.4f)\n', lastAUC, lastErr);
    end
    fprintf('Optimized CSV saved to: %s\n', optimizedCsvPath);
end

function T = applyScaling(T, scaleFactor)
    if ismember('dose_amount', T.Properties.VariableNames)
        doseVals = double(T.dose_amount);
        doseMask = isfinite(doseVals) & doseVals > 0;
        T.dose_amount(doseMask) = doseVals(doseMask) * scaleFactor;
    end

    if ismember('infusion_rate', T.Properties.VariableNames)
        rateVals = double(T.infusion_rate);
        rateMask = isfinite(rateVals) & rateVals > 0;
        T.infusion_rate(rateMask) = rateVals(rateMask) * scaleFactor;
    end

    if ismember('mean_rate', T.Properties.VariableNames)
        meanVals = double(T.mean_rate);
        meanMask = isfinite(meanVals) & meanVals > 0;
        T.mean_rate(meanMask) = meanVals(meanMask) * scaleFactor;
    end

    if ismember('amplitude', T.Properties.VariableNames)
        ampVals = double(T.amplitude);
        ampMask = isfinite(ampVals) & ampVals > 0;
        T.amplitude(ampMask) = ampVals(ampMask) * scaleFactor;
    end

    if ismember('rate_mg_per_min', T.Properties.VariableNames)
        directVals = double(T.rate_mg_per_min);
        directMask = isfinite(directVals) & directVals > 0;
        T.rate_mg_per_min(directMask) = directVals(directMask) * scaleFactor;
    end
end
