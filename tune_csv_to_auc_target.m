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
%   'DisablePlots'   : disable simulation plot export during tuning
%                      (default true; faster and avoids export failures)
%   'CleanupOnConverged': if true, keep only converged outputs and rename
%                        to clean base names without '_AUCtuned' (default true)
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
    p.addParameter('DisablePlots', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    p.addParameter('CleanupOnConverged', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    p.parse(inputCsv, varargin{:});

    targetAUC = p.Results.TargetAUC;
    tolerance = p.Results.Tolerance;
    maxIterations = round(p.Results.MaxIterations);
    damping = p.Results.Damping;
    maxStepFraction = p.Results.MaxStepFraction;
    disablePlots = logical(p.Results.DisablePlots);
    cleanupOnConverged = logical(p.Results.CleanupOnConverged);

    inputCsv = char(inputCsv);
    if ~isfile(inputCsv)
        error('Input CSV not found: %s', inputCsv);
    end

    [inFolder, inBase, ~] = fileparts(inputCsv);

    if strlength(string(p.Results.OutputPrefix)) == 0
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
    iterCsvPaths = cell(maxIterations,1);
    iterOutPrefixes = cell(maxIterations,1);

    for iter = 1:maxIterations
        iterationsRun = iter;

        iterCsv = sprintf('%s_iter%02d.csv', outputPrefix, iter);
        iterOutPrefix = sprintf('%s_sim%02d', outputPrefix, iter);
        iterCsvPaths{iter} = iterCsv;
        iterOutPrefixes{iter} = iterOutPrefix;
        writetable(currentData, iterCsv);

        try
            runOverrides = struct();
            if disablePlots
                runOverrides.generatePlots = false;
            end

            results = run5FU_PBPK_Simulation(iterCsv, iterOutPrefix, runOverrides);
            auc = results.metrics.AUC_central_mg_h_L;
        catch ME
            exportLikeFailure = contains(lower(ME.message), 'export failed') || ...
                                contains(lower(ME.message), 'problem while generating the output') || ...
                                contains(lower(ME.message), 'saveas');

            if exportLikeFailure && ~disablePlots
                retryOverrides = struct('generatePlots', false);
                retryOutPrefix = sprintf('%s_noplot', iterOutPrefix);
                try
                    results = run5FU_PBPK_Simulation(iterCsv, retryOutPrefix, retryOverrides);
                    auc = results.metrics.AUC_central_mg_h_L;
                    iterOutPrefixes{iter} = retryOutPrefix;
                    fprintf('Iter %2d | export failed with plots; retried with plotting disabled.\n', iter);
                catch MEretry
                    error('Simulation failed at iteration %d: %s', iter, MEretry.message);
                end
            else
                error('Simulation failed at iteration %d: %s', iter, ME.message);
            end
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

    if converged && cleanupOnConverged
        [optimizedCsvPath, cleanupSummary] = cleanupAndRenameConvergedArtifacts( ...
            inFolder, inBase, outputPrefix, iterationsRun, iterCsvPaths, iterOutPrefixes, optimizedCsvPath);
        summary.cleanup = cleanupSummary;
        summary.optimizedCsv = optimizedCsvPath;
    end

    if converged
        fprintf('\nConverged: final AUC %.4f is within ±%.4f of %.4f\n', lastAUC, tolerance, targetAUC);
    else
        fprintf('\nStopped at max iterations. Final AUC %.4f (error %+.4f)\n', lastAUC, lastErr);
    end
    fprintf('Optimized CSV saved to: %s\n', optimizedCsvPath);
end

function [finalCsvPath, cleanupSummary] = cleanupAndRenameConvergedArtifacts(inFolder, inBase, outputPrefix, convergedIter, iterCsvPaths, iterOutPrefixes, optimizedCsvPath)
    cleanupSummary = struct();
    cleanupSummary.deletedFiles = {};
    cleanupSummary.keptFiles = {};
    cleanupSummary.renamedFiles = {};

    finalCsvPath = fullfile(inFolder, sprintf('%s.csv', inBase));

    convergedIterCsv = iterCsvPaths{convergedIter};
    if isfile(convergedIterCsv)
        movefile(convergedIterCsv, finalCsvPath, 'f');
    elseif isfile(optimizedCsvPath)
        movefile(optimizedCsvPath, finalCsvPath, 'f');
    end

    if isfile(finalCsvPath)
        cleanupSummary.keptFiles{end+1} = finalCsvPath; %#ok<AGROW>
    end

    [~, outBase, ~] = fileparts(outputPrefix);
    convergedOutPrefix = iterOutPrefixes{convergedIter};
    [~, convergedOutBase, ~] = fileparts(convergedOutPrefix);

    convFiles = dir([convergedOutPrefix '*']);
    for i = 1:numel(convFiles)
        src = fullfile(convFiles(i).folder, convFiles(i).name);
        newName = regexprep(convFiles(i).name, ['^' regexptranslate('escape', convergedOutBase)], inBase);
        dst = fullfile(convFiles(i).folder, newName);
        if ~strcmp(src, dst)
            movefile(src, dst, 'f');
            cleanupSummary.renamedFiles{end+1} = struct('from', src, 'to', dst); %#ok<AGROW>
            cleanupSummary.keptFiles{end+1} = dst; %#ok<AGROW>
        else
            cleanupSummary.keptFiles{end+1} = src; %#ok<AGROW>
        end
    end

    % Delete all non-converged iteration CSVs and simulation artifacts
    for i = 1:numel(iterCsvPaths)
        if i == convergedIter
            continue;
        end
        if ~isempty(iterCsvPaths{i}) && isfile(iterCsvPaths{i})
            delete(iterCsvPaths{i});
            cleanupSummary.deletedFiles{end+1} = iterCsvPaths{i}; %#ok<AGROW>
        end
        if ~isempty(iterOutPrefixes{i})
            oldFiles = dir([iterOutPrefixes{i} '*']);
            for j = 1:numel(oldFiles)
                f = fullfile(oldFiles(j).folder, oldFiles(j).name);
                if isfile(f)
                    delete(f);
                    cleanupSummary.deletedFiles{end+1} = f; %#ok<AGROW>
                end
            end
        end
    end

    % Remove temporary optimized file if it is different from final path
    if isfile(optimizedCsvPath) && ~strcmp(optimizedCsvPath, finalCsvPath)
        delete(optimizedCsvPath);
        cleanupSummary.deletedFiles{end+1} = optimizedCsvPath; %#ok<AGROW>
    end

    % Remove any stray files with the AUC-tuned base prefix
    stray = dir(fullfile(inFolder, [outBase '*']));
    for i = 1:numel(stray)
        f = fullfile(stray(i).folder, stray(i).name);
        if isfile(f)
            delete(f);
            cleanupSummary.deletedFiles{end+1} = f; %#ok<AGROW>
        end
    end
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
