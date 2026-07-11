function [results] = run5FU_PBPK_Simulation(inputFile, outputPrefix, paramOverrides)
%run5FU_PBPK_Simulation  Simulate 5-FU PBPK/PD for a dosing regimen.
%   RESULTS = run5FU_PBPK_Simulation(INPUTFILE) integrates the circadian-
%   modulated, multi-compartment 5-Fluorouracil model for the dosing schedule
%   in INPUTFILE (CSV) and returns a struct of concentration time-series,
%   PK/PD metrics, mass balance and a dosing recommendation. CSV outputs and
%   (optionally) figures are written using OUTPUTPREFIX.
%
%   RESULTS = run5FU_PBPK_Simulation(INPUTFILE, OUTPUTPREFIX, PARAMOVERRIDES)
%   overrides default parameters (any field of initialise5FUParameters, e.g.
%   'BW', 'Vmax_DPD' [umol/min], 'solver_method', 'generatePlots') for
%   sensitivity / population runs.
%
%   Dosing CSV columns: start_time_min, end_time_min, dosing_type, plus a dose
%   column (dose_amount or dose_mg) or infusion_rate. Supported types: bolus,
%   continuous, sinusoidal, custom (see readDosingRegimen, calculateDosingRate).
%
%   Numerics: integrate5FU drives the shared RHS (odes5FU_rhs) with the stiff
%   solver ode15s by default; params.solver_method = 'fixed' selects a
%   fixed-step Euler reference instead.
%
%   Key outputs:
%     results.time_min / .time_hr        uniform time grid
%     results.concentrations.*           per-compartment series (uM)
%     results.metrics.AUC_central_mg_h_L primary exposure metric
%     results.mass_balance               closed-system accounting (umol)
%     results.recommendation             dose guidance string
%
%   References: Diasio 2001 (DPD kinetics), Harris 1990 (circadian DPD),
%   Konings 2010 (tumor partitioning), Schalhorn 1992 / Fang 2016 (clinical PK).
%
%   For RESEARCH & EDUCATIONAL use only; validate against clinical data before
%   any inference. See also integrate5FU, odes5FU_rhs, initialise5FUParameters.

    startup();  % put model subdirectories on the MATLAB path

    if nargin < 2 || isempty(outputPrefix)
        outputPrefix = '5FU_PBPK_output';
    end
    if nargin < 3 || isempty(paramOverrides)
        paramOverrides = struct();
    end

    % --- structured logger (overridable via paramOverrides.logging) ---
    if isfield(paramOverrides, 'logging'), logConfig = paramOverrides.logging; else, logConfig = struct(); end
    if ~isfield(logConfig, 'logDir'),      logConfig.logDir = fullfile(pwd, 'logs'); end
    if ~isfield(logConfig, 'level'),       logConfig.level = 'INFO'; end
    if ~isfield(logConfig, 'enableDebug'), logConfig.enableDebug = false; end
    logger = Logger(logConfig.logDir, logConfig.level, logConfig.enableDebug);
    logger.info('simulation_initializing', struct('inputFile', inputFile, 'outputPrefix', outputPrefix));

    % --- figure window policy (desktop vs headless/HPC) ---
    wants_plots = true;
    if isfield(paramOverrides, 'generatePlots'), wants_plots = logical(paramOverrides.generatePlots); end
    if wants_plots && usejava('desktop')
        set(0, 'DefaultFigureWindowStyle', 'docked'); set(0, 'DefaultFigureVisible', 'on');
    else
        set(0, 'DefaultFigureWindowStyle', 'normal'); set(0, 'DefaultFigureVisible', 'off');
    end

    % --- dosing regimen ---
    try
        dosingRegimen = readDosingRegimen(inputFile);
        logger.info('dosing_regimen_loaded', struct('rows', height(dosingRegimen)));
    catch ME
        logger.fatal('dosing_regimen_read_failed', ...
            struct('file', inputFile, 'message', ME.message, 'identifier', ME.identifier));
        rethrow(ME);
    end

    % --- parameters (+ overrides for Monte Carlo / sensitivity) ---
    params = initialise5FUParameters();
    override_fields = fieldnames(paramOverrides);
    for i = 1:numel(override_fields)
        f = override_fields{i};
        if isfield(params, f)
            params.(f) = paramOverrides.(f);
        end
    end
    validateClearanceParameters(params, logger);

    % --- integrate, package, emit ---
    try
        [time_min, concentrations] = integrate5FU(params, dosingRegimen, logger);
        results = packageResults(time_min, concentrations, params, dosingRegimen);

        saveDetailedCSVOutputs(results, outputPrefix, logger);
        if wants_plots
            generateComprehensivePlots(results, outputPrefix, logger);
        end
        printSummaryStatistics(results);
        logger.info('simulation_output_packaged', struct('outputPrefix', outputPrefix));
    catch ME
        logger.fatal('simulation_failed', ...
            struct('message', ME.message, 'identifier', ME.identifier, 'time', datetime("now")));
        rethrow(ME);
    end
end
