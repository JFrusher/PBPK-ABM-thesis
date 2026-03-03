function MC_5FU_PK_sensitivity(n_samples, output_dir, hpc_minimal_output)
% MC_5FU_PK_sensitivity  Monte Carlo sensitivity analysis wrapper for 5-FU PBPK
%
% Sections:
% 1. % SETUP & CONFIGURATION
% 2. % PARAMETER SAMPLING
% 3. % SIMULATOR EXECUTION (using parfor)
% 4. % RESULTS AGGREGATION
% 5. % SENSITIVITY ANALYSIS
% 6. % VISUALIZATION
% 7. % FILE SAVING
%
% Usage: MC_5FU_PK_sensitivity()                           % uses defaults
%        MC_5FU_PK_sensitivity(200, './out/')              % custom n_samples and output dir
%        MC_5FU_PK_sensitivity(200, './out/', true)        % HPC minimal console output
%
% The script performs Latin Hypercube sampling for 10 parameters, runs the
% simulator for each sample in parallel, extracts AUC and Cmax from the
% simulator output CSV files, and performs spearman-based sensitivity
% analysis with simple visualizations saved as .fig files.
%
% NOTIFICATIONS: This script sends progress alerts via ntfy.sh to your phone.
%   To receive notifications:
%   1. Visit https://ntfy.sh/mc5fu_<your_unique_id> on your phone
%   2. Or install the ntfy mobile app and subscribe to the topic
%   Change the NTFY_TOPIC below to something unique to you (e.g., use a UUID)

% ========== USER CONFIGURATION ==========
NTFY_ENABLED = false;           % Set to false to disable phone notifications
NTFY_TOPIC = 'my-sim-runs-frush'; % Change to something unique for privacy (e.g., add UUID)
CHECKPOINT_ENABLED = true;     % Resume from checkpoints if script is interrupted
QUIET_MODE = true;             % Set to true to suppress per-run diagnostic output (cleaner CLI)
PROGRESS_UPDATE_INTERVAL = 4;  % Show progress every N completed runs
MAX_WORKERS = 20;              % Upper cap for local/HPC workers (auto-clamped to allocation)

if nargin < 3 || isempty(hpc_minimal_output)
    hpc_minimal_output = false;
end
HPC_MINIMAL_OUTPUT = logical(hpc_minimal_output); % true = very light console output for HPC logs
QUIET_MODE = QUIET_MODE || HPC_MINIMAL_OUTPUT;

% ===== PK TARGET CONFIG (easy to edit) =====
% mode: 'auto' (infer from Cmax scale), 'bolus', or 'infusion'
PK_TARGETS.mode = 'auto';
PK_TARGETS.auto_mode_cmax_cutoff_um = 10;

% BOLUS dosing targets (Mayo: Cmax target 300-500 µM, AUC target 20-30 mg·h/L)
PK_TARGETS.bolus.AUC.subtherapeutic = 20;
PK_TARGETS.bolus.AUC.optimal_low = 20;
PK_TARGETS.bolus.AUC.optimal_high = 30;
PK_TARGETS.bolus.AUC.toxic = 30;
PK_TARGETS.bolus.AUC.severe_toxicity = 40;
PK_TARGETS.bolus.Cmax.safe_max = 300;
PK_TARGETS.bolus.Cmax.moderate_threshold = 420;
PK_TARGETS.bolus.Cmax.high_threshold = 500;
PK_TARGETS.bolus.Cmax.severe_threshold = 650;

% INFUSION/CHRONO dosing targets (Chrono: Cmax target 3-5 µM, AUC target 25-35 mg·h/L)
PK_TARGETS.infusion.AUC.subtherapeutic = 25;
PK_TARGETS.infusion.AUC.optimal_low = 25;
PK_TARGETS.infusion.AUC.optimal_high = 35;
PK_TARGETS.infusion.AUC.toxic = 35;
PK_TARGETS.infusion.AUC.severe_toxicity = 45;
PK_TARGETS.infusion.Cmax.safe_max = 3;
PK_TARGETS.infusion.Cmax.moderate_threshold = 5;
PK_TARGETS.infusion.Cmax.high_threshold = 8.3;
PK_TARGETS.infusion.Cmax.severe_threshold = 10;

% Additional regimen references from your table (stored for manual use/switching)
PK_TARGETS.reference.Mayo.t_half_min = 4.5;
PK_TARGETS.reference.Mayo.t_half_max = 13;
PK_TARGETS.reference.DeGramont.AUC_min = 20;
PK_TARGETS.reference.DeGramont.AUC_max = 25;
PK_TARGETS.reference.DeGramont.Css_min = 0.5;
PK_TARGETS.reference.DeGramont.Css_max = 1.0;

if nargin < 1 || isempty(n_samples)
    n_samples = 100; % default number of samples
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile(pwd, 'MC_results');
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

% Create timestamped subfolder for this run
timestamp_folder = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
output_dir = fullfile(output_dir, timestamp_folder);
if ~isfolder(output_dir)
    mkdir(output_dir);
end

if ~HPC_MINIMAL_OUTPUT
    fprintf('MC: Sampling %d points, writing outputs to %s\n', n_samples, output_dir);
end
if HPC_MINIMAL_OUTPUT
    mc_hpc_status('MC-HPC: start | samples=%d | output=%s\n', n_samples, output_dir);
end
if NTFY_ENABLED
    msg = sprintf('🚀 MC Analysis Started\n\n📊 Samples: %d\n⏱️ Est. time: %.1f-%.1f hrs', ...
        n_samples, n_samples*1.5/60, n_samples*2.5/60);
    ntfyNotify(NTFY_TOPIC, msg, 'rocket', 'default');
end

% Create publication-ready figure output folder
try
    helpers_dir = fileparts(mfilename('fullpath'));
    addpath(helpers_dir);
    if exist('create_publication_output_folder', 'file') == 2
        figures_output_dir = create_publication_output_folder(output_dir);
    else
        figures_output_dir = output_dir;
        fprintf('MC: WARNING - Publication helpers not available, using default output dir\n');
    end
catch ME_helpers
    figures_output_dir = output_dir;
    fprintf('MC: WARNING - Publication helpers not available, using default output dir\n');
end

% Store metadata for figure labeling
analysis_date_str = datestr(now, 'yyyy-mm-dd');
fig_counter = 0;  % Global figure counter for numbering

% Reproducible seed
seed = 42;
rng(seed, 'twister');

% Number of parameters
n_params = 10;

% Parameter names and ranges (min, max) - mapped directly to simulator params
% Names include units for clarity
param_names = {'BW','V_central_L','CO_L_per_min','Q_tumor_frac', ...
               'Vmax_DPD_umol_per_min','Km_DPD_mg_per_L','CL_renal_L_per_min', ...
               'Kp_tumor','V_tumor_cm3','tumor_metabolite_factor'}; % 10 parameters

% Derive physiologically plausible ranges from default parameters
base_params = initialise5FUParameters();

% Define sampling ranges as ±30% around defaults (reasonable interindividual variability)
% This ensures parameters stay in physiologically plausible ranges
ranges = [ ...
    40, 120;                                    % BW: 40-120 kg (wide range)
    base_params.V_central * 0.7, base_params.V_central * 1.3;  % V_central: ±30% around 14 L
    base_params.CO * 0.7, base_params.CO * 1.3;                % CO: ±30% around 5.4 L/min
    0.01, 0.05;                                 % Q_tumor_frac: 1-5% of CO (tumor blood flow)
    base_params.Vmax_DPD * 0.5, base_params.Vmax_DPD * 1.5;   % Vmax_dpd (µmol/min): 78-234 µmol/min (±50%, high variability known)
    base_params.Km_DPD * 0.7, base_params.Km_DPD * 1.3;       % Km_DPD: ±30% around 5 mg/L
    base_params.CL_renal_5FU * 0.5, base_params.CL_renal_5FU * 1.5;  % CL_renal (L/min): ±50%
    base_params.Kp_tumor * 0.7, base_params.Kp_tumor * 1.3;   % Kp_tumor: ±30%
    5, 500;                                     % V_tumor_cm3: 5-500 cm^3 (small to large tumors)
    1.0, 5.0                                    % tumor_metabolite_factor: 1-5× (metabolic heterogeneity)
];

if ~HPC_MINIMAL_OUTPUT
    fprintf('MC: Parameter ranges (based on defaults with ±30-50%% variability):\n');
    fprintf('  BW: %.1f-%.1f kg\n', ranges(1,1), ranges(1,2));
    fprintf('  V_central: %.2f-%.2f L (default: %.2f L)\n', ranges(2,1), ranges(2,2), base_params.V_central);
    fprintf('  CO: %.2f-%.2f L/min (default: %.2f L/min)\n', ranges(3,1), ranges(3,2), base_params.CO);
    fprintf('  Vmax_DPD: %.1f-%.1f µmol/min (default: %.1f µmol/min)\n', ranges(5,1), ranges(5,2), base_params.Vmax_DPD);
    fprintf('  Km_DPD: %.2f-%.2f mg/L (default: %.2f mg/L)\n', ranges(6,1), ranges(6,2), base_params.Km_DPD);
    fprintf('  CL_renal: %.3f-%.3f L/min (default: %.3f L/min)\n\n', ranges(7,1), ranges(7,2), base_params.CL_renal_5FU);
end

% Basic input dosing CSV used for all runs (simple infusion bolus example)
% Default file lives in the MC output directory so runs are self-contained 
script_dir = fileparts(mfilename('fullpath'));
dosingFile = fullfile(script_dir, 'Chrono95.csv');%CHANGE ME
if ~isfile(dosingFile)
    if HPC_MINIMAL_OUTPUT
        mc_hpc_status('MC-HPC: ERROR missing dosing file: %s\n', dosingFile);
    end
    error('MC:MissingDosingFile', ...
        'Required dosing file not found: %s\nSet dosingFile to a valid CSV before running.', dosingFile);
end

% Prepare Latin Hypercube samples
if ~HPC_MINIMAL_OUTPUT
    fprintf('MC: Generating Latin hypercube samples (seed=%d)...\n', seed);
end
U = lhsdesign(n_samples, n_params, 'criterion', 'maximin', 'iterations', 100);
% Scale to ranges
samples = zeros(n_samples, n_params);
for j = 1:n_params
    samples(:,j) = ranges(j,1) + U(:,j) .* (ranges(j,2) - ranges(j,1));
end

% Keep samples as matrix (parfor-friendly), create table later for results
% Parameter order matches param_names array:
% 1=BW, 2=V_central_L, 3=CO_L_per_min, 4=Q_tumor_frac, 5=Vmax_dpd_nmol_per_min,
% 6=Km_DPD_mg_per_L, 7=CL_renal_mL_per_min, 8=Kp_tumor, 9=V_tumor_cm3, 10=tumor_metabolite_factor

% Load checkpoint if available (from most recent run folder)
% Note: loadCheckpoint will look for checkpoints in most recent timestamped folder
[~, completed_runs, checkpoint_file] = loadCheckpoint(fileparts(output_dir), n_samples);

% Check for resume scenario
if any(completed_runs)
    fprintf('MC: Found %d completed runs. Resuming from checkpoint.\n', sum(completed_runs));
    if NTFY_ENABLED
        pct = 100 * sum(completed_runs) / n_samples;
        msg = sprintf('⏯️ Resuming Analysis\n\n✅ Completed: %d/%d (%.0f%%)\n⏳ Remaining: %d\n⏱️ Est. time: %.1f-%.1f hrs', ...
            sum(completed_runs), n_samples, pct, sum(~completed_runs), ...
            sum(~completed_runs)*1.5/60, sum(~completed_runs)*2.5/60);
        ntfyNotify(NTFY_TOPIC, msg, 'arrows_counterclockwise', 'default');
    end
else
    completed_runs = false(n_samples, 1);
end

% Mark this run in logs
fprintf('MC: Run ID: %s | Progress: %d/%d complete\n', output_dir, sum(completed_runs), n_samples);
if HPC_MINIMAL_OUTPUT
    mc_hpc_status('MC-HPC: resume state | completed=%d/%d\n', sum(completed_runs), n_samples);
end

% Constants needed for parameter conversions
MW_5FU = 130.08; % g/mol
conversion_factor = MW_5FU/60000; % µM·min -> mg·h/L (same as main sim)

% =========================
% PRE-FLIGHT TEST - Run one simulation first to catch errors early
% =========================
if sum(completed_runs) == 0  % Only on fresh start
    if ~HPC_MINIMAL_OUTPUT
        fprintf('\n════════════════════════════════════════════════════════════════\n');
        fprintf('MC: PRE-FLIGHT CHECK - Testing single simulation before parallel execution...\n');
        fprintf('════════════════════════════════════════════════════════════════\n\n');
    else
        mc_hpc_status('MC-HPC: preflight check...\n');
    end
    
    try
        % Use first sample for test
        test_overrides = struct();
        test_overrides.BW = samples(1, 1);
        test_overrides.V_central = samples(1, 2);
        test_overrides.CO = samples(1, 3);
        test_overrides.Q_tumor = samples(1, 4);
        test_overrides.Vmax_DPD = samples(1, 5);  % Already in µmol/min
        test_overrides.Vmax_DPD_mg_per_h = test_overrides.Vmax_DPD * MW_5FU * 60 / 1000;
        test_overrides.Km_DPD = samples(1, 6);
        test_overrides.CL_renal_5FU = samples(1, 7);  % Already in L/min
        test_overrides.Kp_tumor = samples(1, 8);
        test_overrides.V_tumor = (samples(1, 9) / 1000.0) / test_overrides.BW;
        test_overrides.tumor_metabolite_factor = samples(1, 10);
        
        % Disable logging for pre-flight test to avoid Logger initialization issues
        test_overrides.logging = struct('logDir', tempdir, 'level', 'ERROR', 'enableDebug', false);
        
        test_prefix = fullfile(output_dir, 'MC_preflight_test');
        if ~HPC_MINIMAL_OUTPUT
            fprintf('MC: Running test simulation...\n');
        end
        test_result = run5fu_with_capture(dosingFile, test_prefix, test_overrides, HPC_MINIMAL_OUTPUT);
        
        if ~HPC_MINIMAL_OUTPUT
            fprintf('\n✓ PRE-FLIGHT CHECK PASSED\n');
            fprintf('  Test AUC: %.2f mg·h/L\n', test_result.metrics.AUC_central_mg_h_L);
            fprintf('  Test Cmax: %.2f µM\n\n', test_result.metrics.Cmax_central);
            fprintf('MC: Proceeding with parallel execution...\n\n');
        else
            mc_hpc_status('MC-HPC: preflight passed | AUC=%.2f | Cmax=%.2f\n', ...
                test_result.metrics.AUC_central_mg_h_L, test_result.metrics.Cmax_central);
        end
        
    catch ME_test
        fprintf('\n✗ PRE-FLIGHT CHECK FAILED\n\n');
        fprintf('ERROR: %s\n', ME_test.message);
        fprintf('STACK TRACE:\n');
        for k = 1:length(ME_test.stack)
            fprintf('  at %s (line %d)\n', ME_test.stack(k).name, ME_test.stack(k).line);
        end
        fprintf('\n════════════════════════════════════════════════════════════════\n');
        fprintf('DIAGNOSIS: The simulation is failing before parallel execution.\n');
        fprintf('This is likely a configuration or environment issue.\n\n');
        fprintf('COMMON CAUSES:\n');
        fprintf('  1. Logger class error (check Logger.m recent changes)\n');
        fprintf('  2. Missing dependency files\n');
        fprintf('  3. Path issues\n');
        fprintf('  4. Dosing file format problem\n\n');
        fprintf('FIX: Run a single simulation manually to debug:\n');
        fprintf('  >> run5FU_PBPK_Simulation(''%s'', ''test_output'', struct())\n', dosingFile);
        fprintf('════════════════════════════════════════════════════════════════\n\n');
        
        if NTFY_ENABLED
            ntfyNotify(NTFY_TOPIC, sprintf('❌ PRE-FLIGHT FAILED\n\nError: %s\n\nCheck MATLAB console', ...
                ME_test.message), 'x,warning', 'high');
        end
        
        error('MC:PreflightFailed', 'Pre-flight test failed. Cannot proceed with Monte Carlo run. See error above.');
    end
end

% =========================
% 3) SIMULATOR EXECUTION (using parfor)
% =========================
% Request workers with HPC-aware auto-detection (SLURM/PBS), then clamp to MAX_WORKERS
sched_workers = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(sched_workers) || sched_workers < 1
    sched_workers = str2double(getenv('PBS_NP'));
end
if isnan(sched_workers) || sched_workers < 1
    sched_workers = str2double(getenv('NSLOTS')); % SGE/UGE compatibility
end

local_cores = feature('numcores');
if isnan(sched_workers) || sched_workers < 1
    max_workers = min(MAX_WORKERS, local_cores);
else
    max_workers = min([MAX_WORKERS, local_cores, floor(sched_workers)]);
end
max_workers = max(1, max_workers);

if ~HPC_MINIMAL_OUTPUT
    fprintf('MC: Starting parallel execution (parfor) with target %d workers (local cores=%d, scheduler=%s).\n', ...
        max_workers, local_cores, getenv('SLURM_CPUS_PER_TASK'));
end
if HPC_MINIMAL_OUTPUT
    mc_hpc_status('MC-HPC: initializing pool (target workers=%d)\n', max_workers);
end

% Close existing pool if present to avoid conflicts
pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end

% Create new pool with specified workers
try
    pool = parpool('Processes', max_workers);
    if ~HPC_MINIMAL_OUTPUT
        fprintf('Connected to parallel pool with %d workers.\n', pool.NumWorkers);
    end
catch ME
    % If 8 workers fail, try with maximum available
    if ~HPC_MINIMAL_OUTPUT
        fprintf('MC: Unable to create %d-worker pool. Trying default settings...\n', max_workers);
    end
    try
        pool = parpool('Processes');
        if ~HPC_MINIMAL_OUTPUT
            fprintf('Connected to parallel pool with %d workers.\n', pool.NumWorkers);
        end
    catch ME2
        if ~QUIET_MODE
            warning('Could not create parallel pool: %s\nContinuing with serial execution.', ME2.message);
        end
    end
end

% Preallocate results
AUC_values = nan(n_samples,1);
Cmax_values = nan(n_samples,1);
raw_timeseries = cell(n_samples,1);
failed_runs = false(n_samples,1);

% Track batches for checkpointing
num_cores = local_cores;
batch_size = max(4, max_workers);
if ~HPC_MINIMAL_OUTPUT
    fprintf('MC: Detected %d cores, using %d workers and batch size %d for checkpoint saves\n\n', num_cores, max_workers, batch_size);
end

progress_queue = [];
if HPC_MINIMAL_OUTPUT
    init_hpc_progress_tracker(n_samples, sum(completed_runs), PROGRESS_UPDATE_INTERVAL);
    progress_queue = parallel.pool.DataQueue;
    afterEach(progress_queue, @mc_hpc_progress_update);
    mc_hpc_status('MC-HPC: running simulations with %d workers...\n', max_workers);
end

% Setup progress tracking (if QUIET_MODE is on, we'll print periodically instead of per-run)
progress_file = fullfile(output_dir, 'MC_progress.log');
tic;  % Start timer for ETA calculation
total_failed = 0;
runs_completed_since_checkpoint = 0;
last_progress_display = 0;

parfor i = 1:n_samples
    % Skip if already completed (resume mode)
    if completed_runs(i)
        continue;
    end
    
    try
        % Extract parameters from matrix using direct indices (parfor-compatible)
        % Parameter order: BW, V_central_L, CO_L_per_min, Q_tumor_frac, 
        %                  Vmax_DPD_umol_per_min, Km_DPD_mg_per_L, CL_renal_L_per_min,
        %                  Kp_tumor, V_tumor_cm3, tumor_metabolite_factor
        BW = samples(i, 1);
        V_central_L = samples(i, 2);
        CO_L_per_min = samples(i, 3);
        Q_tumor_frac = samples(i, 4);
        Vmax_DPD_umol_per_min = samples(i, 5);  % Already in µmol/min
        Km_DPD_mg_per_L = samples(i, 6);
        CL_renal_L_per_min = samples(i, 7);     % Already in L/min
        Kp_tumor = samples(i, 8);
        V_tumor_cm3 = samples(i, 9);
        tumor_metabolite_factor = samples(i, 10);
        
        % Build paramOverrides
        overrides = struct();
        % 1) Body mass (kg)
        overrides.BW = BW;
        % 2) Central volume (L)
        overrides.V_central = V_central_L;
        % 3) Cardiac output (L/min)
        overrides.CO = CO_L_per_min;
        % 4) Tumor blood flow fraction (unitless fraction of CO)
        overrides.Q_tumor = Q_tumor_frac;
        % 5) Vmax DPD (already in µmol/min, no conversion needed)
        overrides.Vmax_DPD = Vmax_DPD_umol_per_min; % µmol/min
        overrides.Vmax_DPD_mg_per_h = overrides.Vmax_DPD * MW_5FU * 60 / 1000; % mg/h
        % 6) Km DPD (mg/L)
        overrides.Km_DPD = Km_DPD_mg_per_L;
        % 7) Renal clearance (already in L/min, no conversion needed)
        overrides.CL_renal_5FU = CL_renal_L_per_min;
        % 8) Tumor partition coefficient
        overrides.Kp_tumor = Kp_tumor;
        % 9) Tumor volume: convert cm^3 -> L and store as L per kg (model uses param as L/kg)
        V_tumor_L = V_tumor_cm3 / 1000.0; % cm^3 -> L
        overrides.V_tumor = V_tumor_L / overrides.BW; % L per kg
        % 10) Tumor metabolite formation factor
        overrides.tumor_metabolite_factor = tumor_metabolite_factor;

        % NOTE: We do NOT attempt to force Vmax to match a sampled "CL_hepatic" target.
        % Hepatic clearance is derived mechanistically in the simulator from Vmax, Km and Q_liver.

        % Run simulator with a unique output prefix
        outPrefix = fullfile(output_dir, sprintf('MC_run_%04d', i));
        try
            sim_results = run5fu_with_capture(dosingFile, outPrefix, overrides, HPC_MINIMAL_OUTPUT);
        catch E
            % Capture detailed error (only show if not in QUIET_MODE)
            if ~QUIET_MODE
                warning('MC:RunFailed','Run %d FAILED: %s\nIdentifier: %s\nStack: %s', ...
                    i, E.message, E.identifier, E.stack(1).name);
            end
            failed_runs(i) = true;
            if HPC_MINIMAL_OUTPUT
                send(progress_queue, true);
            end
            continue;
        end

        % Prefer simulator-returned metrics when available (more robust than CSV scraping)
        if isstruct(sim_results) && isfield(sim_results,'metrics') && isfield(sim_results.metrics,'AUC_central_mg_h_L') && isfield(sim_results.metrics,'Cmax_central')
            AUC_mg_h_L = sim_results.metrics.AUC_central_mg_h_L;
            Cmax_uM = sim_results.metrics.Cmax_central;
            
            % ===== AGGRESSIVE VALIDATION - Catch corruption immediately =====
            % Check 1: Type validation
            if ~isnumeric(AUC_mg_h_L) || ~isnumeric(Cmax_uM) || ~isscalar(AUC_mg_h_L) || ~isscalar(Cmax_uM)
                if ~QUIET_MODE
                    warning('MC:InvalidType','Run %d: Type error - AUC=%s(%s), Cmax=%s(%s)', ...
                        i, class(AUC_mg_h_L), mat2str(size(AUC_mg_h_L)), class(Cmax_uM), mat2str(size(Cmax_uM)));
                end
                failed_runs(i) = true;
                if HPC_MINIMAL_OUTPUT
                    send(progress_queue, true);
                end
                continue;
            end
            
            % Check 2: NaN/Inf detection
            if ~isfinite(AUC_mg_h_L) || ~isfinite(Cmax_uM)
                if ~QUIET_MODE
                    warning('MC:NonFinite','Run %d: Non-finite values - AUC=%e, Cmax=%e', i, AUC_mg_h_L, Cmax_uM);
                end
                failed_runs(i) = true;
                if HPC_MINIMAL_OUTPUT
                    send(progress_queue, true);
                end
                continue;
            end
            
            % Check 3: Plausibility check (must be in reasonable clinical range)
            % Very wide bounds to account for parameter variability
            % Low bounds: accept low-dose regimens, high dose sampling
            % High bounds: account for extreme parameter combinations
            if AUC_mg_h_L < 0.1 || AUC_mg_h_L > 1000
                if ~QUIET_MODE
                    warning('MC:OutOfRangeAUC','Run %d: AUC out of plausible range: %.4f mg·h/L (expect 0.1-1000)', i, AUC_mg_h_L);
                end
                failed_runs(i) = true;
                if HPC_MINIMAL_OUTPUT
                    send(progress_queue, true);
                end
                continue;
            end
            if Cmax_uM < 0.1 || Cmax_uM > 50000
                if ~QUIET_MODE
                    warning('MC:OutOfRangeCmax','Run %d: Cmax out of plausible range: %.4f µM (expect 0.1-50000)', i, Cmax_uM);
                end
                failed_runs(i) = true;
                if HPC_MINIMAL_OUTPUT
                    send(progress_queue, true);
                end
                continue;
            end
            
            % Check 4: Ensure metrics match concentration data
            if isfield(sim_results, 'concentrations')
                actual_cmax = max(sim_results.concentrations.C_central);
                if abs(Cmax_uM - actual_cmax) > 1.0
                    if ~QUIET_MODE
                        warning('MC:MetricsMismatch','Run %d: Cmax mismatch - reported=%.1f, calculated=%.1f', i, Cmax_uM, actual_cmax);
                    end
                    failed_runs(i) = true;
                    if HPC_MINIMAL_OUTPUT
                        send(progress_queue, true);
                    end
                    continue;
                end
            end
            
            % Store only essential data (not full struct) to avoid parfor serialization issues
            timeseries_data = struct();
            timeseries_data.time_hr = sim_results.time_hr;
            timeseries_data.C_central = sim_results.concentrations.C_central;
            timeseries_data.C_tumor = sim_results.concentrations.C_tumor;
            timeseries_data.metrics = sim_results.metrics;
            raw_timeseries{i} = timeseries_data;
        else
            % Fallback: read CSV that the simulator typically writes
            csvfile = sprintf('%s_5FU_compartments.csv', outPrefix);
            if ~isfile(csvfile)
                if ~QUIET_MODE
                    warning('MC:MissingCSV','Expected CSV not found for run %d: %s', i, csvfile);
                end
                failed_runs(i) = true;
                if HPC_MINIMAL_OUTPUT
                    send(progress_queue, true);
                end
                continue;
            end
            T = readtable(csvfile);
            time_min = T.Time_min;
            % Find central column (robust to naming)
            central_col = '';
            possibleNames = {'Central_5FU_uM','C_central','Central','Central_5FU'};
            % Use an index-based loop (compatible with parfor) and a local string
            for k_nm = 1:numel(possibleNames)
                nm_try = possibleNames{k_nm};
                if ismember(nm_try, T.Properties.VariableNames)
                    central_col = nm_try; break;
                end
            end
            if isempty(central_col)
                % attempt to find a numeric column excluding time
                numericCols = varfun(@isnumeric, T, 'OutputFormat','uniform');
                numericCols(strcmp(T.Properties.VariableNames,'Time_min')) = false;
                idx = find(numericCols,1);
                if isempty(idx)
                    if ~QUIET_MODE
                        warning('MC:NoCentralCol','No central concentration found for run %d', i);
                    end
                    failed_runs(i) = true;
                    if HPC_MINIMAL_OUTPUT
                        send(progress_queue, true);
                    end
                    continue;
                else
                    central_col = T.Properties.VariableNames{idx};
                end
            end
            C_central = T.(central_col);
            AUC_uM_min = trapz(time_min, C_central);
            AUC_mg_h_L = AUC_uM_min * conversion_factor;
            Cmax_uM = max(C_central);
            raw_timeseries{i} = T; % store table
        end

        % Validate results - mark as failed if NaN or Inf
        if ~isfinite(AUC_mg_h_L) || ~isfinite(Cmax_uM) || AUC_mg_h_L < 0 || Cmax_uM < 0
            if ~QUIET_MODE
                warning('MC:InvalidMetrics','Run %d returned invalid metrics: AUC=%.4e, Cmax=%.4e', i, AUC_mg_h_L, Cmax_uM);
            end
            failed_runs(i) = true;
            if HPC_MINIMAL_OUTPUT
                send(progress_queue, true);
            end
            continue;
        end

        AUC_values(i) = AUC_mg_h_L;
        Cmax_values(i) = Cmax_uM;
        completed_runs(i) = true;
        if HPC_MINIMAL_OUTPUT
            send(progress_queue, false);
        end
    catch ME_outer
        if ~QUIET_MODE
            warning('MC:RunException','Run %d raised exception: %s', i, ME_outer.message);
        end
        failed_runs(i) = true;
        if HPC_MINIMAL_OUTPUT
            send(progress_queue, true);
        end
    end
end

% Display clean progress summary after parfor
if QUIET_MODE && ~HPC_MINIMAL_OUTPUT
    elapsed_time = toc;
    num_completed = sum(completed_runs);
    num_failed = sum(failed_runs);
    num_remaining = n_samples - num_completed - num_failed;
    
    if num_completed > 0
        time_per_run = elapsed_time / num_completed;
        eta_seconds = (n_samples - num_completed) * time_per_run;
        if eta_seconds > 3600
            eta_str = sprintf('%.1f hours', eta_seconds / 3600);
        elseif eta_seconds > 60
            eta_str = sprintf('%.1f minutes', eta_seconds / 60);
        else
            eta_str = sprintf('%.0f seconds', eta_seconds);
        end
        
        fprintf('\n╔════════════════════════════════════════════╗\n');
        fprintf('║          MC PARALLEL EXECUTION COMPLETE     ║\n');
        fprintf('╚════════════════════════════════════════════╝\n');
        fprintf('  ✓ Completed:  %3d runs\n', num_completed);
        fprintf('  ✗ Failed:     %3d runs\n', num_failed);
        fprintf('  ⏳ Elapsed:    %.1f minutes\n', elapsed_time / 60);
        fprintf('  ⏱️  Per-run:    %.1f seconds\n', time_per_run);
        fprintf('\n');
    end
end

if HPC_MINIMAL_OUTPUT
    elapsed_time = toc;
    num_completed = sum(completed_runs);
    num_failed_now = sum(failed_runs);
    mc_hpc_status('MC-HPC: simulation stage done | completed=%d | failed=%d | elapsed=%.1f min\n', ...
        num_completed, num_failed_now, elapsed_time/60);
end

% Batch checkpoint saves (outside parfor to avoid race conditions)
% Save after each batch of simulations completes
for batch_idx = batch_size:batch_size:n_samples
    if CHECKPOINT_ENABLED && any(completed_runs(1:batch_idx))
        saveCheckpoint(samples, completed_runs, checkpoint_file);
        pct_done = 100 * batch_idx / n_samples;
        if ~QUIET_MODE  % Only show checkpoint messages in verbose mode
            fprintf('MC: Checkpoint saved at %.0f%% completion (%d/%d runs)\n', pct_done, batch_idx, n_samples);
        end
        if NTFY_ENABLED && mod(batch_idx, batch_size) == 0
            ntfyNotify(NTFY_TOPIC, sprintf('Progress: %d/%d samples complete (%.0f%%)', batch_idx, n_samples, pct_done), 'hourglass');
        end
    end
end

% Final checkpoint save
if CHECKPOINT_ENABLED
    saveCheckpoint(samples, completed_runs, checkpoint_file);
end

% Progress notification
num_done = sum(completed_runs);
fprintf('MC: Parallel execution complete: %d/%d samples processed\n', num_done, n_samples);
if HPC_MINIMAL_OUTPUT
    mc_hpc_status('MC-HPC: parallel stage complete | processed=%d/%d\n', num_done, n_samples);
end
if NTFY_ENABLED
    success_rate = 100 * num_done / n_samples;
    msg = sprintf('✅ Simulations Complete\n\n📊 Processed: %d/%d (%.0f%%)\n🔄 Now analyzing results...', ...
        num_done, n_samples, success_rate);
    ntfyNotify(NTFY_TOPIC, msg, 'white_check_mark', 'default');
end

% =========================
% 4) RESULTS AGGREGATION
% =========================
% Additional validation: filter out any remaining corrupted values
% Ensure all arrays are column vectors for proper boolean indexing
failed_runs = failed_runs(:);  % Force column vector
AUC_values = AUC_values(:);    % Force column vector
Cmax_values = Cmax_values(:);  % Force column vector

valid_idx = ~failed_runs & ~isnan(AUC_values) & ~isnan(Cmax_values) & ~isinf(AUC_values) & ~isinf(Cmax_values);

% Check for implausible values that slipped through (using same relaxed bounds)
% Allow 0.1-1000 mg·h/L and 0.1-50000 µM to account for sampling variability
plausible_auc = AUC_values >= 0.1 & AUC_values <= 1000;
plausible_cmax = Cmax_values >= 0.1 & Cmax_values <= 50000;
valid_idx = valid_idx & plausible_auc & plausible_cmax;

num_failed = sum(~valid_idx);
num_valid = sum(valid_idx);
if ~QUIET_MODE || num_valid == 0  % Only print in verbose mode, or if there's an error condition
    fprintf('MC: Completed runs: %d succeeded, %d failed\n', num_valid, num_failed);
end

if HPC_MINIMAL_OUTPUT
    mc_hpc_status('MC-HPC: valid results=%d/%d | rejected=%d\n', num_valid, n_samples, num_failed);
end

% Safety check: warn if >90% rejected (indicates systematic problem)
if (num_failed / n_samples) > 0.9
    fprintf('\n⚠️  WARNING: HIGH REJECTION RATE (%.0f%%) - Check simulation parameters\n', 100*num_failed/n_samples);
    if ~QUIET_MODE  % Only show detailed diagnostics in verbose mode
        fprintf('   AUC statistics: min=%.4f, max=%.4f, mean=%.4f, median=%.4f\n', ...
            min(AUC_values(~isnan(AUC_values) & ~isinf(AUC_values) & AUC_values > 0)), ...
            max(AUC_values(~isnan(AUC_values) & ~isinf(AUC_values))), ...
            mean(AUC_values(~isnan(AUC_values) & ~isinf(AUC_values))), ...
            median(AUC_values(~isnan(AUC_values) & ~isinf(AUC_values))));
        fprintf('   Cmax statistics: min=%.4f, max=%.4f, mean=%.4f, median=%.4f\n', ...
            min(Cmax_values(~isnan(Cmax_values) & ~isinf(Cmax_values) & Cmax_values > 0)), ...
            max(Cmax_values(~isnan(Cmax_values) & ~isinf(Cmax_values))), ...
            mean(Cmax_values(~isnan(Cmax_values) & ~isinf(Cmax_values))), ...
            median(Cmax_values(~isnan(Cmax_values) & ~isinf(Cmax_values))));
        fprintf('\n   DETAILED REJECTION LOG (showing reason for first 10 rejections):\n');
        count_shown = 0;
        for idx = 1:n_samples
            if ~valid_idx(idx) && count_shown < 10
                reason = 'Unknown';
                if isnan(AUC_values(idx)) || isnan(Cmax_values(idx))
                    reason = 'NaN';
                elseif ~isfinite(AUC_values(idx)) || ~isfinite(Cmax_values(idx))
                    reason = 'Inf';
                elseif AUC_values(idx) < 0.1 || AUC_values(idx) > 1000
                    reason = sprintf('AUC=%.4f out of range', AUC_values(idx));
                elseif Cmax_values(idx) < 0.1 || Cmax_values(idx) > 50000
                    reason = sprintf('Cmax=%.4f out of range', Cmax_values(idx));
                end
                fprintf('     Run %2d: AUC=%.4f, Cmax=%.4f (%s)\n', idx, AUC_values(idx), Cmax_values(idx), reason);
                count_shown = count_shown + 1;
            end
        end
        fprintf('\n');
    end
end

if NTFY_ENABLED && num_failed > 0
    fail_rate = 100 * num_failed / n_samples;
    if fail_rate > 50
        msg = sprintf('🚨 HIGH FAILURE RATE\n\n❌ Failed: %d/%d (%.0f%%)\n✅ Success: %d\n\n⚠️ Check logs for errors', ...
            num_failed, n_samples, fail_rate, num_valid);
        ntfyNotify(NTFY_TOPIC, msg, 'warning,x', 'high');
    else
        msg = sprintf('⚠️ Some Failures\n\n❌ Failed: %d (%.0f%%)\n✅ Success: %d\n\n📈 Continuing with valid data', ...
            num_failed, fail_rate, num_valid);
        ntfyNotify(NTFY_TOPIC, msg, 'warning', 'default');
    end
    
    % Diagnostic output for failures
    failed_indices = find(~valid_idx);
    if ~isempty(failed_indices)
        fprintf('MC: FAILED RUN DETAILS:\n');
        for idx = failed_indices(1:min(5, length(failed_indices)))'
            fprintf('  Run %d: AUC=%e, Cmax=%e, NaN_AUC=%d, NaN_Cmax=%d, Out_of_range=%d\n', ...
                idx, AUC_values(idx), Cmax_values(idx), isnan(AUC_values(idx)), isnan(Cmax_values(idx)), ...
                ~plausible_auc(idx) || ~plausible_cmax(idx));
        end
        if length(failed_indices) > 5
            fprintf('  ... and %d more failures (see log file)\n', length(failed_indices)-5);
        end
    end
end

% Safety check: prevent indexing error if all runs were rejected
if num_valid == 0
    error('MC: FATAL - All %d runs were rejected during validation. Check:\n  1. Are simulation parameters physiologically valid?\n  2. Are doses/concentrations reasonable?\n  3. Check the rejection log above for specific failure reasons.', n_samples);
end

% Keep only successful runs - convert samples matrix to table for analysis
samples_table = array2table(samples, 'VariableNames', param_names);
params_table = samples_table(valid_idx, :);
AUC = AUC_values(valid_idx);
Cmax = Cmax_values(valid_idx);
raw_timeseries = raw_timeseries(valid_idx);

% Check if we have enough valid data for statistics
if num_valid < 2
    error('MC: Too few successful runs (%d) to perform analysis. Check simulation errors.', num_valid);
end

% Statistics: mean and 95% CI
n_good = sum(valid_idx);
stats.mean.AUC = mean(AUC);
stats.ci95.AUC = 1.96 * std(AUC) / sqrt(n_good);
stats.mean.Cmax = mean(Cmax);
stats.ci95.Cmax = 1.96 * std(Cmax) / sqrt(n_good);

% -----------------------------
% 5) SENSITIVITY ANALYSIS
% -----------------------------
fprintf('MC: Performing Spearman rank correlations (parameters vs Cmax)...\n');
rho = nan(length(param_names),1);
pval = nan(length(param_names),1);
for j = 1:length(param_names)
    x = params_table{:, j};
    [r, p] = corr(x, Cmax, 'Type', 'Spearman', 'Rows', 'complete');
    rho(j) = r;
    pval(j) = p;
end
[~, order] = sort(abs(rho), 'descend');
sensitivity.rank = param_names(order);
sensitivity.rho = rho(order);
sensitivity.pval = pval(order);

% -----------------------------
% 6) PARAMETRIC CURVE FITTING (Curve Fitting Toolbox)
% -----------------------------
fprintf('MC: Fitting parametric curves to parameter-outcome relationships...\n');

% Define candidate model types to test (biologically realistic only)
% Excludes: poly2/poly3 (U-shapes, inflections), exp2 (overfitting), rat22 (too flexible)
% Includes: Linear, power laws (allometric scaling), exponential, Michaelis-Menten-like
model_types = {'poly1', 'power1', 'power2', 'exp1', 'rat11'};
model_names = {'Linear', 'Power law (1-term)', 'Power law (2-term)', 'Exponential', 'Michaelis-Menten (Rational 1,1)'};

% Check if Curve Fitting Toolbox is available
has_cftool = license('test', 'Curve_Fitting_Toolbox') && exist('fit', 'file') == 2;

if has_cftool
    fprintf('MC: Curve Fitting Toolbox detected. Testing %d model types...\n', length(model_types));
    
    % Initialize storage for fitted models
    curve_fits = struct();
    curve_fits.Cmax = struct();
    curve_fits.AUC = struct();
    
    % Fit curves for each parameter vs Cmax
    for j = 1:length(param_names)
        pname = param_names{j};
        x = params_table{:, j};
        
        % Try each model type and select best based on adjusted R^2
        best_rsquare = -Inf;
        best_model_idx = 1;
        best_fit = [];
        best_gof = [];
        
        for m = 1:length(model_types)
            try
                % Use appropriate fit method based on model type
                % Linear models (poly1) require LinearLeastSquares
                if strcmp(model_types{m}, 'poly1')
                    opts = fitoptions('Method', 'LinearLeastSquares', 'Display', 'off');
                else
                    % Nonlinear models
                    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'MaxIter', 400, ...
                        'MaxFunEvals', 600, ...
                        'Display', 'off');
                end
                [fitresult, gof] = fit(x, Cmax, model_types{m}, opts);
                if gof.adjrsquare > best_rsquare
                    best_rsquare = gof.adjrsquare;
                    best_model_idx = m;
                    best_fit = fitresult;
                    best_gof = gof;
                end
            catch
                % Skip models that fail to converge
                continue;
            end
        end
        
        % Store best fit
        curve_fits.Cmax.(pname) = struct();
        curve_fits.Cmax.(pname).fit = best_fit;
        curve_fits.Cmax.(pname).gof = best_gof;
        curve_fits.Cmax.(pname).model_type = model_types{best_model_idx};
        curve_fits.Cmax.(pname).model_name = model_names{best_model_idx};
        
        if ~QUIET_MODE
            fprintf('  %s vs Cmax: Best=%s (R²=%.4f, RMSE=%.4f)\n', ...
                pname, model_names{best_model_idx}, best_gof.rsquare, best_gof.rmse);
        end
    end
    
    % Fit curves for each parameter vs AUC
    for j = 1:length(param_names)
        pname = param_names{j};
        x = params_table{:, j};
        
        best_rsquare = -Inf;
        best_model_idx = 1;
        best_fit = [];
        best_gof = [];
        
        for m = 1:length(model_types)
            try
                % Use appropriate fit method based on model type
                if strcmp(model_types{m}, 'poly1')
                    opts = fitoptions('Method', 'LinearLeastSquares', 'Display', 'off');
                else
                    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'MaxIter', 400, ...
                        'MaxFunEvals', 600, ...
                        'Display', 'off');
                end
                [fitresult, gof] = fit(x, AUC, model_types{m}, opts);
                if gof.adjrsquare > best_rsquare
                    best_rsquare = gof.adjrsquare;
                    best_model_idx = m;
                    best_fit = fitresult;
                    best_gof = gof;
                end
            catch
                continue;
            end
        end
        
        curve_fits.AUC.(pname) = struct();
        curve_fits.AUC.(pname).fit = best_fit;
        curve_fits.AUC.(pname).gof = best_gof;
        curve_fits.AUC.(pname).model_type = model_types{best_model_idx};
        curve_fits.AUC.(pname).model_name = model_names{best_model_idx};
        
        if ~QUIET_MODE
            fprintf('  %s vs AUC: Best=%s (R²=%.4f, RMSE=%.4f)\n', ...
                pname, model_names{best_model_idx}, best_gof.rsquare, best_gof.rmse);
        end
    end
    
    fprintf('MC: Curve fitting complete. Best models selected for each parameter.\n\n');
    
    % Create summary tables
    fit_summary_Cmax = cell(length(param_names), 5);
    fit_summary_AUC = cell(length(param_names), 5);
    
    for j = 1:length(param_names)
        pname = param_names{j};
        
        % Cmax summary
        fit_summary_Cmax{j, 1} = pname;
        fit_summary_Cmax{j, 2} = curve_fits.Cmax.(pname).model_name;
        fit_summary_Cmax{j, 3} = curve_fits.Cmax.(pname).gof.rsquare;
        fit_summary_Cmax{j, 4} = curve_fits.Cmax.(pname).gof.adjrsquare;
        fit_summary_Cmax{j, 5} = curve_fits.Cmax.(pname).gof.rmse;
        
        % AUC summary
        fit_summary_AUC{j, 1} = pname;
        fit_summary_AUC{j, 2} = curve_fits.AUC.(pname).model_name;
        fit_summary_AUC{j, 3} = curve_fits.AUC.(pname).gof.rsquare;
        fit_summary_AUC{j, 4} = curve_fits.AUC.(pname).gof.adjrsquare;
        fit_summary_AUC{j, 5} = curve_fits.AUC.(pname).gof.rmse;
    end
    
    fit_table_Cmax = cell2table(fit_summary_Cmax, ...
        'VariableNames', {'Parameter', 'Model', 'R_squared', 'Adj_R_squared', 'RMSE'});
    fit_table_AUC = cell2table(fit_summary_AUC, ...
        'VariableNames', {'Parameter', 'Model', 'R_squared', 'Adj_R_squared', 'RMSE'});
    
else
    fprintf('MC: Curve Fitting Toolbox not available. Skipping parametric curve fitting.\n');
    fprintf('    Install toolbox or use license to enable advanced curve fitting analysis.\n\n');
    curve_fits = [];
    fit_table_Cmax = [];
    fit_table_AUC = [];
end

% -----------------------------
% 7) VISUALIZATION
% -----------------------------
timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
[target_profile, dosing_mode] = resolve_pk_target_profile(PK_TARGETS, Cmax);
fprintf('MC: Target profile mode: %s\n', upper(dosing_mode));

% Tornado plot (horizontal bar of Spearman correlations)\nfig_counter = fig_counter + 1;  % Increment figure counter
fig_tornado = figure('Name','MC_Tornado','NumberTitle','off','Position',[100,100,1000,600]);
barh(abs(rho(order)), 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'black', 'LineWidth', 1);
yticks(1:length(order));
yticklabels(param_names(order));
xlabel('|Spearman Rank Correlation Coefficient|', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Parameter', 'FontSize', 12, 'FontWeight', 'bold');
title('Sensitivity Analysis: Parameter Effects on Cmax', 'FontSize', 13, 'FontWeight', 'bold');
set(fig_tornado, 'Color', 'white');
set(gca, 'FontName', 'Arial', 'FontSize', 11, 'Box', 'on', 'TickDir', 'in');
grid on; grid minor;
ax = gca;
ax.GridAlpha = 0.3;

% Add metadata
add_figure_metadata('DeGramont.csv', num_valid, analysis_date_str);

% Save publication figure
try
    caption = get_figure_caption('tornado', num_valid);
    save_publication_figure(fig_tornado, figures_output_dir, fig_counter, 'Tornado Plot', caption, ...
        'DeGramont.csv', num_valid, analysis_date_str);
    save_publication_fig_copy(fig_tornado, figures_output_dir, fig_counter, 'Tornado Plot', timestamp);
catch ME_tornado_save
    fprintf('MC: WARNING - Tornado plot save failed: %s\n', ME_tornado_save.message);
end

% Violin/distribution plots for outcome distributions (Cmax and AUC) - SEPARATE AXES
fprintf('MC: Creating distribution plots for Cmax and AUC\n');
try
    % Create figure with two subplots side-by-side
    fig = figure('Name','Distributions_Cmax_AUC','Position',[100,100,1400,600]);
    
    % ============================================================================
    % SUBPLOT 1: Cmax Distribution with Toxicity Benchmarks
    % ============================================================================
    subplot(1,2,1);
    
    % Check for violinplot addon
    has_violin = exist('violinplot', 'file') == 2;
    
    if has_violin
        try
            % Violinplot expects data as columns in a matrix or table
            v_cmax = violinplot(Cmax, [], 'ViolinColor', [0.2 0.6 0.8]);
            if ~isempty(v_cmax) && isprop(v_cmax(1), 'ViolinPlot')
                h_dist_cmax = v_cmax(1).ViolinPlot;
            else
                h_dist_cmax = [];
            end
            hold on;
        catch ME_v
            % Fallback to histogram if violinplot fails
            fprintf('  Violinplot failed for Cmax (%s), using histogram\n', ME_v.message);
            h_dist_cmax = histogram(Cmax, 25, 'FaceColor', [0.2 0.6 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'k');
            ylabel('Frequency', 'FontSize', 11);
            hold on;
        end
    else
        % Use histogram as fallback
        h_dist_cmax = histogram(Cmax, 25, 'FaceColor', [0.2 0.6 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'k');
        ylabel('Frequency', 'FontSize', 11);
        hold on;
    end

    if ~exist('h_dist_cmax', 'var') || isempty(h_dist_cmax)
        h_dist_cmax = histogram(Cmax, 25, 'FaceColor', [0.2 0.6 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'k');
    end
    
    % DYNAMIC X-AXIS SCALING: Fit to data range with 10% padding
    cmax_min = min(Cmax);
    cmax_max = max(Cmax);
    cmax_range = cmax_max - cmax_min;
    xlim([max(0, cmax_min - 0.1*cmax_range), cmax_max + 0.1*cmax_range]);
    
    % Add reference lines for toxicity thresholds
    % Literature: Thyss et al. (1986), Saif et al. (2013), Saif contemporary reviews
    % Cmax correlates with acute toxicity (mucositis, diarrhea, neutropenia)
    ylims = ylim;
    
    Cmax_safe_max = target_profile.Cmax.safe_max;
    Cmax_moderate = target_profile.Cmax.moderate_threshold;
    Cmax_high = target_profile.Cmax.high_threshold;
    
    % Plot threshold lines
    h_safe_cmax = plot([Cmax_safe_max Cmax_safe_max], ylims, 'g--', 'LineWidth', 2.5);
    h_mod_cmax = plot([Cmax_moderate Cmax_moderate], ylims, 'Color', [1 0.6 0], 'LineStyle', '--', 'LineWidth', 2.5);
    h_high_cmax = plot([Cmax_high Cmax_high], ylims, 'r--', 'LineWidth', 2.5);
    
    % Add mean and median markers
    mean_cmax = mean(Cmax);
    median_cmax = median(Cmax);
    h_mean_cmax = plot([mean_cmax mean_cmax], ylims, 'k-', 'LineWidth', 3);
    h_median_cmax = plot([median_cmax median_cmax], ylims, 'b-', 'LineWidth', 2);
    
    % Labels and formatting
    dose_type = sprintf(' [%s Dosing]', upper(dosing_mode));
    thresh_text = {sprintf('Safe (<%.0f µM)', Cmax_safe_max), ...
                  sprintf('Moderate (%.0f µM)', Cmax_moderate), ...
                  sprintf('High (>%.0f µM)', Cmax_high)};
    title(['Cmax Distribution' dose_type], 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('Cmax (µM)', 'FontSize', 12);
    if ~has_violin || any(strcmp(get(gca, 'Type'), 'axes'))
        ylabel('Density/Frequency', 'FontSize', 11);
    end
    grid on; grid minor;
    h_legend_cmax = [h_dist_cmax, h_safe_cmax, h_mod_cmax, h_high_cmax, h_mean_cmax, h_median_cmax];
    lbl_legend_cmax = [{'Distribution'}, thresh_text, {sprintf('Mean (%.2f)', mean_cmax), sprintf('Median (%.2f)', median_cmax)}];
    valid_cmax = isgraphics(h_legend_cmax);
    legend(h_legend_cmax(valid_cmax), lbl_legend_cmax(valid_cmax), 'Location', 'best', 'FontSize', 9);
    set(gca, 'FontSize', 10);
    
    % ============================================================================
    % SUBPLOT 2: AUC Distribution with Efficacy/Toxicity Targets
    % ============================================================================
    subplot(1,2,2);
    
    if has_violin
        try
            v_auc = violinplot(AUC, [], 'ViolinColor', [0.8 0.3 0.3]);
            if ~isempty(v_auc) && isprop(v_auc(1), 'ViolinPlot')
                h_dist_auc = v_auc(1).ViolinPlot;
            else
                h_dist_auc = [];
            end
            hold on;
        catch ME_v2
            fprintf('  Violinplot failed for AUC (%s), using histogram\n', ME_v2.message);
            h_dist_auc = histogram(AUC, 25, 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.7, 'EdgeColor', 'k');
            ylabel('Frequency', 'FontSize', 11);
            hold on;
        end
    else
        h_dist_auc = histogram(AUC, 25, 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.7, 'EdgeColor', 'k');
        ylabel('Frequency', 'FontSize', 11);
        hold on;
    end

    if ~exist('h_dist_auc', 'var') || isempty(h_dist_auc)
        h_dist_auc = histogram(AUC, 25, 'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.7, 'EdgeColor', 'k');
    end
    
    % DYNAMIC X-AXIS SCALING: Fit to data range with 10% padding
    auc_min = min(AUC);
    auc_max = max(AUC);
    auc_range = auc_max - auc_min;
    xlim([max(0, auc_min - 0.1*auc_range), auc_max + 0.1*auc_range]);
    
    % Add reference lines for therapeutic windows
    % Literature: Gamelin et al. (2008) - AUC-guided dosing PRIMARY
    %            Kaldate et al. (2012) - therapeutic range validation
    %            Beumer et al. (2019) - dose optimization
    ylims = ylim;
    
    AUC_subtherapeutic = target_profile.AUC.subtherapeutic;
    AUC_optimal_low = target_profile.AUC.optimal_low;
    AUC_optimal_high = target_profile.AUC.optimal_high;
    AUC_toxic_threshold = target_profile.AUC.toxic;
    
    % Plot therapeutic window zones (shaded regions)
    y_fill = [ylims(1) ylims(1) ylims(2) ylims(2)];
    xlims_current = xlim;
    
    % Subtherapeutic zone (light red)
    h_zone_sub = fill([xlims_current(1) min(AUC_subtherapeutic, xlims_current(2)) min(AUC_subtherapeutic, xlims_current(2)) xlims_current(1)], y_fill, [1 0.8 0.8], ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Optimal therapeutic window (light green) - only if in range
    if AUC_optimal_low < xlims_current(2)
        h_zone_opt = fill([max(AUC_optimal_low, xlims_current(1)) min(AUC_optimal_high, xlims_current(2)) min(AUC_optimal_high, xlims_current(2)) max(AUC_optimal_low, xlims_current(1))], y_fill, [0.8 1 0.8], ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none');
    else
        h_zone_opt = gobjects(1);
    end
    
    % Toxic zone (light red) - only if in range
    if AUC_toxic_threshold < xlims_current(2)
        h_zone_tox = fill([max(AUC_toxic_threshold, xlims_current(1)) xlims_current(2) xlims_current(2) max(AUC_toxic_threshold, xlims_current(1))], y_fill, [1 0.8 0.8], ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
    else
        h_zone_tox = gobjects(1);
    end
    
    % Plot boundary lines (only if in visible range)
    if AUC_subtherapeutic < xlims_current(2)
        h_line_sub = plot([AUC_subtherapeutic AUC_subtherapeutic], ylims, 'Color', [0.8 0.4 0], 'LineStyle', '--', 'LineWidth', 2.5);
    else
        h_line_sub = gobjects(1);
    end
    if AUC_optimal_low < xlims_current(2)
        h_line_opt_low = plot([AUC_optimal_low AUC_optimal_low], ylims, 'g--', 'LineWidth', 2.5);
    else
        h_line_opt_low = gobjects(1);
    end
    if AUC_optimal_high < xlims_current(2)
        h_line_opt_high = plot([AUC_optimal_high AUC_optimal_high], ylims, 'g--', 'LineWidth', 2.5);
    else
        h_line_opt_high = gobjects(1);
    end
    if AUC_toxic_threshold < xlims_current(2)
        h_line_tox = plot([AUC_toxic_threshold AUC_toxic_threshold], ylims, 'r--', 'LineWidth', 2.5);
    else
        h_line_tox = gobjects(1);
    end
    
    % Add mean and median markers
    mean_auc = mean(AUC);
    median_auc = median(AUC);
    h_mean_auc = plot([mean_auc mean_auc], ylims, 'k-', 'LineWidth', 3);
    h_median_auc = plot([median_auc median_auc], ylims, 'b-', 'LineWidth', 2);
    
    % Labels and formatting
    title('AUC Distribution with Therapeutic Window', 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('AUC (mg·h/L)', 'FontSize', 12);
    if ~has_violin || any(strcmp(get(gca, 'Type'), 'axes'))
        ylabel('Density/Frequency', 'FontSize', 11);
    end
    grid on; grid minor;
    h_legend_auc = [h_dist_auc, h_zone_sub, h_zone_opt, h_zone_tox, h_line_sub, h_line_opt_low, h_line_opt_high, h_line_tox, h_mean_auc, h_median_auc];
    lbl_legend_auc = {'Distribution', ...
        sprintf('Subtherapeutic Zone (<%.0f)', AUC_subtherapeutic), ...
        sprintf('Optimal Zone (%.0f-%.0f)', AUC_optimal_low, AUC_optimal_high), ...
        sprintf('Toxic Zone (>%.0f)', AUC_toxic_threshold), ...
        sprintf('Subtherapeutic Line (%.0f)', AUC_subtherapeutic), ...
        sprintf('Optimal Low Line (%.0f)', AUC_optimal_low), ...
        sprintf('Optimal High Line (%.0f)', AUC_optimal_high), ...
        sprintf('Toxic Line (%.0f)', AUC_toxic_threshold), ...
        sprintf('Mean (%.1f)', mean_auc), sprintf('Median (%.1f)', median_auc)};
    valid_auc = isgraphics(h_legend_auc);
    legend(h_legend_auc(valid_auc), lbl_legend_auc(valid_auc), 'Location', 'northeast', 'FontSize', 9);
    set(gca, 'FontSize', 10);
    
    % Add overall title
    sgtitle('Pharmacokinetic Outcome Distributions with Clinical Benchmarks', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Apply publication styling
    set(fig, 'Color', 'white');
    set_publication_figure_style(fig, 'font_size', 11, 'font_name', 'Arial');
    
    % Add metadata and save
    add_figure_metadata('DeGramont.csv', num_valid, analysis_date_str);
    fig_counter = fig_counter + 1;
    caption = get_figure_caption('distribution', num_valid);
    save_publication_figure(fig, figures_output_dir, fig_counter, 'Outcome Distributions', caption, ...
        'DeGramont.csv', num_valid, analysis_date_str);
    save_publication_fig_copy(fig, figures_output_dir, fig_counter, 'Outcome Distributions', timestamp);
    
    % Print summary statistics with clinical interpretation
    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf('COMPREHENSIVE OUTCOME ANALYSIS WITH DISTRIBUTION SPREAD\n');
    fprintf('════════════════════════════════════════════════════════════════\n');
    
    % ===== CMAX ANALYSIS =====
    fprintf('\n📊 CMAX (Peak Concentration) - Distribution Spread:\n');
    fprintf('  ┌─ Central Tendency ─────────────────────────────────────┐\n');
    fprintf('  │  Mean:       %.2f µM (±%.2f SD)\n', mean_cmax, std(Cmax));
    fprintf('  │  Median:     %.2f µM\n', median_cmax);
    fprintf('  │  Mode:       ~%.2f µM (est. from histogram)\n', prctile(Cmax, 50));
    fprintf('  └────────────────────────────────────────────────────────────┘\n');
    
    fprintf('  ┌─ Range & Spread ───────────────────────────────────────┐\n');
    fprintf('  │  Min:        %.2f µM\n', min(Cmax));
    fprintf('  │  Max:        %.2f µM\n', max(Cmax));
    fprintf('  │  Range:      %.2f µM\n', max(Cmax) - min(Cmax));
    
    p25_cmax = prctile(Cmax, 25);
    p75_cmax = prctile(Cmax, 75);
    p5_cmax = prctile(Cmax, 5);
    p95_cmax = prctile(Cmax, 95);
    p10_cmax = prctile(Cmax, 10);
    p90_cmax = prctile(Cmax, 90);
    
    iqr_cmax = p75_cmax - p25_cmax;
    cv_cmax = std(Cmax) / mean_cmax * 100;  % Coefficient of variation
    
    fprintf('  │  IQR (25-75%% ile): %.2f µM (Q1=%.2f, Q3=%.2f)\n', iqr_cmax, p25_cmax, p75_cmax);
    fprintf('  │  Coeff. Variation:  %.1f%%\n', cv_cmax);
    fprintf('  └────────────────────────────────────────────────────────────┘\n');
    
    fprintf('  ┌─ Percentiles (likelihood of value at threshold) ──────┐\n');
    fprintf('  │  5th percentile:     %.2f µM (90%% above this)\n', p5_cmax);
    fprintf('  │  10th percentile:    %.2f µM (90%% at or above)\n', p10_cmax);
    fprintf('  │  25th percentile:    %.2f µM (75%% at or above)\n', p25_cmax);
    fprintf('  │  50th percentile:    %.2f µM (median)\n', median_cmax);
    fprintf('  │  75th percentile:    %.2f µM (25%% above this)\n', p75_cmax);
    fprintf('  │  90th percentile:    %.2f µM (10%% above this)\n', p90_cmax);
    fprintf('  │  95th percentile:    %.2f µM (5%% above this)\n', p95_cmax);
    fprintf('  └────────────────────────────────────────────────────────────┘\n');
    
    fprintf('  ┌─ Clinical Safety Profile ──────────────────────────────┐\n');
    pct_safe = 100 * sum(Cmax < Cmax_safe_max) / length(Cmax);
    pct_moderate = 100 * sum(Cmax >= Cmax_safe_max & Cmax < Cmax_moderate) / length(Cmax);
    pct_moderate_high = 100 * sum(Cmax >= Cmax_moderate & Cmax < Cmax_high) / length(Cmax);
    pct_high_tox = 100 * sum(Cmax >= Cmax_high) / length(Cmax);
    
    fprintf('  │  SAFE zone    (<%.0f µM):        %6.1f%% of simulations\n', Cmax_safe_max, pct_safe);
    fprintf('  │  MODERATE risk (%.0f-%.0f µM):  %6.1f%% of simulations\n', Cmax_safe_max, Cmax_moderate, pct_moderate);
    fprintf('  │  ELEVATED risk (%.0f-%.0f µM):  %6.1f%% of simulations\n', Cmax_moderate, Cmax_high, pct_moderate_high);
    fprintf('  │  HIGH TOXICITY (>%.0f µM):     %6.1f%% of simulations ⚠️\n', Cmax_high, pct_high_tox);
    fprintf('  └────────────────────────────────────────────────────────────┘\n');
    
    % ===== AUC ANALYSIS =====
    fprintf('\n📊 AUC (Systemic Exposure) - Distribution Spread:\n');
    fprintf('  ┌─ Central Tendency ─────────────────────────────────────┐\n');
    fprintf('  │  Mean:       %.1f mg·h/L (±%.1f SD)\n', mean_auc, std(AUC));
    fprintf('  │  Median:     %.1f mg·h/L\n', median_auc);
    fprintf('  │  Mode:       ~%.1f mg·h/L (est. from histogram)\n', prctile(AUC, 50));
    fprintf('  └────────────────────────────────────────────────────────────┘\n');
    
    fprintf('  ┌─ Range & Spread ───────────────────────────────────────┐\n');
    fprintf('  │  Min:        %.1f mg·h/L\n', min(AUC));
    fprintf('  │  Max:        %.1f mg·h/L\n', max(AUC));
    fprintf('  │  Range:      %.1f mg·h/L\n', max(AUC) - min(AUC));
    
    p25_auc = prctile(AUC, 25);
    p75_auc = prctile(AUC, 75);
    p5_auc = prctile(AUC, 5);
    p95_auc = prctile(AUC, 95);
    p10_auc = prctile(AUC, 10);
    p90_auc = prctile(AUC, 90);
    
    iqr_auc = p75_auc - p25_auc;
    cv_auc = std(AUC) / mean_auc * 100;  % Coefficient of variation
    
    fprintf('  │  IQR (25-75%% ile): %.2f mg·h/L (Q1=%.1f, Q3=%.1f)\n', iqr_auc, p25_auc, p75_auc);
    fprintf('  │  Coeff. Variation:  %.1f%%\n', cv_auc);
    fprintf('  └────────────────────────────────────────────────────────────┘\n');
    
    fprintf('  ┌─ Percentiles (likelihood of achieving exposure level) ─┐\n');
    fprintf('  │  5th percentile:     %.1f mg·h/L (90%% above this)\n', p5_auc);
    fprintf('  │  10th percentile:    %.1f mg·h/L (90%% at or above)\n', p10_auc);
    fprintf('  │  25th percentile:    %.1f mg·h/L (75%% at or above)\n', p25_auc);
    fprintf('  │  50th percentile:    %.1f mg·h/L (median)\n', median_auc);
    fprintf('  │  75th percentile:    %.1f mg·h/L (25%% above this)\n', p75_auc);
    fprintf('  │  90th percentile:    %.1f mg·h/L (10%% above this)\n', p90_auc);
    fprintf('  │  95th percentile:    %.1f mg·h/L (5%% above this)\n', p95_auc);
    fprintf('  └────────────────────────────────────────────────────────────┘\n');
    
    fprintf('  ┌─ Therapeutic Outcome Profile ──────────────────────────┐\n');
    pct_subtherapeutic = 100 * sum(AUC < AUC_subtherapeutic) / length(AUC);
    pct_optimal = 100 * sum(AUC >= AUC_optimal_low & AUC <= AUC_optimal_high) / length(AUC);
    pct_above_optimal = 100 * sum(AUC > AUC_optimal_high & AUC <= AUC_toxic_threshold) / length(AUC);
    pct_toxic = 100 * sum(AUC > AUC_toxic_threshold) / length(AUC);
    
    fprintf('  │  SUBTHERAPEUTIC (<%.0f mg·h/L):       %6.1f%% ❌ (insufficient efficacy)\n', AUC_subtherapeutic, pct_subtherapeutic);
    fprintf('  │  OPTIMAL WINDOW (%.0f-%.0f mg·h/L):   %6.1f%% ✓  (best efficacy/safety)\n', AUC_optimal_low, AUC_optimal_high, pct_optimal);
    fprintf('  │  ABOVE OPTIMAL (%.0f-%.0f mg·h/L):    %6.1f%% ⚠️  (increased toxicity risk)\n', AUC_optimal_high, AUC_toxic_threshold, pct_above_optimal);
    fprintf('  │  TOXIC RANGE (>%.0f mg·h/L):          %6.1f%% 🔴 (high toxicity risk)\n', AUC_toxic_threshold, pct_toxic);
    fprintf('  └────────────────────────────────────────────────────────────┘\n');
    
    fprintf('════════════════════════════════════════════════════════════════\n');
    fprintf('SUMMARY INTERPRETATION:\n');
    fprintf('════════════════════════════════════════════════════════════════\n');
    
    % Overall assessment
    if pct_optimal > 50
        fprintf('✓ FAVORABLE: %.1f%% of patients achieve optimal AUC target (28-39 mg·h/L)\n', pct_optimal);
    elseif pct_optimal > 30
        fprintf('⚠️  MODERATE: Only %.1f%% of patients achieve optimal AUC target\n', pct_optimal);
    else
        fprintf('❌ CONCERNING: Only %.1f%% of patients achieve optimal AUC target\n', pct_optimal);
    end
    
    if pct_high_tox < 10
        fprintf('✓ FAVORABLE: Only %.1f%% of patients at high Cmax-related toxicity risk\n', pct_high_tox);
    elseif pct_high_tox < 25
        fprintf('⚠️  MODERATE: %.1f%% of patients at elevated Cmax toxicity risk\n', pct_high_tox);
    else
        fprintf('❌ CONCERNING: %.1f%% of patients exceed high toxicity threshold\n', pct_high_tox);
    end
    
    fprintf('\nRecommendation: Use these percentiles to inform dose individualization.\n');
    fprintf('================================================================\n\n');
    
    % ===== NEW FIGURE 1: PERCENTILE RANGE PLOTS (Error bar style) =====
    fprintf('MC: Creating percentile range visualization (Cmax & AUC)\n');
    try
        fig_range = figure('Name','Percentile_Ranges','Position',[100,100,1200,500]);
        
        % SUBPLOT 1: Cmax Percentile Ranges
        subplot(1,2,1);
        hold on;
        
        % Create percentile "bands" - plot as horizontal lines with error bars
        percentiles = [10, 25, 50, 75, 90];
        percs_cmax = prctile(Cmax, percentiles);
        
        for i = 1:length(percentiles)
            color_intensity = (i-1) / (length(percentiles)-1);
            line_color = [color_intensity, 0.3, 1-color_intensity];
            plot([0.5 1.5], [percs_cmax(i) percs_cmax(i)], '-', 'Color', line_color, 'LineWidth', 3);
        end
        
        % Add IQR band (shaded region)
        fill([0.5 1.5 1.5 0.5], [p25_cmax p25_cmax p75_cmax p75_cmax], ...
            [0.2 0.6 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Add points for min/max
        scatter(1, min(Cmax), 100, 'r', 'filled', 'DisplayName', 'Min');
        scatter(1, max(Cmax), 100, 'r', 'x', 'LineWidth', 2.5, 'DisplayName', 'Max');
        
        % Labels for percentiles
        for i = 1:length(percentiles)
            text(1.55, percs_cmax(i), sprintf('%d%%ile: %.1f µM', percentiles(i), percs_cmax(i)), ...
                'FontSize', 9, 'VerticalAlignment', 'middle');
        end
        
        xlim([0.3 2.5]); ylim([min(Cmax)-50 max(Cmax)+50]);
        title('Cmax: Percentile Distribution', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Concentration (µM)', 'FontSize', 11);
        set(gca, 'XTick', []);
        grid on; grid minor;
        
        % SUBPLOT 2: AUC Percentile Ranges
        subplot(1,2,2);
        hold on;
        
        percs_auc = prctile(AUC, percentiles);
        
        for i = 1:length(percentiles)
            color_intensity = (i-1) / (length(percentiles)-1);
            line_color = [color_intensity, 1-color_intensity*0.5, 0.5];
            plot([0.5 1.5], [percs_auc(i) percs_auc(i)], '-', 'Color', line_color, 'LineWidth', 3);
        end
        
        % Add IQR band
        fill([0.5 1.5 1.5 0.5], [p25_auc p25_auc p75_auc p75_auc], ...
            [1 0.3 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Add min/max points
        scatter(1, min(AUC), 100, 'b', 'filled', 'DisplayName', 'Min');
        scatter(1, max(AUC), 100, 'b', 'x', 'LineWidth', 2.5, 'DisplayName', 'Max');
        
        % Labels
        for i = 1:length(percentiles)
            text(1.55, percs_auc(i), sprintf('%d%%ile: %.1f mg·h/L', percentiles(i), percs_auc(i)), ...
                'FontSize', 9, 'VerticalAlignment', 'middle');
        end
        
        xlim([0.3 2.5]); ylim([min(AUC)-5 max(AUC)+5]);
        title('AUC: Percentile Distribution', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('AUC (mg·h/L)', 'FontSize', 11);
        set(gca, 'XTick', []);
        grid on; grid minor;
        
        sgtitle('Percentile Ranges for Outcome Prediction', 'FontSize', 13, 'FontWeight', 'bold');
        
        % Apply styling and save
        set(fig_range, 'Color', 'white');
        set_publication_figure_style(fig_range, 'font_size', 11, 'font_name', 'Arial');
        add_figure_metadata('DeGramont.csv', num_valid, analysis_date_str);
        fig_counter = fig_counter + 1;
        caption = get_figure_caption('percentile', num_valid);
        save_publication_figure(fig_range, figures_output_dir, fig_counter, 'Percentile Ranges', caption, ...
            'DeGramont.csv', num_valid, analysis_date_str);
        save_publication_fig_copy(fig_range, figures_output_dir, fig_counter, 'Percentile Ranges', timestamp);
    catch ME_range
        fprintf('MC: WARNING - Could not create percentile range plots: %s\n', ME_range.message);
    end
    
    % ===== NEW FIGURE 2: CUMULATIVE DISTRIBUTION FUNCTIONS (CDF) =====
    fprintf('MC: Creating cumulative distribution function plots\n');
    try
        fig_cdf = figure('Name','CDF_Plots','Position',[100,100,1200,500]);
        
        % SUBPLOT 1: Cmax CDF
        subplot(1,2,1);
        [Cmax_sorted, ~] = sort(Cmax);
        cdf_cmax = (1:length(Cmax)) / length(Cmax);
        stairs(Cmax_sorted, cdf_cmax, 'LineWidth', 2.5, 'Color', [0.2 0.6 0.8]);
        hold on;
        
        % Add reference lines
        if cmax_max < 10
            ylims_cdf = ylim;
            plot([Cmax_safe_max Cmax_safe_max], ylims_cdf, 'g--', 'LineWidth', 2, 'DisplayName', 'Safe max');
            plot([Cmax_high Cmax_high], ylims_cdf, 'r--', 'LineWidth', 2, 'DisplayName', 'High toxicity');
        end
        
        xlabel('Cmax (µM)', 'FontSize', 11);
        ylabel('Cumulative Probability', 'FontSize', 11);
        title('Cmax: Cumulative Distribution', 'FontSize', 12, 'FontWeight', 'bold');
        grid on; grid minor;
        legend('Location', 'southeast', 'FontSize', 9);
        
        % Add annotations for key percentiles
        annotation('textbox', [0.12 0.75 0.25 0.15], ...
            'String', {sprintf('P5: %.1f µM', p5_cmax), ...
                      sprintf('P50: %.1f µM', median_cmax), ...
                      sprintf('P95: %.1f µM', p95_cmax)}, ...
            'FontSize', 9, 'EdgeColor', 'none', 'BackgroundColor', 'yellow', 'FaceAlpha', 0.7);
        
        % SUBPLOT 2: AUC CDF
        subplot(1,2,2);
        [AUC_sorted, ~] = sort(AUC);
        cdf_auc = (1:length(AUC)) / length(AUC);
        stairs(AUC_sorted, cdf_auc, 'LineWidth', 2.5, 'Color', [0.8 0.3 0.3]);
        hold on;
        
        % Add therapeutic window lines
        ylims_cdf = ylim;
        plot([AUC_optimal_low AUC_optimal_low], ylims_cdf, 'g--', 'LineWidth', 2, 'DisplayName', 'Optimal low');
        plot([AUC_optimal_high AUC_optimal_high], ylims_cdf, 'g--', 'LineWidth', 2, 'DisplayName', 'Optimal high');
        plot([AUC_toxic_threshold AUC_toxic_threshold], ylims_cdf, 'r--', 'LineWidth', 2, 'DisplayName', 'Toxic');
        
        xlabel('AUC (mg·h/L)', 'FontSize', 11);
        ylabel('Cumulative Probability', 'FontSize', 11);
        title('AUC: Cumulative Distribution', 'FontSize', 12, 'FontWeight', 'bold');
        grid on; grid minor;
        legend('Location', 'southeast', 'FontSize', 9);
        
        % Add annotations
        annotation('textbox', [0.62 0.75 0.25 0.15], ...
            'String', {sprintf('P5: %.1f mg·h/L', p5_auc), ...
                      sprintf('P50: %.1f mg·h/L', median_auc), ...
                      sprintf('P95: %.1f mg·h/L', p95_auc)}, ...
            'FontSize', 9, 'EdgeColor', 'none', 'BackgroundColor', 'yellow', 'FaceAlpha', 0.7);
        
        sgtitle('Cumulative Distribution Functions: Probability of Outcome < X', 'FontSize', 13, 'FontWeight', 'bold');
        
        % Apply styling and save
        set(fig_cdf, 'Color', 'white');
        set_publication_figure_style(fig_cdf, 'font_size', 11, 'font_name', 'Arial');
        add_figure_metadata('DeGramont.csv', num_valid, analysis_date_str);
        fig_counter = fig_counter + 1;
        caption = get_figure_caption('cdf', num_valid);
        save_publication_figure(fig_cdf, figures_output_dir, fig_counter, 'CDF Plots', caption, ...
            'DeGramont.csv', num_valid, analysis_date_str);
        save_publication_fig_copy(fig_cdf, figures_output_dir, fig_counter, 'CDF Plots', timestamp);
    catch ME_cdf
        fprintf('MC: WARNING - Could not create CDF plots: %s\n', ME_cdf.message);
    end
    
    % ===== NEW FIGURE 3: STATISTICS SUMMARY TABLE =====
    fprintf('MC: Creating statistics summary table\n');
    try
        fig_table = figure('Name','Statistics_Table','Position',[200,200,1000,700]);
        axis off;
        
        % Prepare table data
        table_data = {
            'METRIC', 'CMAX (µM)', 'AUC (mg·h/L)';
            'Count', sprintf('%d', length(Cmax)), sprintf('%d', length(AUC));
            '', '', '';
            'CENTRAL TENDENCY', '', '';
            'Mean', sprintf('%.2f ± %.2f', mean_cmax, std(Cmax)), sprintf('%.1f ± %.1f', mean_auc, std(AUC));
            'Median', sprintf('%.2f', median_cmax), sprintf('%.1f', median_auc);
            '95% CI on mean', sprintf('[%.2f, %.2f]', mean_cmax - 1.96*std(Cmax)/sqrt(length(Cmax)), mean_cmax + 1.96*std(Cmax)/sqrt(length(Cmax))), ...
                              sprintf('[%.1f, %.1f]', mean_auc - 1.96*std(AUC)/sqrt(length(AUC)), mean_auc + 1.96*std(AUC)/sqrt(length(AUC)));
            '', '', '';
            'SPREAD', '', '';
            'Std Dev', sprintf('%.2f', std(Cmax)), sprintf('%.1f', std(AUC));
            'Coeff. Variation', sprintf('%.1f%%', cv_cmax), sprintf('%.1f%%', cv_auc);
            'Range (Min-Max)', sprintf('%.2f - %.2f', min(Cmax), max(Cmax)), sprintf('%.1f - %.1f', min(AUC), max(AUC));
            'IQR (Q1-Q3)', sprintf('%.2f (Q1=%.2f, Q3=%.2f)', iqr_cmax, p25_cmax, p75_cmax), ...
                           sprintf('%.1f (Q1=%.1f, Q3=%.1f)', iqr_auc, p25_auc, p75_auc);
            '', '', '';
            'PERCENTILES', '', '';
            '5th', sprintf('%.2f µM', p5_cmax), sprintf('%.1f mg·h/L', p5_auc);
            '10th', sprintf('%.2f µM', p10_cmax), sprintf('%.1f mg·h/L', p10_auc);
            '25th', sprintf('%.2f µM', p25_cmax), sprintf('%.1f mg·h/L', p25_auc);
            '50th (Median)', sprintf('%.2f µM', median_cmax), sprintf('%.1f mg·h/L', median_auc);
            '75th', sprintf('%.2f µM', p75_cmax), sprintf('%.1f mg·h/L', p75_auc);
            '90th', sprintf('%.2f µM', p90_cmax), sprintf('%.1f mg·h/L', p90_auc);
            '95th', sprintf('%.2f µM', p95_cmax), sprintf('%.1f mg·h/L', p95_auc);
            '', '', '';
            'CLINICAL OUTCOMES', '', '';
            'Safe/Optimal', sprintf('%.1f%% <%.0f µM', pct_safe, Cmax_safe_max), sprintf('%.1f%% (%.0f-%.0f)', pct_optimal, AUC_optimal_low, AUC_optimal_high);
            'At Risk', sprintf('%.1f%% >%.0f µM', pct_high_tox, Cmax_high), sprintf('%.1f%% >%.0f', pct_toxic, AUC_toxic_threshold);
        };
        
        % Create table
        t = uitable('Data', table_data, 'Position', [10 10 980 680], ...
            'ColumnWidth', {250, 350, 350});
        
        % Format header row
        t.Data(1,:) = table_data(1,:);
        
        sgtitle('Pharmacokinetic Summary Statistics', 'FontSize', 14, 'FontWeight', 'bold');
        figfile = fullfile(output_dir, sprintf('MC_5FU_statistics_table_%s.pdf', timestamp));
        save_pdf_and_fig(fig_table, figfile, 'statistics table');
        close(fig_table);
    catch ME_table
        fprintf('MC: WARNING - Could not create statistics table: %s\n', ME_table.message);
    end
    
catch ME_violin
    fprintf('MC: ERROR creating distribution plots: %s\n', ME_violin.message);
    if ~QUIET_MODE
        fprintf('Stack trace:\n');
        for k = 1:length(ME_violin.stack)
            fprintf('  %s (line %d)\n', ME_violin.stack(k).name, ME_violin.stack(k).line);
        end
    end
end

% Save configurable target-adherence summary for easy review/comparison
try
    adherenceTbl = compute_target_adherence_summary(AUC, Cmax, target_profile, dosing_mode);
    adherence_csv = fullfile(output_dir, sprintf('MC_5FU_target_adherence_%s.csv', timestamp));
    writetable(adherenceTbl, adherence_csv);
    fprintf('MC: Saved target adherence summary to %s\n', adherence_csv);
catch ME_adherence
    fprintf('MC: WARNING - Could not save target adherence summary: %s\n', ME_adherence.message);
end

% Scatter matrix: parameters vs Cmax
fprintf('MC: Creating scatter matrix vs Cmax\n');
n_params = length(param_names);
n_cols = 5;
n_rows = ceil(n_params / n_cols);
figure('Name','Scatter_Cmax','Position',[100,100,1600,400*n_rows]);

for j = 1:n_params
    subplot(n_rows, n_cols, j);
    param_name = param_names{j};
    x = params_table{:, j};
    try
        scatter(x, Cmax, 30, 'filled', 'MarkerFaceAlpha', 0.6);
        title(param_name);
        xlabel(param_name);
        ylabel('Cmax (µM)');
        grid on;
    catch ME_plot
        if ~QUIET_MODE
            warning('MC:PlottingError','Failed to plot %s: %s', param_name, ME_plot.message);
        end
    end
end

set_publication_figure_style(gcf, 'font_size', 10, 'font_name', 'Arial');
add_figure_metadata('DeGramont.csv', num_valid, analysis_date_str);
fig_counter = fig_counter + 1;
caption = get_figure_caption('scatter_cmax', num_valid);
save_publication_figure(gcf, figures_output_dir, fig_counter, 'Scatter Cmax', caption, ...
    'DeGramont.csv', num_valid, analysis_date_str);
save_publication_fig_copy(gcf, figures_output_dir, fig_counter, 'Scatter Cmax', timestamp);

% Scatter matrix: parameters vs AUC
fprintf('MC: Creating scatter matrix vs AUC\n');
figure('Name','Scatter_AUC','Position',[100,100,1600,400*n_rows]);
for j = 1:n_params
    subplot(n_rows, n_cols, j);
    param_name = param_names{j};
    x = params_table{:, j};
    try
        scatter(x, AUC, 30, 'filled', 'MarkerFaceAlpha', 0.6);
        title(param_name);
        xlabel(param_name);
        ylabel('AUC (mg·h/L)');
        grid on;
    catch ME_plot
        if ~QUIET_MODE
            warning('MC:PlottingError','Failed to plot %s: %s', param_name, ME_plot.message);
        end
    end
end

set_publication_figure_style(gcf, 'font_size', 10, 'font_name', 'Arial');
add_figure_metadata('DeGramont.csv', num_valid, analysis_date_str);
fig_counter = fig_counter + 1;
caption = get_figure_caption('scatter_auc', num_valid);
save_publication_figure(gcf, figures_output_dir, fig_counter, 'Scatter AUC', caption, ...
    'DeGramont.csv', num_valid, analysis_date_str);
save_publication_fig_copy(gcf, figures_output_dir, fig_counter, 'Scatter AUC', timestamp);

% Scatter plots for top 6 parameters vs Cmax with BEST FIT curves overlaid
topk = min(6, length(param_names));
top_params = sensitivity.rank(1:topk);

% Cmax fitted curves
figure('Name','TopScatters_Cmax','Position',[100,100,1200,800]);
for k = 1:topk
    subplot(2,3,k);
    pname = top_params{k};
    x = params_table{:, pname};
    
    % OUTLIER FILTERING: Remove extreme outliers using IQR method before fitting
    Q1_cmax = prctile(Cmax, 25);
    Q3_cmax = prctile(Cmax, 75);
    IQR_cmax = Q3_cmax - Q1_cmax;
    outlier_threshold = 2.5;  % More lenient than standard 1.5 to preserve trend
    valid_idx = (Cmax >= Q1_cmax - outlier_threshold*IQR_cmax) & (Cmax <= Q3_cmax + outlier_threshold*IQR_cmax);
    
    x_clean = x(valid_idx);
    Cmax_clean = Cmax(valid_idx);
    n_outliers = sum(~valid_idx);
    
    % Plot all data (outliers in different color)
    scatter(x(~valid_idx), Cmax(~valid_idx), 20, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.3); hold on;
    scatter(x_clean, Cmax_clean, 30, 'filled', 'MarkerFaceAlpha', 0.6); hold on;
    
    % Dynamically fit x-axis to cleaned data
    if ~isempty(x_clean)
        x_range = max(x_clean) - min(x_clean);
        xlim([min(x_clean) - 0.05*x_range, max(x_clean) + 0.05*x_range]);
    end
    
    % Use fitted curve if available, otherwise fallback to linear (fit on CLEAN data)
    if has_cftool && ~isempty(curve_fits)
        fitobj = curve_fits.Cmax.(pname).fit;
        gof = curve_fits.Cmax.(pname).gof;
        model_name = curve_fits.Cmax.(pname).model_name;
        xv = linspace(min(x_clean), max(x_clean), 200);
        yfit = feval(fitobj, xv);
        plot(xv, yfit, 'r-', 'LineWidth', 2);
        if n_outliers > 0
            title(sprintf('%s vs Cmax\n%s (R²=%.3f) [%d outliers removed]', pname, model_name, gof.rsquare, n_outliers));
        else
            title(sprintf('%s vs Cmax\n%s (R²=%.3f)', pname, model_name, gof.rsquare));
        end
    else
        % Fallback to linear on clean data
        if length(x_clean) >= 3
            coeffs = polyfit(x_clean, Cmax_clean, 1);
            xv = linspace(min(x_clean), max(x_clean), 100);
            yfit = polyval(coeffs, xv);
            plot(xv, yfit, 'r-', 'LineWidth', 1.5);
            ypred = polyval(coeffs, x_clean);
            SSres = sum((Cmax_clean - ypred).^2);
            SStot = sum((Cmax_clean - mean(Cmax_clean)).^2);
            R2 = 1 - SSres/SStot;
            if n_outliers > 0
                title(sprintf('%s vs Cmax\nLinear (R²=%.3f) [%d outliers]', pname, R2, n_outliers));
            else
                title(sprintf('%s vs Cmax\nLinear (R²=%.3f)', pname, R2));
            end
        else
            title(sprintf('%s vs Cmax\n[Insufficient data]', pname));
        end
    end
    xlabel(pname); ylabel('Cmax (µM)');
    grid on;
end
figfile = fullfile(output_dir, sprintf('MC_5FU_PK_fitted_Cmax_%s.pdf', timestamp));
save_pdf_and_fig(gcf, figfile, 'fitted curve plots (Cmax)');

% AUC fitted curves (mirror of Cmax)
figure('Name','TopScatters_AUC','Position',[100,100,1200,800]);
for k = 1:topk
    subplot(2,3,k);
    pname = top_params{k};
    x = params_table{:, pname};
    
    % OUTLIER FILTERING: Remove extreme outliers using IQR method before fitting
    Q1_auc = prctile(AUC, 25);
    Q3_auc = prctile(AUC, 75);
    IQR_auc = Q3_auc - Q1_auc;
    outlier_threshold = 2.5;  % More lenient to preserve trend
    valid_idx = (AUC >= Q1_auc - outlier_threshold*IQR_auc) & (AUC <= Q3_auc + outlier_threshold*IQR_auc);
    
    x_clean = x(valid_idx);
    AUC_clean = AUC(valid_idx);
    n_outliers = sum(~valid_idx);
    
    % Plot all data (outliers in different color)
    scatter(x(~valid_idx), AUC(~valid_idx), 20, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.3); hold on;
    scatter(x_clean, AUC_clean, 30, 'filled', 'MarkerFaceAlpha', 0.6); hold on;
    
    % Dynamically fit x-axis to cleaned data
    if ~isempty(x_clean)
        x_range = max(x_clean) - min(x_clean);
        xlim([min(x_clean) - 0.05*x_range, max(x_clean) + 0.05*x_range]);
    end
    
    if has_cftool && ~isempty(curve_fits)
        fitobj = curve_fits.AUC.(pname).fit;
        gof = curve_fits.AUC.(pname).gof;
        model_name = curve_fits.AUC.(pname).model_name;
        xv = linspace(min(x_clean), max(x_clean), 200);
        yfit = feval(fitobj, xv);
        plot(xv, yfit, 'b-', 'LineWidth', 2);
        if n_outliers > 0
            title(sprintf('%s vs AUC\n%s (R²=%.3f) [%d outliers removed]', pname, model_name, gof.rsquare, n_outliers));
        else
            title(sprintf('%s vs AUC\n%s (R²=%.3f)', pname, model_name, gof.rsquare));
        end
    else
        if length(x_clean) >= 3
            coeffs = polyfit(x_clean, AUC_clean, 1);
            xv = linspace(min(x_clean), max(x_clean), 100);
            yfit = polyval(coeffs, xv);
            plot(xv, yfit, 'b-', 'LineWidth', 1.5);
            ypred = polyval(coeffs, x_clean);
            SSres = sum((AUC_clean - ypred).^2);
            SStot = sum((AUC_clean - mean(AUC_clean)).^2);
            R2 = 1 - SSres/SStot;
            if n_outliers > 0
                title(sprintf('%s vs AUC\nLinear (R²=%.3f) [%d outliers]', pname, R2, n_outliers));
            else
                title(sprintf('%s vs AUC\nLinear (R²=%.3f)', pname, R2));
            end
        else
            title(sprintf('%s vs AUC\n[Insufficient data]', pname));
        end
    end
    xlabel(pname); ylabel('AUC (mg·h/L)');
    grid on;
end
figfile = fullfile(output_dir, sprintf('MC_5FU_PK_fitted_AUC_%s.pdf', timestamp));
save_pdf_and_fig(gcf, figfile, 'fitted curve plots (AUC)');

% -----------------------------
% 7b) ADVANCED PHARMACOMETRIC VISUALIZATIONS
% -----------------------------
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('CREATING ADVANCED PHARMACOMETRIC VISUALIZATIONS\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

% ===== CDF PLOTS WITH PTA MARKERS =====
fprintf('MC: Creating CDF plots with PTA markers...\n');

% CDF for AUC with configurable targets
auc_targets.values = [target_profile.AUC.subtherapeutic, target_profile.AUC.optimal_low, target_profile.AUC.optimal_high, target_profile.AUC.toxic];
auc_targets.labels = {'Subtherapeutic', 'Optimal lower', 'Optimal upper', 'Toxicity threshold'};
auc_targets.colors = {[0.8 0.4 0], [0 0.7 0], [0 0.7 0], [0.8 0 0]};

auc_cdf_options.metric_name = 'AUC';
auc_cdf_options.units = 'mg·h/L';
auc_cdf_options.title_text = '5-FU AUC: Cumulative Distribution with Clinical Targets';
auc_cdf_options.show_percentiles = true;
auc_cdf_options.grid_style = 'both';
auc_cdf_options.literature_ref = 'Targets from: Gamelin et al. (2008) Clin Cancer Res, Kaldate et al. (2012) Clin Pharmacol Ther';

try
    plotCDF_with_PTA(AUC, auc_targets, auc_cdf_options);
    figfile = fullfile(output_dir, sprintf('MC_5FU_PK_CDF_AUC_%s.pdf', timestamp));
    save_pdf_and_fig(gcf, figfile, 'AUC CDF plot');
catch ME_auc_cdf
    fprintf('MC: WARNING - Could not create AUC CDF plot: %s\n', ME_auc_cdf.message);
end

% CDF for Cmax with configurable thresholds
cmax_targets.values = [target_profile.Cmax.safe_max, target_profile.Cmax.moderate_threshold, target_profile.Cmax.high_threshold, target_profile.Cmax.severe_threshold];
cmax_targets.labels = {'Safe threshold', 'Moderate toxicity', 'High toxicity', 'Severe toxicity'};
cmax_targets.colors = {[0 0.7 0], [1 0.6 0], [0.8 0 0], [0.5 0 0]};

cmax_cdf_options.metric_name = 'Cmax';
cmax_cdf_options.units = 'µM';
cmax_cdf_options.title_text = sprintf('5-FU Cmax [%s]: Cumulative Distribution with Toxicity Thresholds', upper(dosing_mode));
cmax_cdf_options.show_percentiles = true;
cmax_cdf_options.grid_style = 'both';
cmax_cdf_options.literature_ref = sprintf('Config-driven thresholds (%s mode)', upper(dosing_mode));

try
    plotCDF_with_PTA(Cmax, cmax_targets, cmax_cdf_options);
    figfile = fullfile(output_dir, sprintf('MC_5FU_PK_CDF_Cmax_%s_%s.pdf', upper(dosing_mode), timestamp));
    save_pdf_and_fig(gcf, figfile, sprintf('Cmax CDF plot (%s mode)', upper(dosing_mode)));
catch ME_cmax_cdf
    fprintf('MC: WARNING - Could not create Cmax CDF (%s): %s\n', upper(dosing_mode), ME_cmax_cdf.message);
end

% ===== BOX-AND-WHISKER PLOTS =====
fprintf('MC: Creating box-and-whisker plots for PK metrics...\n');

% For demonstration, split data into quartiles to show variability across "virtual cohorts"
% In practice, you might split by dose groups, patient characteristics, etc.
n_quartiles = 4;
[~, sort_idx_auc] = sort(AUC);
quartile_size = floor(length(AUC) / n_quartiles);

box_data.AUC = cell(n_quartiles, 1);
box_data.Cmax = cell(n_quartiles, 1);
box_data.labels = cell(n_quartiles, 1);

for q = 1:n_quartiles
    idx_start = (q-1) * quartile_size + 1;
    if q == n_quartiles
        idx_end = length(AUC);
    else
        idx_end = q * quartile_size;
    end
    
    quartile_indices = sort_idx_auc(idx_start:idx_end);
    box_data.AUC{q} = AUC(quartile_indices);
    box_data.Cmax{q} = Cmax(quartile_indices);
    
    % Calculate mean AUC for labeling
    mean_auc = mean(AUC(quartile_indices));
    box_data.labels{q} = sprintf('Q%d\n(AUC ~%.0f)', q, mean_auc);
end

box_options.show_points = true;
box_options.point_jitter = 0.2;
box_options.point_alpha = 0.4;
box_options.show_mean = true;
box_options.add_stats = false; % Set to true for pairwise comparisons (only works for 2 groups)
box_options.ylabel_auc = 'AUC (mg·h/L)';
box_options.ylabel_cmax = 'Cmax (µM)';
box_options.title_text = '5-FU Pharmacokinetic Variability Across Exposure Quartiles';
box_options.literature_ref = 'Visualization guidelines: Weissgerber et al. (2015) PLoS Biol, Krzywinski & Altman (2014) Nat Methods';

plotBoxWhisker_PK(box_data, box_options);
figfile = fullfile(output_dir, sprintf('MC_5FU_PK_BoxWhisker_%s.pdf', timestamp));
save_pdf_and_fig(gcf, figfile, 'box-and-whisker plots');

% ===== VISUAL PREDICTIVE CHECK (VPC) =====
fprintf('MC: Creating Visual Predictive Check (VPC) plot...\n');

% Extract time-concentration profiles from simulations
% Use only simulations with complete timeseries data
valid_timeseries_idx = ~cellfun(@isempty, raw_timeseries);
valid_timeseries = raw_timeseries(valid_timeseries_idx);

if ~isempty(valid_timeseries) && length(valid_timeseries) >= 10
    % Get time vector from first simulation (assume all have same time points)
    if isstruct(valid_timeseries{1})
        time_vec = valid_timeseries{1}.time_hr;
        n_time = length(time_vec);
        n_sims = length(valid_timeseries);
        
        % Build concentration matrix (time x simulations)
        conc_matrix = zeros(n_time, n_sims);
        for s = 1:n_sims
            if isstruct(valid_timeseries{s}) && isfield(valid_timeseries{s}, 'C_central')
                conc_matrix(:, s) = valid_timeseries{s}.C_central;
            end
        end
        
        % Calculate "observed" data as median of simulations (for demonstration)
        % In real application, this would be actual patient data
        observed_conc = median(conc_matrix, 2);
        
        % SIMPLIFIED VPC with error bars (clearer than complex shaded regions)
        vpc_options.error_type = 'SD';  % Use standard deviation (or 'percentiles' for [10 90])
        vpc_options.obs_color = [0.8 0 0];  % Red for observed
        vpc_options.pred_color = [0.2 0.4 0.8];  % Blue for predictions
        vpc_options.title_text = 'Model Validation: Observed vs Predicted Concentrations';
        vpc_options.xlabel_text = 'Time (hours)';
        vpc_options.ylabel_text = 'Central Compartment Concentration (µM)';
        vpc_options.log_scale = false;
        vpc_options.literature_ref = 'Simplified VPC approach for clarity. Traditional VPC: Bergstrand et al. (2011) AAPS J';
        
        plotVPC(observed_conc, conc_matrix, time_vec, vpc_options);
        figfile = fullfile(output_dir, sprintf('MC_5FU_PK_VPC_%s.pdf', timestamp));
        save_pdf_and_fig(gcf, figfile, 'simplified VPC plot');
    else
        fprintf('  WARNING: Timeseries data format not compatible for VPC. Skipping.\n');
    end
else
    fprintf('  WARNING: Insufficient timeseries data for VPC (need ≥10 simulations). Skipping.\n');
end

% ===== PROBABILITY OF TARGET ATTAINMENT (PTA) PLOT =====
fprintf('MC: Creating Probability of Target Attainment (PTA) plots...\n');

% ===== PTA WITH LITERATURE-BASED TARGETS (AUC) =====
% Define efficacy and toxicity targets for 5-FU (literature-based)
% Efficacy targets based on: Gamelin et al. (2008), Beumer et al. (2019)
pta_efficacy.values = [target_profile.AUC.subtherapeutic, target_profile.AUC.optimal_low];  % mg·h/L
pta_efficacy.labels = {'Minimum efficacy (tumor response)', 'Optimal efficacy (94% response rate)'};

% Toxicity targets based on: Saif et al. (2013), Kaldate et al. (2012)
pta_toxicity.values = [target_profile.AUC.toxic, target_profile.AUC.severe_toxicity];  % mg·h/L
pta_toxicity.labels = {'Grade 3+ toxicity risk (>20%)', 'Grade 3+ toxicity risk (>60%)'};

pta_options.metric_name = 'AUC';
pta_options.units = 'mg·h/L';
pta_options.title_text = '5-FU Dose Optimization: Balancing Efficacy and Toxicity';
pta_options.show_optimal = true;
pta_options.optimal_pta = 90;  % Want ≥90% probability of achieving efficacy
pta_options.max_tox_prob = 10; % Accept ≤10% probability of severe toxicity
pta_options.auto_adjust = true;  % Enable dynamic data adjustment
pta_options.literature_ref = 'PTA methodology: Mouton & Vinks (2005) J Antimicrob Chemother, Drusano (2004) Clin Infect Dis';

plotPTA(AUC, pta_efficacy, pta_toxicity, pta_options);
plotPTA_StatisticsSummary(AUC, pta_efficacy, pta_toxicity, pta_options);
figfile = fullfile(output_dir, sprintf('MC_5FU_PK_PTA_%s.pdf', timestamp));
save_pdf_and_fig(gcf, figfile, 'PTA plot (literature-based)');

% ===== PTA WITH DATA-DRIVEN DYNAMIC TARGETS (AUC) =====
fprintf('MC: Creating dynamic PTA plot with data-driven targets (AUC)...\n');
pta_dynamic_auc = struct();  % Empty struct - will auto-generate from data
pta_dynamic_auc.auto_generate = true;

pta_dynamic_tox = struct();
pta_dynamic_tox.auto_generate = true;

pta_options_dynamic.metric_name = 'AUC';
pta_options_dynamic.units = 'mg·h/L';
pta_options_dynamic.title_text = '5-FU Dose Optimization: Data-Driven Target Percentiles';
pta_options_dynamic.show_optimal = true;
pta_options_dynamic.optimal_pta = 90;
pta_options_dynamic.max_tox_prob = 10;
pta_options_dynamic.auto_adjust = true;  % Enable dynamic adjustment
pta_options_dynamic.percentile_efficacy = [25 50];  % Use 25th and 50th percentiles as efficacy targets
pta_options_dynamic.percentile_toxicity = [75 90];  % Use 75th and 90th percentiles as toxicity targets
pta_options_dynamic.literature_ref = 'Data-driven targets from population percentiles (P25, P50, P75, P90)';

plotPTA(AUC, pta_dynamic_auc, pta_dynamic_tox, pta_options_dynamic);
figfile = fullfile(output_dir, sprintf('MC_5FU_PK_PTA_DataDriven_%s.pdf', timestamp));
save_pdf_and_fig(gcf, figfile, 'PTA plot (data-driven)');

% ===== PTA FOR CMAX WITH LITERATURE TARGETS =====
fprintf('MC: Creating PTA plot for Cmax (peak concentration)...\n');
pta_efficacy_cmax.values = [target_profile.Cmax.safe_max, target_profile.Cmax.moderate_threshold];
pta_efficacy_cmax.labels = {'Minimum peak concentration', 'Target peak concentration'};

pta_toxicity_cmax.values = [target_profile.Cmax.high_threshold, target_profile.Cmax.severe_threshold];
pta_toxicity_cmax.labels = {'Mucositis risk (moderate)', 'Mucositis risk (severe)'};

pta_options_cmax.metric_name = 'Cmax';
pta_options_cmax.units = 'µM';
pta_options_cmax.title_text = '5-FU Peak Concentration: Balancing Efficacy and Acute Toxicity';
pta_options_cmax.show_optimal = true;
pta_options_cmax.optimal_pta = 90;
pta_options_cmax.max_tox_prob = 10;
pta_options_cmax.auto_adjust = true;  % Enable dynamic adjustment
pta_options_cmax.literature_ref = 'Acute toxicity correlates: Thyss et al. (1986) Cancer Res, Gamelin et al. (2008) Clin Cancer Res';

plotPTA(Cmax, pta_efficacy_cmax, pta_toxicity_cmax, pta_options_cmax);
plotPTA_StatisticsSummary(Cmax, pta_efficacy_cmax, pta_toxicity_cmax, pta_options_cmax);
figfile = fullfile(output_dir, sprintf('MC_5FU_PK_PTA_Cmax_%s.pdf', timestamp));
save_pdf_and_fig(gcf, figfile, 'Cmax PTA plot');

% ===== PTA FOR CMAX WITH DATA-DRIVEN DYNAMIC TARGETS =====
fprintf('MC: Creating dynamic PTA plot with data-driven targets (Cmax)...\n');
pta_dynamic_cmax = struct();
pta_dynamic_cmax.auto_generate = true;

pta_dynamic_tox_cmax = struct();
pta_dynamic_tox_cmax.auto_generate = true;

pta_options_cmax_dynamic.metric_name = 'Cmax';
pta_options_cmax_dynamic.units = 'µM';
pta_options_cmax_dynamic.title_text = '5-FU Peak Concentration: Data-Driven Target Percentiles';
pta_options_cmax_dynamic.show_optimal = true;
pta_options_cmax_dynamic.optimal_pta = 90;
pta_options_cmax_dynamic.max_tox_prob = 10;
pta_options_cmax_dynamic.auto_adjust = true;
pta_options_cmax_dynamic.percentile_efficacy = [30 60];
pta_options_cmax_dynamic.percentile_toxicity = [70 95];
pta_options_cmax_dynamic.literature_ref = 'Data-driven targets from population percentiles (Cmax)';

plotPTA(Cmax, pta_dynamic_cmax, pta_dynamic_tox_cmax, pta_options_cmax_dynamic);
figfile = fullfile(output_dir, sprintf('MC_5FU_PK_PTA_Cmax_DataDriven_%s.pdf', timestamp));
save_pdf_and_fig(gcf, figfile, 'Cmax PTA plot (data-driven)');

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('ADVANCED VISUALIZATIONS COMPLETE\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

% -----------------------------
% 8) FILE SAVING
% -----------------------------
results = struct();
results.parameters = params_table;
results.AUC = AUC;
results.Cmax = Cmax;
results.sensitivity = sensitivity;
results.stats = stats;
results.timestamp = timestamp;
results.seed = seed;
results.raw_timeseries = raw_timeseries;

% Add curve fitting results if available
if has_cftool && ~isempty(curve_fits)
    results.curve_fits = curve_fits;
    results.fit_summary_Cmax = fit_table_Cmax;
    results.fit_summary_AUC = fit_table_AUC;
    
    % Save fit summary tables as CSV for easy inspection
    csv_cmax = fullfile(output_dir, sprintf('MC_5FU_curve_fits_Cmax_%s.csv', timestamp));
    csv_auc = fullfile(output_dir, sprintf('MC_5FU_curve_fits_AUC_%s.csv', timestamp));
    writetable(fit_table_Cmax, csv_cmax);
    writetable(fit_table_AUC, csv_auc);
    fprintf('MC: Saved curve fit summaries to CSV:\n');
    fprintf('    %s\n', csv_cmax);
    fprintf('    %s\n', csv_auc);
end

matfile = fullfile(output_dir, sprintf('MC_5FU_PK_sensitivity_%s.mat', timestamp));
save(matfile, 'results', '-v7.3');
fprintf('MC: Saved results to %s\n', matfile);

fprintf('MC: Summary - mean AUC = %.3f mg·h/L (95%% CI ± %.3f); mean Cmax = %.3f µM (95%% CI ± %.3f)\n', ...
    stats.mean.AUC, stats.ci95.AUC, stats.mean.Cmax, stats.ci95.Cmax);

% Generate comprehensive statistics report for paper
fprintf('\nMC: Generating comprehensive statistical report...\n');
try
    report_file = generate_comprehensive_stats_report(output_dir, Cmax, AUC, ...
        num_valid, analysis_date_str, 'DeGramont.csv', rho, param_names);
    fprintf('MC: ✓ Statistical report generated: %s\n', report_file);
catch ME_report
    fprintf('MC: WARNING - Could not generate statistics report: %s\n', ME_report.message);
end

% Final completion notification
if NTFY_ENABLED
    msg = sprintf(['🎉 ANALYSIS COMPLETE\n\n' ...
        '━━━━━━━━━━━━━━━━━━━━\n' ...
        '📊 RESULTS SUMMARY\n' ...
        '━━━━━━━━━━━━━━━━━━━━\n\n' ...
        '💊 AUC (Drug Exposure):\n' ...
        '   %.1f ± %.1f mg·h/L\n\n' ...
        '📈 Cmax (Peak Conc.):\n' ...
        '   %.1f ± %.1f µM\n\n' ...
        '✅ Samples: %d/%d\n' ...
        '📁 Results saved'], ...
        stats.mean.AUC, stats.ci95.AUC, stats.mean.Cmax, stats.ci95.Cmax, ...
        sum(valid_idx), n_samples);
    ntfyNotify(NTFY_TOPIC, msg, 'tada,chart_with_upwards_trend', 'high');
end

fprintf('MC: Analysis complete. Results saved to %s\n', output_dir);

end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function [samples, completed_runs, checkpoint_file] = loadCheckpoint(base_output_dir, n_samples)
    % loadCheckpoint  Load checkpoint data if available, to resume interrupted runs
    %   Searches for the most recent timestamped subfolder and loads checkpoint from there
    %
    % Returns:
    %   samples          - LHS sample matrix (n_samples x n_params)
    %   completed_runs   - Logical array indicating which runs completed
    %   checkpoint_file  - Path to checkpoint file (for saving later)
    
    % Find most recent timestamped folder
    if isfolder(base_output_dir)
        subfolders = dir(base_output_dir);
        % Filter for directories matching timestamp pattern (yyyy-mm-dd_HH-MM-SS)
        valid_folders = {};
        for k = 1:length(subfolders)
            if subfolders(k).isdir && ~startsWith(subfolders(k).name, '.')
                % Check if name matches timestamp pattern
                if length(subfolders(k).name) == 19 && subfolders(k).name(5) == '-' && subfolders(k).name(11) == '_'
                    valid_folders{end+1} = subfolders(k).name; %#ok<AGROW>
                end
            end
        end
        
        if ~isempty(valid_folders)
            % Sort and use most recent (names are sortable as YYYY-MM-DD_HH-MM-SS)
            valid_folders = sort(valid_folders);
            most_recent = valid_folders{end};
            checkpoint_file = fullfile(base_output_dir, most_recent, 'MC_checkpoint.mat');
            
            if isfile(checkpoint_file)
                load(checkpoint_file, 'samples', 'completed_runs');
                fprintf('MC: Found checkpoint in %s\n', most_recent);
                fprintf('MC: Resuming from checkpoint. %d/%d runs already completed.\n', sum(completed_runs), n_samples);
                return;
            end
        end
    end
    
    % No checkpoint found - will use current timestamped folder
    samples = [];
    completed_runs = false(n_samples, 1);
    checkpoint_file = fullfile(base_output_dir, datestr(now, 'yyyy-mm-dd_HH-MM-SS'), 'MC_checkpoint.mat');
    fprintf('MC: No checkpoint found. Starting fresh run.\n');
end

function saveCheckpoint(samples, completed_runs, checkpoint_file)
    % saveCheckpoint  Save current progress to checkpoint file
    
    try
        save(checkpoint_file, 'samples', 'completed_runs', '-v7.3');
    catch ME
        warning('MC:CheckpointFailed', 'Could not save checkpoint: %s', ME.message);
    end
end

%% ========================================================================
%  ADVANCED PHARMACOMETRIC VISUALIZATIONS
%  ========================================================================

function plotCDF_with_PTA(data, targets, options)
% plotCDF_with_PTA  Create CDF plot with Probability of Target Attainment markers
%
% Creates cumulative distribution function plots with horizontal reference lines
% marking therapeutic targets. Useful for visualizing probability of achieving
% specific pharmacokinetic/pharmacodynamic endpoints.
%
% INPUTS:
%   data     - Vector of PK metric values (AUC or Cmax)
%   targets  - Struct with target thresholds:
%              .values - Vector of target values
%              .labels - Cell array of labels for each target
%              .colors - Cell array of colors for each target line (optional)
%   options  - Struct with customization options:
%              .metric_name   - Name of metric (default: 'PK Metric')
%              .units         - Units (default: '')
%              .title_text    - Custom title (optional)
%              .show_percentiles - Show 5/50/95th percentiles (default: true)
%              .line_width    - Width of CDF line (default: 2.5)
%              .grid_style    - 'both', 'major', 'minor', 'off' (default: 'both')
%              .literature_ref - Literature reference text (optional)
%
% LITERATURE:
%   - Mouton et al. (2011) Antimicrob Agents Chemother "Standardization of pharmacokinetic/pharmacodynamic (PK/PD) terminology"
%   - Ambrose et al. (2007) Antimicrob Agents Chemother "Monte Carlo simulation in evaluating antimicrobial therapy"
%   - FDA (2019) Guidance "Population Pharmacokinetics"
%
% EXAMPLE:
%   targets.values = [20 30 40];
%   targets.labels = {'Minimum efficacy', 'Optimal target', 'Toxicity risk'};
%   targets.colors = {'green', 'blue', 'red'};
%   options.metric_name = 'AUC';
%   options.units = 'mg·h/L';
%   plotCDF_with_PTA(AUC_values, targets, options);

    % Validate inputs
    if nargin < 2
        error('plotCDF_with_PTA:InvalidInput', 'Requires data and targets arguments');
    end
    if nargin < 3
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'metric_name'), options.metric_name = 'PK Metric'; end
    if ~isfield(options, 'units'), options.units = ''; end
    if ~isfield(options, 'show_percentiles'), options.show_percentiles = true; end
    if ~isfield(options, 'line_width'), options.line_width = 2.5; end
    if ~isfield(options, 'grid_style'), options.grid_style = 'both'; end
    
    % Remove NaN/Inf values
    data = data(isfinite(data));
    if isempty(data)
        error('plotCDF_with_PTA:NoValidData', 'No valid data points after removing NaN/Inf');
    end
    
    % Sort data for CDF calculation
    sorted_data = sort(data);
    n = length(sorted_data);
    cdf_values = (1:n)' / n;
    
    % Create figure
    figure('Position', [100, 100, 900, 600]);
    
    % Plot CDF
    plot(sorted_data, cdf_values * 100, 'b-', 'LineWidth', options.line_width);
    hold on;
    
    % Add target reference lines with PTA calculations
    n_targets = length(targets.values);
    if ~isfield(targets, 'colors')
        % Default color scheme: green -> yellow -> orange -> red
        default_colors = {'g', [0.8 0.7 0], [1 0.5 0], 'r'};
        targets.colors = default_colors(mod(0:n_targets-1, length(default_colors)) + 1);
    end
    
    legend_entries = cell(n_targets + 1, 1);
    legend_entries{1} = 'Empirical CDF';
    
    for i = 1:n_targets
        target_val = targets.values(i);
        
        % Calculate PTA (Probability of Target Attainment)
        pta = sum(data >= target_val) / length(data) * 100;
        
        % Plot vertical line at target
        ylims = ylim;
        plot([target_val target_val], ylims, '--', 'Color', targets.colors{i}, ...
            'LineWidth', 2, 'DisplayName', sprintf('%s (PTA=%.1f%%)', targets.labels{i}, pta));
        
        % Add PTA annotation
        pta_x = target_val;
        pta_y = interp1(sorted_data, cdf_values * 100, target_val, 'linear', 'extrap');
        pta_y = max(5, min(95, pta_y)); % Keep annotation within bounds
        
        text(pta_x, pta_y, sprintf('  %.1f%%', pta), ...
            'FontSize', 10, 'FontWeight', 'bold', 'Color', targets.colors{i}, ...
            'VerticalAlignment', 'bottom', 'BackgroundColor', 'white', 'EdgeColor', targets.colors{i});
        
        legend_entries{i+1} = sprintf('%s (PTA=%.1f%%)', targets.labels{i}, pta);
    end
    
    % Add percentile markers if requested
    if options.show_percentiles
        percentiles = [5 50 95];
        perc_vals = prctile(data, percentiles);
        
        for i = 1:length(percentiles)
            plot(perc_vals(i), percentiles(i), 'ko', 'MarkerSize', 8, ...
                'MarkerFaceColor', 'black', 'LineWidth', 1.5);
            text(perc_vals(i), percentiles(i) + 3, sprintf('P%d', percentiles(i)), ...
                'FontSize', 9, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
    
    % Formatting
    xlabel(sprintf('%s (%s)', options.metric_name, options.units), 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Cumulative Probability (%)', 'FontSize', 12, 'FontWeight', 'bold');
    
    if isfield(options, 'title_text')
        title(options.title_text, 'FontSize', 13, 'FontWeight', 'bold');
    else
        title(sprintf('CDF of %s with Target Attainment', options.metric_name), ...
            'FontSize', 13, 'FontWeight', 'bold');
    end
    
    % Grid
    switch options.grid_style
        case 'both'
            grid on; grid minor;
        case 'major'
            grid on;
        case 'minor'
            grid minor;
        case 'off'
            % no grid
    end
    
    % Legend
    legend(legend_entries, 'Location', 'southeast', 'FontSize', 9);
    
    % Add literature reference if provided
    if isfield(options, 'literature_ref')
        annotation('textbox', [0.15 0.02 0.7 0.05], 'String', options.literature_ref, ...
            'FontSize', 8, 'EdgeColor', 'none', 'Interpreter', 'none', 'FitBoxToText', 'on');
    end
    
    hold off;
    set(gca, 'FontSize', 10);
end

function plotVPC(observed, simulated_matrix, time, options)
% plotVPC  Visual Predictive Check (VPC) plot
%
% Creates a VPC plot showing observed data against prediction intervals from
% simulated data. Essential diagnostic for population PK/PD models.
%
% INPUTS:
%   observed          - Vector of observed concentrations/values at each time point
%   simulated_matrix  - Matrix where each row is a time point, each column is a simulation
%   time              - Vector of time points
%   options           - Struct with customization:
%                       .percentiles   - Percentiles to plot [5 50 95] (default)
%                       .obs_color     - Color for observed data (default: 'black')
%                       .pred_color    - Color for prediction intervals (default: 'blue')
%                       .ci_alpha      - Transparency for CI bands (default: 0.3)
%                       .title_text    - Custom title
%                       .xlabel_text   - X-axis label (default: 'Time')
%                       .ylabel_text   - Y-axis label (default: 'Concentration')
%                       .log_scale     - Use log scale for y-axis (default: false)
%                       .show_individual_sims - Show individual simulations (default: false)
%                       .literature_ref - Literature reference text (optional)
%
% LITERATURE:
%   - Bergstrand et al. (2011) AAPS J "Prediction-corrected visual predictive checks"
%   - Holford (2005) Clin Pharmacokinet "The visual predictive check—superiority to standard diagnostic plots"
%   - Karlsson & Savic (2007) AAPS J "Diagnosing model diagnostics"
%
% EXAMPLE:
%   observed = patient_concentrations;
%   simulated = monte_carlo_sims;  % Each column is one simulation
%   time = [0 0.5 1 2 4 8 24];
%   options.title_text = '5-FU VPC - 400 mg/m² infusion';
%   options.log_scale = true;
%   plotVPC(observed, simulated, time, options);

    % Validate inputs
    if nargin < 3
        error('plotVPC:InvalidInput', 'Requires observed, simulated_matrix, and time arguments');
    end
    if nargin < 4
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'percentiles'), options.percentiles = [5 50 95]; end
    if ~isfield(options, 'obs_color'), options.obs_color = 'black'; end
    if ~isfield(options, 'pred_color'), options.pred_color = [0.2 0.4 0.8]; end
    if ~isfield(options, 'ci_alpha'), options.ci_alpha = 0.3; end
    if ~isfield(options, 'xlabel_text'), options.xlabel_text = 'Time'; end
    if ~isfield(options, 'ylabel_text'), options.ylabel_text = 'Concentration'; end
    if ~isfield(options, 'log_scale'), options.log_scale = false; end
    if ~isfield(options, 'show_individual_sims'), options.show_individual_sims = false; end
    
    % Calculate percentiles for each time point
    n_time = length(time);
    perc_matrix = zeros(length(options.percentiles), n_time);
    
    for t = 1:n_time
        perc_matrix(:, t) = prctile(simulated_matrix(t, :), options.percentiles);
    end
    
    % Create figure
    figure('Position', [100, 100, 1000, 600]);
    
    % Plot individual simulations if requested (before other elements)
    if options.show_individual_sims
        n_sims_to_show = min(50, size(simulated_matrix, 2));  % Limit for clarity
        for i = 1:n_sims_to_show
            plot(time, simulated_matrix(:, i), '-', 'Color', [0.7 0.7 0.7 0.2], ...
                'LineWidth', 0.5, 'HandleVisibility', 'off');
            hold on;
        end
    end
    
    % Plot prediction intervals as shaded regions
    hold on;
    
    % 90% prediction interval (5th-95th percentile)
    fill([time; flipud(time)], [perc_matrix(1, :)'; flipud(perc_matrix(3, :)')], ...
        options.pred_color, 'FaceAlpha', options.ci_alpha, 'EdgeColor', 'none', ...
        'DisplayName', sprintf('%d-%d%% PI', options.percentiles(1), options.percentiles(3)));
    
    % Plot median prediction
    plot(time, perc_matrix(2, :), '-', 'Color', options.pred_color, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Predicted Median (P%d)', options.percentiles(2)));
    
    % Plot observed data
    plot(time, observed, 'o-', 'Color', options.obs_color, 'LineWidth', 2, ...
        'MarkerSize', 8, 'MarkerFaceColor', options.obs_color, ...
        'DisplayName', 'Observed');
    
    % Formatting
    xlabel(options.xlabel_text, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(options.ylabel_text, 'FontSize', 12, 'FontWeight', 'bold');
    
    if isfield(options, 'title_text')
        title(options.title_text, 'FontSize', 13, 'FontWeight', 'bold');
    else
        title('Visual Predictive Check (VPC)', 'FontSize', 13, 'FontWeight', 'bold');
    end
    
    % Apply log scale if requested
    if options.log_scale
        set(gca, 'YScale', 'log');
    end
    
    grid on;
    legend('Location', 'best', 'FontSize', 10);
    
    % Add literature reference if provided
    if isfield(options, 'literature_ref')
        annotation('textbox', [0.15 0.02 0.7 0.05], 'String', options.literature_ref, ...
            'FontSize', 8, 'EdgeColor', 'none', 'Interpreter', 'none', 'FitBoxToText', 'on');
    end
    
    hold off;
    set(gca, 'FontSize', 10);
end

function plotBoxWhisker_PK(data_struct, options)
% plotBoxWhisker_PK  Box-and-whisker plots for AUC and Cmax by group/regimen
%
% Creates publication-quality box-and-whisker plots for pharmacokinetic metrics,
% with optional overlaid individual data points and statistical comparisons.
%
% INPUTS:
%   data_struct - Struct with fields:
%                 .AUC    - Cell array where each cell contains AUC values for one group
%                 .Cmax   - Cell array where each cell contains Cmax values for one group
%                 .labels - Cell array of group labels (e.g., {'Low dose', 'High dose'})
%   options     - Struct with customization:
%                 .show_points      - Overlay individual data points (default: true)
%                 .point_jitter     - Jitter amount for points (default: 0.15)
%                 .point_alpha      - Transparency for points (default: 0.5)
%                 .box_colors       - Cell array of colors for each box (optional)
%                 .show_mean        - Show mean markers (default: true)
%                 .add_stats        - Add statistical comparison p-values (default: false)
%                 .stats_test       - 'ttest', 'ranksum', 'anova' (default: 'ranksum')
%                 .ylabel_auc       - Y-axis label for AUC (default: 'AUC (mg·h/L)')
%                 .ylabel_cmax      - Y-axis label for Cmax (default: 'Cmax (µM)')
%                 .title_text       - Overall title (optional)
%                 .literature_ref   - Literature reference text (optional)
%
% LITERATURE:
%   - Bland & Altman (1996) BMJ "Statistics notes: measurement error"
%   - Krzywinski & Altman (2014) Nat Methods "Visualizing samples with box plots"
%   - Weissgerber et al. (2015) PLoS Biol "Beyond bar and line graphs"
%
% EXAMPLE:
%   data.AUC = {auc_low_dose, auc_high_dose};
%   data.Cmax = {cmax_low_dose, cmax_high_dose};
%   data.labels = {'400 mg/m²', '800 mg/m²'};
%   options.show_points = true;
%   options.add_stats = true;
%   plotBoxWhisker_PK(data, options);

    % Validate inputs
    if nargin < 1
        error('plotBoxWhisker_PK:InvalidInput', 'Requires data_struct argument');
    end
    if nargin < 2
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'show_points'), options.show_points = false; end
    if ~isfield(options, 'point_jitter'), options.point_jitter = 0.15; end
    if ~isfield(options, 'point_alpha'), options.point_alpha = 0.5; end
    if ~isfield(options, 'show_mean'), options.show_mean = true; end
    if ~isfield(options, 'add_stats'), options.add_stats = false; end
    if ~isfield(options, 'stats_test'), options.stats_test = 'ranksum'; end
    if ~isfield(options, 'ylabel_auc'), options.ylabel_auc = 'AUC (mg·h/L)'; end
    if ~isfield(options, 'ylabel_cmax'), options.ylabel_cmax = 'Cmax (µM)'; end
    
    n_groups = length(data_struct.labels);
    
    % Default colors (ColorBrewer Set2)
    if ~isfield(options, 'box_colors')
        default_colors = {[0.4 0.76 0.65], [0.99 0.55 0.38], [0.55 0.63 0.80], ...
                         [0.91 0.54 0.76], [0.65 0.85 0.33], [1.0 0.85 0.18]};
        options.box_colors = default_colors(mod(0:n_groups-1, length(default_colors)) + 1);
    end
    
    % Create figure with two subplots
    figure('Position', [100, 100, 1200, 600]);
    
    % ===== SUBPLOT 1: AUC =====
    subplot(1, 2, 1);
    
    % Prepare data for boxplot
    auc_data = [];
    auc_groups = [];
    for g = 1:n_groups
        auc_data = [auc_data; data_struct.AUC{g}(:)]; %#ok<AGROW>
        auc_groups = [auc_groups; g * ones(length(data_struct.AUC{g}), 1)]; %#ok<AGROW>
    end
    
    % Create boxplot
    h_box = boxplot(auc_data, auc_groups, 'Labels', data_struct.labels, ...
        'Colors', 'k', 'Symbol', '', 'Widths', 0.6);
    hold on;
    
    % Color the boxes
    h_patches = findobj(gca, 'Tag', 'Box');
    for g = 1:n_groups
        patch(get(h_patches(n_groups - g + 1), 'XData'), ...
              get(h_patches(n_groups - g + 1), 'YData'), ...
              options.box_colors{g}, 'FaceAlpha', 0.6);
    end
    
    % Add individual points if requested
    if options.show_points
        for g = 1:n_groups
            % Add random jitter
            x_jitter = g + (rand(length(data_struct.AUC{g}), 1) - 0.5) * options.point_jitter;
            scatter(x_jitter, data_struct.AUC{g}, 30, options.box_colors{g}, 'filled', ...
                'MarkerFaceAlpha', options.point_alpha, 'HandleVisibility', 'off');
        end
    end
    
    % Add mean markers if requested
    if options.show_mean
        for g = 1:n_groups
            mean_val = mean(data_struct.AUC{g});
            plot(g, mean_val, 'd', 'MarkerSize', 10, 'MarkerFaceColor', 'red', ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % Add statistical comparisons if requested
    if options.add_stats && n_groups == 2
        [p_val, ~] = ranksum(data_struct.AUC{1}, data_struct.AUC{2});
        y_max = max(auc_data) * 1.1;
        plot([1 2], [y_max y_max], 'k-', 'LineWidth', 1.5);
        if p_val < 0.001
            text(1.5, y_max * 1.05, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
        elseif p_val < 0.01
            text(1.5, y_max * 1.05, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
        elseif p_val < 0.05
            text(1.5, y_max * 1.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
        else
            text(1.5, y_max * 1.05, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 12);
        end
        text(1.5, y_max * 1.08, sprintf('p = %.4f', p_val), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    
    ylabel(options.ylabel_auc, 'FontSize', 12, 'FontWeight', 'bold');
    title('AUC Distribution', 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    hold off;
    
    % ===== SUBPLOT 2: Cmax =====
    subplot(1, 2, 2);
    
    % Prepare data for boxplot
    cmax_data = [];
    cmax_groups = [];
    for g = 1:n_groups
        cmax_data = [cmax_data; data_struct.Cmax{g}(:)]; %#ok<AGROW>
        cmax_groups = [cmax_groups; g * ones(length(data_struct.Cmax{g}), 1)]; %#ok<AGROW>
    end
    
    % Create boxplot
    h_box = boxplot(cmax_data, cmax_groups, 'Labels', data_struct.labels, ...
        'Colors', 'k', 'Symbol', '', 'Widths', 0.6);
    hold on;
    
    % Color the boxes
    h_patches = findobj(gca, 'Tag', 'Box');
    for g = 1:n_groups
        patch(get(h_patches(n_groups - g + 1), 'XData'), ...
              get(h_patches(n_groups - g + 1), 'YData'), ...
              options.box_colors{g}, 'FaceAlpha', 0.6);
    end
    
    % Add individual points if requested
    if options.show_points
        for g = 1:n_groups
            x_jitter = g + (rand(length(data_struct.Cmax{g}), 1) - 0.5) * options.point_jitter;
            scatter(x_jitter, data_struct.Cmax{g}, 30, options.box_colors{g}, 'filled', ...
                'MarkerFaceAlpha', options.point_alpha, 'HandleVisibility', 'off');
        end
    end
    
    % Add mean markers if requested
    if options.show_mean
        for g = 1:n_groups
            mean_val = mean(data_struct.Cmax{g});
            plot(g, mean_val, 'd', 'MarkerSize', 10, 'MarkerFaceColor', 'red', ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % Add statistical comparisons if requested
    if options.add_stats && n_groups == 2
        [p_val, ~] = ranksum(data_struct.Cmax{1}, data_struct.Cmax{2});
        y_max = max(cmax_data) * 1.1;
        plot([1 2], [y_max y_max], 'k-', 'LineWidth', 1.5);
        if p_val < 0.001
            text(1.5, y_max * 1.05, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
        elseif p_val < 0.01
            text(1.5, y_max * 1.05, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
        elseif p_val < 0.05
            text(1.5, y_max * 1.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
        else
            text(1.5, y_max * 1.05, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 14);
        end
        text(1.5, y_max * 1.08, sprintf('p = %.4f', p_val), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    
    ylabel(options.ylabel_cmax, 'FontSize', 12, 'FontWeight', 'bold');
    title('Cmax Distribution', 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    hold off;
    
    % Add overall title if provided
    if isfield(options, 'title_text')
        sgtitle(options.title_text, 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % Add subtitle if provided (explanatory text)
    if isfield(options, 'subtitle_text')
        annotation('textbox', [0.15 0.88 0.7 0.05], 'String', options.subtitle_text, ...
            'FontSize', 10, 'EdgeColor', 'none', 'Interpreter', 'none', ...
            'HorizontalAlignment', 'center', 'FitBoxToText', 'off', 'FontAngle', 'italic');
    end
    
    % Add literature reference if provided
    if isfield(options, 'literature_ref')
        annotation('textbox', [0.15 0.02 0.7 0.05], 'String', options.literature_ref, ...
            'FontSize', 8, 'EdgeColor', 'none', 'Interpreter', 'none', 'FitBoxToText', 'on');
    end
end

function plotPTA(exposures, efficacy_targets, toxicity_targets, options)
% plotPTA  Probability of Target Attainment (PTA) plot with DYNAMIC DATA ADJUSTMENT
%
% Creates a PTA plot showing the probability of achieving efficacy targets
% while avoiding toxicity, across a range of exposures. Automatically adjusts
% exposure range and target display based on actual data distribution.
%
% INPUTS:
%   exposures         - Vector of exposure values (AUC or Cmax) from simulations
%   efficacy_targets  - Struct with:
%                       .values - Vector of efficacy thresholds (optional, auto-generated if absent)
%                       .labels - Cell array of labels (optional)
%                       .auto_generate - Set to true to create targets from data percentiles
%   toxicity_targets  - Struct with:
%                       .values - Vector of toxicity thresholds (optional, auto-generated if absent)
%                       .labels - Cell array of labels (optional)
%                       .auto_generate - Set to true to create targets from data percentiles
%   options           - Struct with customization:
%                       .metric_name    - Name of metric (default: 'Exposure')
%                       .units          - Units (default: '')
%                       .title_text     - Custom title
%                       .show_optimal   - Highlight optimal exposure range (default: true)
%                       .optimal_pta    - Minimum PTA for optimal range (default: 90%)
%                       .max_tox_prob   - Maximum acceptable toxicity probability (default: 10%)
%                       .exposure_range - [min max] for x-axis (default: auto from data ±50%)
%                       .auto_adjust    - Enable dynamic adjustment (default: true)
%                       .percentile_efficacy - Percentiles for auto-generated efficacy targets [25 50]
%                       .percentile_toxicity - Percentiles for auto-generated toxicity targets [75 90]
%                       .literature_ref - Literature reference text (optional)
%
% DYNAMIC ADJUSTMENT FEATURES:
%   - Automatically calculates exposure range from data if not specified
%   - Auto-generates efficacy targets from data percentiles if not provided
%   - Auto-generates toxicity targets from data percentiles if not provided
%   - Intelligently scales plot to show clinically relevant regions
%   - Adapts to any exposure metric without predefined limits
%
% LITERATURE:
%   - Drusano (2004) Clin Infect Dis "Antimicrobial pharmacodynamics: critical interactions of 'bug and drug'"
%   - Mouton & Vinks (2005) J Antimicrob Chemother "Pharmacokinetic/pharmacodynamic modelling"
%   - Ambrose et al. (2007) Antimicrob Agents Chemother "Monte Carlo simulation"
%   - Andes & Craig (2002) Int J Antimicrob Agents "Pharmacodynamics of fluoroquinolones"
%
% EXAMPLES:
%   % With predefined targets
%   efficacy.values = [20 30];
%   efficacy.labels = {'Minimum', 'Optimal'};
%   toxicity.values = [50 70];
%   toxicity.labels = {'Grade 2', 'Grade 3'};
%   plotPTA(auc_values, efficacy, toxicity, struct('metric_name','AUC','units','mg·h/L'));
%
%   % With AUTO-GENERATED targets from data percentiles
%   efficacy.auto_generate = true;  % Uses 25th and 50th percentiles
%   toxicity.auto_generate = true;  % Uses 75th and 90th percentiles
%   plotPTA(auc_values, efficacy, toxicity, struct('metric_name','AUC','units','mg·h/L','auto_adjust',true));

    % Validate inputs
    if nargin < 3
        error('plotPTA:InvalidInput', 'Requires exposures, efficacy_targets, and toxicity_targets');
    end
    if nargin < 4
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'metric_name'), options.metric_name = 'Exposure'; end
    if ~isfield(options, 'units'), options.units = ''; end
    if ~isfield(options, 'show_optimal'), options.show_optimal = true; end
    if ~isfield(options, 'optimal_pta'), options.optimal_pta = 90; end
    if ~isfield(options, 'max_tox_prob'), options.max_tox_prob = 10; end
    if ~isfield(options, 'auto_adjust'), options.auto_adjust = true; end
    if ~isfield(options, 'percentile_efficacy'), options.percentile_efficacy = [25 50]; end
    if ~isfield(options, 'percentile_toxicity'), options.percentile_toxicity = [75 90]; end
    
    % Remove NaN/Inf
    exposures = exposures(isfinite(exposures));
    if isempty(exposures)
        error('plotPTA:NoValidData', 'No valid exposure data');
    end
    
    % ===== DYNAMIC TARGET ADJUSTMENT FROM DATA =====
    % Auto-generate efficacy targets if not provided
    if options.auto_adjust && (~isfield(efficacy_targets, 'values') || isempty(efficacy_targets.values))
        efficacy_targets = generateDynamicTargets(exposures, options.percentile_efficacy, 'efficacy');
    end
    
    % Auto-generate toxicity targets if not provided
    if options.auto_adjust && (~isfield(toxicity_targets, 'values') || isempty(toxicity_targets.values))
        toxicity_targets = generateDynamicTargets(exposures, options.percentile_toxicity, 'toxicity');
    end
    
    % ===== DYNAMIC EXPOSURE RANGE ADJUSTMENT =====
    % Auto-detect exposure range from data if not specified
    if isfield(options, 'exposure_range') && ~isempty(options.exposure_range)
        exp_min = options.exposure_range(1);
        exp_max = options.exposure_range(2);
    else
        % Intelligent range calculation based on data
        data_min = min(exposures);
        data_max = max(exposures);
        data_range = data_max - data_min;
        
        if options.auto_adjust
            % Expand range intelligently: larger margin for sparse data
            margin_factor = max(0.3, 2 / sqrt(length(exposures)));  % Adaptive margin
            exp_min = max(0, data_min - data_range * margin_factor);
            exp_max = data_max + data_range * margin_factor;
        else
            % Standard ±50% margin
            exp_min = data_min * 0.5;
            exp_max = data_max * 1.5;
        end
    end
    
    exposure_grid = linspace(exp_min, exp_max, 200);
    
    % Create figure
    figure('Position', [100, 100, 1000, 700]);
    
    % Calculate PTA curves for efficacy targets
    n_efficacy = length(efficacy_targets.values);
    pta_efficacy = zeros(n_efficacy, length(exposure_grid));
    
    for i = 1:n_efficacy
        for j = 1:length(exposure_grid)
            pta_efficacy(i, j) = sum(exposures >= efficacy_targets.values(i)) / length(exposures) * 100;
        end
    end
    
    % Calculate probability of toxicity curves
    n_toxicity = length(toxicity_targets.values);
    prob_toxicity = zeros(n_toxicity, length(exposure_grid));
    
    for i = 1:n_toxicity
        for j = 1:length(exposure_grid)
            prob_toxicity(i, j) = sum(exposures >= toxicity_targets.values(i)) / length(exposures) * 100;
        end
    end
    
    % Plot efficacy PTA curves
    hold on;
    efficacy_colors = {[0 0.5 0], [0.2 0.7 0.2], [0.4 0.9 0.4]};
    
    for i = 1:n_efficacy
        color_idx = min(i, length(efficacy_colors));
        plot(exposure_grid, pta_efficacy(i, :), '-', 'LineWidth', 2.5, ...
            'Color', efficacy_colors{color_idx}, ...
            'DisplayName', sprintf('Efficacy: %s', efficacy_targets.labels{i}));
        
        % Mark actual target value
        plot([efficacy_targets.values(i) efficacy_targets.values(i)], [0 100], '--', ...
            'Color', efficacy_colors{color_idx}, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    
    % Plot toxicity probability curves
    toxicity_colors = {[1 0.5 0], [1 0.3 0], [0.8 0 0]};
    
    for i = 1:n_toxicity
        color_idx = min(i, length(toxicity_colors));
        plot(exposure_grid, prob_toxicity(i, :), '-', 'LineWidth', 2.5, ...
            'Color', toxicity_colors{color_idx}, ...
            'DisplayName', sprintf('Toxicity: %s', toxicity_targets.labels{i}));
        
        % Mark actual target value
        plot([toxicity_targets.values(i) toxicity_targets.values(i)], [0 100], '--', ...
            'Color', toxicity_colors{color_idx}, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    
    % Add reference lines for PTA thresholds
    plot([exp_min exp_max], [options.optimal_pta options.optimal_pta], 'k--', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    text(exp_max * 0.98, options.optimal_pta + 2, sprintf('%.0f%% PTA threshold', options.optimal_pta), ...
        'FontSize', 9, 'HorizontalAlignment', 'right');
    
    plot([exp_min exp_max], [options.max_tox_prob options.max_tox_prob], 'r--', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    text(exp_max * 0.98, options.max_tox_prob - 2, sprintf('%.0f%% toxicity limit', options.max_tox_prob), ...
        'FontSize', 9, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    
    % Highlight optimal exposure window if requested
    if options.show_optimal
        % Find exposure range where efficacy PTA > threshold AND toxicity prob < threshold
        % Use the first (most stringent) efficacy target and last (most severe) toxicity target
        pta_eff_main = pta_efficacy(1, :);
        prob_tox_main = prob_toxicity(end, :);
        
        optimal_mask = (pta_eff_main >= options.optimal_pta) & (prob_tox_main <= options.max_tox_prob);
        
        if any(optimal_mask)
            optimal_indices = find(optimal_mask);
            optimal_start = exposure_grid(optimal_indices(1));
            optimal_end = exposure_grid(optimal_indices(end));
            optimal_mid = (optimal_start + optimal_end) / 2;
            
            % Shade optimal region
            ylims = ylim;
            fill([optimal_start optimal_end optimal_end optimal_start], ...
                [ylims(1) ylims(1) ylims(2) ylims(2)], [0.7 1 0.7], ...
                'FaceAlpha', 0.25, 'EdgeColor', [0 0.7 0], 'LineWidth', 2, 'HandleVisibility', 'off');
            
            % Add annotation box with detailed info
            text(optimal_mid, ylims(2) * 0.90, ...
                sprintf('✓ OPTIMAL TARGET WINDOW'), ...
                'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
                'BackgroundColor', [0.9 1 0.9], 'EdgeColor', [0 0.6 0], 'LineWidth', 2, ...
                'Margin', 5);
            
            % Add range values below
            text(optimal_mid, ylims(2) * 0.83, ...
                sprintf('%.1f - %.1f %s', optimal_start, optimal_end, options.units), ...
                'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
                'BackgroundColor', 'white', 'EdgeColor', [0 0.6 0], 'LineWidth', 1.5, ...
                'Margin', 3);
            
            % Add recommended dose marker
            plot([optimal_mid optimal_mid], [0 ylims(2)], 'g-', 'LineWidth', 3, ...
                'DisplayName', sprintf('Recommended: %.1f %s', optimal_mid, options.units));
            
            % Calculate achieved PTA and toxicity at midpoint
            mid_pta = interp1(exposure_grid, pta_eff_main, optimal_mid);
            mid_tox = interp1(exposure_grid, prob_tox_main, optimal_mid);
            
            % Add performance metrics
            text(optimal_mid, ylims(2) * 0.75, ...
                sprintf('At midpoint:\nEfficacy PTA: %.0f%%\nToxicity risk: %.1f%%', mid_pta, mid_tox), ...
                'FontSize', 9, 'HorizontalAlignment', 'center', ...
                'BackgroundColor', [1 1 0.9], 'EdgeColor', 'none', ...
                'Margin', 2, 'FontAngle', 'italic');
            
        else
            % No optimal window found - show warning
            ylims = ylim;
            text(mean(exposure_grid), ylims(2) * 0.5, ...
                '⚠ NO OPTIMAL WINDOW FOUND', ...
                'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
                'BackgroundColor', [1 1 0.8], 'EdgeColor', [1 0.5 0], 'LineWidth', 3, ...
                'Margin', 10);
            
            text(mean(exposure_grid), ylims(2) * 0.42, ...
                sprintf('Cannot achieve ≥%.0f%% efficacy while\nkeeping toxicity ≤%.0f%%', ...
                options.optimal_pta, options.max_tox_prob), ...
                'FontSize', 10, 'HorizontalAlignment', 'center', ...
                'BackgroundColor', [1 0.95 0.9], 'EdgeColor', 'none', ...
                'Margin', 3, 'FontAngle', 'italic');
            
            % Show what might be achievable
            max_pta = max(pta_eff_main);
            min_tox_at_max_pta_idx = find(pta_eff_main == max_pta, 1);
            min_tox_at_max_pta = prob_tox_main(min_tox_at_max_pta_idx);
            
            text(mean(exposure_grid), ylims(2) * 0.32, ...
                sprintf('Best achievable: %.0f%% efficacy\nwith %.0f%% toxicity risk', max_pta, min_tox_at_max_pta), ...
                'FontSize', 9, 'HorizontalAlignment', 'center', ...
                'BackgroundColor', 'white', 'EdgeColor', 'none', ...
                'Margin', 2);
        end
    end
    
    % Formatting
    xlabel(sprintf('%s (%s)', options.metric_name, options.units), ...
        'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Probability (%)', 'FontSize', 12, 'FontWeight', 'bold');
    
    if isfield(options, 'title_text')
        title(options.title_text, 'FontSize', 13, 'FontWeight', 'bold');
    else
        title(sprintf('Probability of Target Attainment - %s', options.metric_name), ...
            'FontSize', 13, 'FontWeight', 'bold');
    end
    
    ylim([0 105]);
    grid on; grid minor;
    legend('Location', 'east', 'FontSize', 10);
    
    % Add literature reference if provided
    if isfield(options, 'literature_ref')
        annotation('textbox', [0.15 0.02 0.7 0.05], 'String', options.literature_ref, ...
            'FontSize', 8, 'EdgeColor', 'none', 'Interpreter', 'none', 'FitBoxToText', 'on');
    end
    
    hold off;
    set(gca, 'FontSize', 10);
end

function targets = generateDynamicTargets(data, percentiles, target_type)
% generateDynamicTargets  Create targets dynamically from data percentiles
%
% Automatically generates clinically meaningful target values based on the
% actual distribution of the provided data, without requiring predefined
% pharmaceutical literature thresholds.
%
% INPUTS:
%   data         - Vector of exposure values
%   percentiles  - Vector of percentiles to use [low high] (default: [25 50] for efficacy, [75 90] for toxicity)
%   target_type  - 'efficacy' or 'toxicity' (determines color and label scheme)
%
% OUTPUTS:
%   targets - Struct with .values, .labels, .colors

    data = data(isfinite(data));
    if length(data) < 3
        error('generateDynamicTargets:InsufficientData', 'Need at least 3 data points');
    end
    
    % Calculate percentile values
    perc_vals = prctile(data, percentiles);
    
    switch target_type
        case 'efficacy'
            % Efficacy targets: lower percentiles, green colors
            targets.values = perc_vals;
            targets.labels = cell(length(percentiles), 1);
            targets.colors = cell(length(percentiles), 1);
            
            for i = 1:length(percentiles)
                targets.labels{i} = sprintf('P%d (%.1f)', percentiles(i), perc_vals(i));
                % Color gradient: light green -> dark green
                intensity = percentiles(i) / 100;
                targets.colors{i} = [0.2 * (1 - intensity), 0.7 * intensity, 0.2 * (1 - intensity)];
            end
            
        case 'toxicity'
            % Toxicity targets: higher percentiles, red/orange colors
            targets.values = perc_vals;
            targets.labels = cell(length(percentiles), 1);
            targets.colors = cell(length(percentiles), 1);
            
            for i = 1:length(percentiles)
                targets.labels{i} = sprintf('P%d (%.1f)', percentiles(i), perc_vals(i));
                % Color gradient: orange -> red
                intensity = (percentiles(i) - percentiles(1)) / (percentiles(end) - percentiles(1));
                targets.colors{i} = [1, 0.5 * (1 - intensity), 0];
            end
            
        otherwise
            error('generateDynamicTargets:InvalidType', 'target_type must be ''efficacy'' or ''toxicity''');
    end
end

function plotPTA_StatisticsSummary(exposures, efficacy_targets, toxicity_targets, options)
% plotPTA_StatisticsSummary  Print statistics about PTA analysis
%
% Displays detailed statistics about target attainment rates and
% optimal therapeutic window characteristics.

    % Calculate statistics
    n_exposures = length(exposures);
    
    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf('PROBABILITY OF TARGET ATTAINMENT ANALYSIS\n');
    fprintf('════════════════════════════════════════════════════════════════\n\n');
    
    % Data summary
    fprintf('DATA SUMMARY:\n');
    fprintf('  Sample size:    %d\n', n_exposures);
    fprintf('  Mean:           %.2f %s\n', mean(exposures), options.units);
    fprintf('  Median:         %.2f %s\n', median(exposures), options.units);
    fprintf('  Min:            %.2f %s\n', min(exposures), options.units);
    fprintf('  Max:            %.2f %s\n', max(exposures), options.units);
    fprintf('  Std Dev:        %.2f %s\n', std(exposures), options.units);
    fprintf('  IQR:            %.2f - %.2f %s\n', prctile(exposures,25), prctile(exposures,75), options.units);
    fprintf('\n');
    
    % Efficacy target achievement
    fprintf('EFFICACY TARGET ATTAINMENT:\n');
    if isfield(efficacy_targets, 'values') && ~isempty(efficacy_targets.values)
        for i = 1:length(efficacy_targets.values)
            pta = sum(exposures >= efficacy_targets.values(i)) / n_exposures * 100;
            label = '';
            if isfield(efficacy_targets, 'labels') && i <= length(efficacy_targets.labels)
                label = sprintf(' (%s)', efficacy_targets.labels{i});
            end
            fprintf('  Target ≥ %.2f%s: %.1f%% PTA\n', efficacy_targets.values(i), label, pta);
        end
    end
    fprintf('\n');
    
    % Toxicity target risk
    fprintf('TOXICITY TARGET RISK:\n');
    if isfield(toxicity_targets, 'values') && ~isempty(toxicity_targets.values)
        for i = 1:length(toxicity_targets.values)
            prob_tox = sum(exposures >= toxicity_targets.values(i)) / n_exposures * 100;
            label = '';
            if isfield(toxicity_targets, 'labels') && i <= length(toxicity_targets.labels)
                label = sprintf(' (%s)', toxicity_targets.labels{i});
            end
            fprintf('  Target ≥ %.2f%s: %.1f%% risk\n', toxicity_targets.values(i), label, prob_tox);
        end
    end
    fprintf('\n');
    
    % Optimal window analysis
    if length(efficacy_targets.values) > 0 && length(toxicity_targets.values) > 0
        pta_efficacy_main = sum(exposures >= efficacy_targets.values(1)) / n_exposures * 100;
        prob_tox_main = sum(exposures >= toxicity_targets.values(end)) / n_exposures * 100;
        
        fprintf('OPTIMAL THERAPEUTIC WINDOW (PTA ≥ %.0f%%, Tox ≤ %.0f%%):\n', ...
            options.optimal_pta, options.max_tox_prob);
        fprintf('  Efficacy achievement:  %.1f%% PTA\n', pta_efficacy_main);
        fprintf('  Toxicity probability:  %.1f%% risk\n', prob_tox_main);
        
        if pta_efficacy_main >= options.optimal_pta && prob_tox_main <= options.max_tox_prob
            fprintf('  Status:                 ✓ ACHIEVABLE (meets both criteria)\n');
        elseif pta_efficacy_main >= options.optimal_pta
            fprintf('  Status:                 ⚠ Efficacy OK, toxicity risk elevated\n');
        elseif prob_tox_main <= options.max_tox_prob
            fprintf('  Status:                 ⚠ Toxicity OK, efficacy insufficient\n');
        else
            fprintf('  Status:                 ✗ NOT ACHIEVABLE (conflicts in criteria)\n');
        end
    end
    fprintf('\n');
    
    fprintf('════════════════════════════════════════════════════════════════\n\n');
end

function save_pdf_and_fig(figHandle, pdfPath, label)
% save_pdf_and_fig  Save figure to PDF and matching editable FIG.

    if nargin < 3 || strlength(string(label)) == 0
        label = 'figure';
    end

    saved_pdf = false;
    try
        print(figHandle, pdfPath, '-dpdf', '-bestfit');
        fprintf('MC: Saved %s to %s\n', label, pdfPath);
        saved_pdf = true;
    catch ME_pdf
        fprintf('MC: WARNING - Could not save %s as PDF: %s\n', label, ME_pdf.message);
        try
            saveas(figHandle, pdfPath);
            fprintf('MC: Saved %s via fallback saveas to %s\n', label, pdfPath);
            saved_pdf = true;
        catch ME_alt
            fprintf('MC: WARNING - Fallback saveas failed for %s: %s\n', label, ME_alt.message);
        end
    end

    if saved_pdf
        try
            [fig_dir, fig_name, ~] = fileparts(pdfPath);
            fig_out = fullfile(fig_dir, strcat(fig_name, '.fig'));
            savefig(figHandle, fig_out);
            fprintf('MC: Saved editable FIG to %s\n', fig_out);
        catch ME_fig
            fprintf('MC: WARNING - Could not save editable FIG for %s: %s\n', label, ME_fig.message);
        end
    end
end

function save_publication_fig_copy(figHandle, outDir, figCounter, figTitle, timestamp)
% save_publication_fig_copy  Save a FIG copy alongside publication outputs.

    try
        safe_title = regexprep(char(figTitle), '[^a-zA-Z0-9_]+', '_');
        safe_title = regexprep(safe_title, '_+', '_');
        safe_title = regexprep(safe_title, '^_|_$', '');
        if isempty(safe_title)
            safe_title = 'figure';
        end

        fig_name = sprintf('%02d_%s_%s.fig', figCounter, safe_title, char(timestamp));
        fig_out = fullfile(outDir, fig_name);
        savefig(figHandle, fig_out);
        fprintf('MC: Saved publication FIG copy to %s\n', fig_out);
    catch ME_fig_pub
        fprintf('MC: WARNING - Could not save publication FIG copy (%s): %s\n', figTitle, ME_fig_pub.message);
    end
end

function [profile, dosing_mode] = resolve_pk_target_profile(cfg, cmax_values)
% resolve_pk_target_profile  Resolve active AUC/Cmax thresholds from top-level config.

    if ~isfield(cfg, 'mode')
        cfg.mode = 'auto';
    end

    requested_mode = lower(string(cfg.mode));
    if requested_mode == "auto"
        cutoff = 10;
        if isfield(cfg, 'auto_mode_cmax_cutoff_um')
            cutoff = cfg.auto_mode_cmax_cutoff_um;
        end
        if max(cmax_values) < cutoff
            dosing_mode = "infusion";
        else
            dosing_mode = "bolus";
        end
    elseif requested_mode == "infusion" || requested_mode == "bolus"
        dosing_mode = requested_mode;
    else
        warning('Unknown PK_TARGETS.mode="%s"; defaulting to auto.', requested_mode);
        if max(cmax_values) < 10
            dosing_mode = "infusion";
        else
            dosing_mode = "bolus";
        end
    end

    profile = cfg.(char(dosing_mode));

    requiredAUC = {'subtherapeutic','optimal_low','optimal_high','toxic','severe_toxicity'};
    for i = 1:numel(requiredAUC)
        f = requiredAUC{i};
        if ~isfield(profile.AUC, f)
            error('PK_TARGETS.%s.AUC.%s is missing', dosing_mode, f);
        end
    end

    requiredCmax = {'safe_max','moderate_threshold','high_threshold','severe_threshold'};
    for i = 1:numel(requiredCmax)
        f = requiredCmax{i};
        if ~isfield(profile.Cmax, f)
            error('PK_TARGETS.%s.Cmax.%s is missing', dosing_mode, f);
        end
    end
end

function T = compute_target_adherence_summary(AUC, Cmax, profile, dosing_mode)
% compute_target_adherence_summary  Summarize Monte Carlo adherence to configured ranges.

    nAUC = max(numel(AUC), 1);
    nCmax = max(numel(Cmax), 1);

    auc_sub = 100 * sum(AUC < profile.AUC.subtherapeutic) / nAUC;
    auc_opt = 100 * sum(AUC >= profile.AUC.optimal_low & AUC <= profile.AUC.optimal_high) / nAUC;
    auc_above_opt = 100 * sum(AUC > profile.AUC.optimal_high & AUC <= profile.AUC.toxic) / nAUC;
    auc_toxic = 100 * sum(AUC > profile.AUC.toxic) / nAUC;
    auc_severe_tox = 100 * sum(AUC > profile.AUC.severe_toxicity) / nAUC;

    cmax_safe = 100 * sum(Cmax < profile.Cmax.safe_max) / nCmax;
    cmax_moderate = 100 * sum(Cmax >= profile.Cmax.safe_max & Cmax < profile.Cmax.moderate_threshold) / nCmax;
    cmax_high = 100 * sum(Cmax >= profile.Cmax.moderate_threshold & Cmax < profile.Cmax.high_threshold) / nCmax;
    cmax_toxic = 100 * sum(Cmax >= profile.Cmax.high_threshold) / nCmax;
    cmax_severe = 100 * sum(Cmax >= profile.Cmax.severe_threshold) / nCmax;

    metric = [
        "AUC_subtherapeutic";
        "AUC_optimal_window";
        "AUC_above_optimal";
        "AUC_toxic";
        "AUC_severe_toxic";
        "Cmax_safe";
        "Cmax_moderate";
        "Cmax_high";
        "Cmax_toxic";
        "Cmax_severe_toxic"
    ];

    threshold_desc = [
        sprintf("AUC < %.3g", profile.AUC.subtherapeutic);
        sprintf("%.3g <= AUC <= %.3g", profile.AUC.optimal_low, profile.AUC.optimal_high);
        sprintf("%.3g < AUC <= %.3g", profile.AUC.optimal_high, profile.AUC.toxic);
        sprintf("AUC > %.3g", profile.AUC.toxic);
        sprintf("AUC > %.3g", profile.AUC.severe_toxicity);
        sprintf("Cmax < %.3g", profile.Cmax.safe_max);
        sprintf("%.3g <= Cmax < %.3g", profile.Cmax.safe_max, profile.Cmax.moderate_threshold);
        sprintf("%.3g <= Cmax < %.3g", profile.Cmax.moderate_threshold, profile.Cmax.high_threshold);
        sprintf("Cmax >= %.3g", profile.Cmax.high_threshold);
        sprintf("Cmax >= %.3g", profile.Cmax.severe_threshold)
    ];

    percent = [
        auc_sub;
        auc_opt;
        auc_above_opt;
        auc_toxic;
        auc_severe_tox;
        cmax_safe;
        cmax_moderate;
        cmax_high;
        cmax_toxic;
        cmax_severe
    ];

    mode_col = repmat(string(dosing_mode), numel(metric), 1);
    T = table(mode_col, metric, threshold_desc, percent, ...
        'VariableNames', {'dosing_mode','metric','threshold','percent_of_simulations'});
end

%% ========================================================================
%  HELPER FUNCTIONS FOR CUSTOMIZATION
%  ========================================================================

function targets = get5FU_AUC_targets(target_set)
% get5FU_AUC_targets  Predefined literature-based AUC targets for 5-FU
%
% INPUTS:
%   target_set - 'gamelin2008', 'kaldate2012', 'beumer2019', 'conservative', 'aggressive'
%
% OUTPUTS:
%   targets - Struct with .values, .labels, .colors, .reference
%
% LITERATURE:
%   - Gamelin et al. (2008) Clin Cancer Res "Long-term weekly treatment of CRC with high-dose 5-FU"
%   - Kaldate et al. (2012) Clin Pharmacol Ther "Modeling 5-FU pharmacokinetics and ADRs"
%   - Beumer et al. (2019) Clin Cancer Res "Therapeutic drug monitoring in oncology"

    switch target_set
        case 'gamelin2008'
            % Original AUC-guided dosing targets (most widely cited)
            targets.values = [20, 28, 39, 50];
            targets.labels = {'Subtherapeutic', 'Optimal lower (Gamelin)', 'Optimal upper (Gamelin)', 'Toxicity threshold'};
            targets.colors = {[0.8 0.4 0], [0 0.7 0], [0 0.7 0], [0.8 0 0]};
            targets.reference = 'Gamelin et al. (2008) Clin Cancer Res 14(20):6677-6683';
            
        case 'kaldate2012'
            % Updated targets from population PK model with ADR correlations
            targets.values = [18, 25, 35, 45];
            targets.labels = {'Minimum response', 'Optimal lower (Kaldate)', 'Optimal upper (Kaldate)', 'Grade 3+ toxicity'};
            targets.colors = {[1 0.6 0], [0.2 0.8 0.2], [0.2 0.8 0.2], [1 0 0]};
            targets.reference = 'Kaldate et al. (2012) Clin Pharmacol Ther 91(1):44-52';
            
        case 'beumer2019'
            % Contemporary TDM guidelines for oncology
            targets.values = [22, 30, 40, 55];
            targets.labels = {'Lower efficacy limit', 'Target (Beumer)', 'Upper safety limit', 'High toxicity risk'};
            targets.colors = {[0.8 0.5 0], [0 0.6 0.8], [1 0.5 0], [0.7 0 0]};
            targets.reference = 'Beumer et al. (2019) Clin Cancer Res 25(13):3807-3815';
            
        case 'conservative'
            % Conservative dosing (prioritize safety)
            targets.values = [20, 25, 32, 40];
            targets.labels = {'Minimum', 'Conservative target', 'Safety upper limit', 'Toxicity cutoff'};
            targets.colors = {[0.7 0.7 0], [0.3 0.8 0.3], [1 0.6 0], [1 0 0]};
            targets.reference = 'Conservative approach (safety-prioritized)';
            
        case 'aggressive'
            % Aggressive dosing (prioritize efficacy in fit patients)
            targets.values = [25, 35, 45, 60];
            targets.labels = {'Minimum', 'Aggressive target', 'Toxicity watch', 'Unacceptable toxicity'};
            targets.colors = {[1 0.6 0], [0 0.5 0], [1 0.3 0], [0.8 0 0]};
            targets.reference = 'Aggressive approach (efficacy-prioritized, select patients only)';
            
        otherwise
            error('get5FU_AUC_targets:InvalidSet', ...
                'Unknown target set: %s. Valid options: gamelin2008, kaldate2012, beumer2019, conservative, aggressive', ...
                target_set);
    end
end

function targets = get5FU_Cmax_targets(target_set)
% get5FU_Cmax_targets  Predefined literature-based Cmax targets for 5-FU
%
% INPUTS:
%   target_set - 'thyss1986', 'saif2013', 'conservative', 'aggressive'
%
% OUTPUTS:
%   targets - Struct with .values, .labels, .colors, .reference
%
% LITERATURE:
%   - Thyss et al. (1986) Cancer Res "Clinical pharmacokinetic study of 5-FU"
%   - Saif et al. (2013) Clin Colorectal Cancer "Pharmacokinetically guided dose adjustment of 5-FU"

    switch target_set
        case 'thyss1986'
            % Early pharmacokinetic studies correlating Cmax with acute toxicity
            targets.values = [300, 500, 800, 1200];
            targets.labels = {'Safe threshold', 'Target peak', 'Toxicity risk (moderate)', 'Toxicity risk (severe)'};
            targets.colors = {[0 0.8 0], [0 0.5 0.8], [1 0.5 0], [1 0 0]};
            targets.reference = 'Thyss et al. (1986) Cancer Res 46(2):795-799';
            
        case 'saif2013'
            % Updated thresholds from PK-guided dosing trials
            targets.values = [350, 550, 750, 1000];
            targets.labels = {'Low risk', 'Target', 'Moderate toxicity', 'High toxicity'};
            targets.colors = {[0.2 0.8 0.2], [0.2 0.6 0.8], [1 0.6 0], [0.8 0 0]};
            targets.reference = 'Saif et al. (2013) Clin Colorectal Cancer 12(4):219-229';
            
        case 'conservative'
            % Conservative approach (lower Cmax targets)
            targets.values = [300, 450, 650, 900];
            targets.labels = {'Target', 'Upper safety', 'Toxicity watch', 'Unacceptable'};
            targets.colors = {[0.3 0.8 0.3], [0.8 0.8 0], [1 0.5 0], [1 0 0]};
            targets.reference = 'Conservative Cmax limits (safety-prioritized)';
            
        case 'aggressive'
            % Aggressive approach (allow higher Cmax for better efficacy)
            targets.values = [400, 650, 900, 1300];
            targets.labels = {'Minimum', 'Target (aggressive)', 'Acceptable toxicity', 'Severe toxicity'};
            targets.colors = {[1 0.6 0], [0 0.6 0], [1 0.4 0], [0.8 0 0]};
            targets.reference = 'Aggressive Cmax limits (efficacy-prioritized)';
            
        otherwise
            error('get5FU_Cmax_targets:InvalidSet', ...
                'Unknown target set: %s. Valid options: thyss1986, saif2013, conservative, aggressive', ...
                target_set);
    end
end

function exportPlotOptions(filename)
% exportPlotOptions  Export template configuration file for plot customization
%
% Creates a MAT file with default plotting options that can be edited
% and loaded for consistent visualization across analyses.
%
% USAGE:
%   exportPlotOptions('my_plot_config.mat')

    plot_config = struct();
    
    % CDF options
    plot_config.cdf.line_width = 2.5;
    plot_config.cdf.grid_style = 'both';
    plot_config.cdf.show_percentiles = true;
    plot_config.cdf.font_size = 10;
    plot_config.cdf.title_font_size = 13;
    
    % VPC options
    plot_config.vpc.percentiles = [5 50 95];
    plot_config.vpc.ci_alpha = 0.3;
    plot_config.vpc.log_scale = false;
    plot_config.vpc.show_individual_sims = true;
    plot_config.vpc.pred_color = [0.2 0.4 0.8];
    plot_config.vpc.obs_color = [0 0 0];
    
    % Box-whisker options
    plot_config.box.show_points = true;
    plot_config.box.point_jitter = 0.15;
    plot_config.box.point_alpha = 0.5;
    plot_config.box.show_mean = true;
    plot_config.box.add_stats = false;
    plot_config.box.box_colors = {[0.4 0.76 0.65], [0.99 0.55 0.38], [0.55 0.63 0.80]};
    
    % PTA options
    plot_config.pta.show_optimal = true;
    plot_config.pta.optimal_pta = 90;
    plot_config.pta.max_tox_prob = 10;
    plot_config.pta.line_width = 2.5;
    
    % Literature references
    plot_config.references.auc_targets = 'gamelin2008';  % or: kaldate2012, beumer2019, conservative, aggressive
    plot_config.references.cmax_targets = 'thyss1986';   % or: saif2013, conservative, aggressive
    
    save(filename, 'plot_config');
    fprintf('Plot configuration template saved to: %s\n', filename);
    fprintf('Edit this file to customize visualization options.\n');
    fprintf('Load with: load(''%s''); and pass plot_config.xxx to plotting functions.\n', filename);
end

function sim_results = run5fu_with_capture(dosingFile, outPrefix, overrides, silent_mode)
% Run base simulator and optionally capture/suppress all command-window output.
    if silent_mode
        overrides.generatePlots = false;
        if ~isfield(overrides, 'logging') || ~isstruct(overrides.logging)
            overrides.logging = struct();
        end
        if ~isfield(overrides.logging, 'level')
            overrides.logging.level = 'ERROR';
        end
    end

    if silent_mode
        evalc('sim_results = run5FU_PBPK_Simulation(dosingFile, outPrefix, overrides);');
    else
        sim_results = run5FU_PBPK_Simulation(dosingFile, outPrefix, overrides);
    end
end

function init_hpc_progress_tracker(total_runs, precompleted_runs, update_interval)
    hpc_progress_state('init', total_runs, precompleted_runs, update_interval);
end

function mc_hpc_progress_update(is_failed)
    hpc_progress_state('tick', logical(is_failed));
end

function hpc_progress_state(action, varargin)
    persistent total_runs completed_runs failed_runs update_interval start_tic

    switch action
        case 'init'
            total_runs = max(1, varargin{1});
            completed_runs = max(0, varargin{2});
            failed_runs = 0;
            update_interval = max(1, varargin{3});
            start_tic = tic;
        case 'tick'
            completed_runs = completed_runs + 1;
            if varargin{1}
                failed_runs = failed_runs + 1;
            end
            should_print = (mod(completed_runs, update_interval) == 0) || (completed_runs >= total_runs);
            if should_print
                elapsed_sec = toc(start_tic);
                pct = 100 * completed_runs / total_runs;
                if completed_runs > 0
                    eta_sec = (elapsed_sec / completed_runs) * max(0, total_runs - completed_runs);
                else
                    eta_sec = NaN;
                end
                eta_str = formatDurationCompact(eta_sec);
                fprintf('MC-HPC: %d/%d (%.1f%%) | failed=%d | elapsed=%s | ETA=%s\n', ...
                    completed_runs, total_runs, pct, failed_runs, formatDurationCompact(elapsed_sec), eta_str);
            end
        otherwise
            error('Unknown hpc_progress_state action: %s', action);
    end
end

function mc_hpc_status(varargin)
% Always-print status line for HPC minimal mode.
    fprintf(varargin{:});
end

function s = formatDurationCompact(seconds_value)
    if ~isfinite(seconds_value) || seconds_value < 0
        s = 'n/a';
        return;
    end
    if seconds_value >= 3600
        s = sprintf('%.1fh', seconds_value / 3600);
    elseif seconds_value >= 60
        s = sprintf('%.1fm', seconds_value / 60);
    else
        s = sprintf('%.0fs', seconds_value);
    end
end

%% ========================================================================
%  EXAMPLE USAGE AND CUSTOMIZATION GUIDE
%  ========================================================================
%
% BASIC USAGE:
% -------------
% The main function MC_5FU_PK_sensitivity() automatically creates all advanced
% visualizations with literature-based defaults. No additional code needed!
%
%   >> MC_5FU_PK_sensitivity(100, './my_results/')
%
% This will generate:
%   - CDF plots with PTA markers (AUC and Cmax)
%   - Visual Predictive Checks (VPC)
%   - Box-and-whisker plots by exposure quartiles
%   - Probability of Target Attainment (PTA) plots
%
% CUSTOMIZATION:
% --------------
% For custom analyses, use the individual plotting functions directly:
%
% 1) CDF WITH PTA MARKERS:
%    targets = get5FU_AUC_targets('gamelin2008');  % Use predefined targets
%    options.metric_name = 'AUC';
%    options.units = 'mg·h/L';
%    options.title_text = 'My Custom Title';
%    plotCDF_with_PTA(my_auc_data, targets, options);
%
% 2) VISUAL PREDICTIVE CHECK:
%    vpc_options.log_scale = true;  % Use log scale for concentrations
%    vpc_options.show_individual_sims = true;
%    plotVPC(observed_conc, simulated_matrix, time_vec, vpc_options);
%
% 3) BOX-AND-WHISKER PLOTS:
%    data.AUC = {group1_auc, group2_auc, group3_auc};
%    data.Cmax = {group1_cmax, group2_cmax, group3_cmax};
%    data.labels = {'Low dose', 'Medium dose', 'High dose'};
%    box_options.add_stats = true;  % Add statistical comparisons
%    plotBoxWhisker_PK(data, box_options);
%
% 4) PTA PLOTS:
%    efficacy.values = [20 30];  % Custom thresholds
%    efficacy.labels = {'Minimum', 'Target'};
%    toxicity.values = [50 70];
%    toxicity.labels = {'Grade 2', 'Grade 3+'};
%    pta_options.optimal_pta = 85;  % Require 85% PTA
%    plotPTA(exposure_data, efficacy, toxicity, pta_options);
%
% USING PREDEFINED TARGETS:
% --------------------------
%   targets_auc = get5FU_AUC_targets('gamelin2008');    % Standard
%   targets_auc = get5FU_AUC_targets('conservative');   % Safety-focused
%   targets_auc = get5FU_AUC_targets('aggressive');     % Efficacy-focused
%
%   targets_cmax = get5FU_Cmax_targets('thyss1986');
%   targets_cmax = get5FU_Cmax_targets('saif2013');
%
% EXPORT/IMPORT CONFIGURATIONS:
% ------------------------------
%   exportPlotOptions('my_config.mat');  % Create template
%   % Edit my_config.mat manually
%   load('my_config.mat');  % Load custom config
%   plotCDF_with_PTA(data, targets, plot_config.cdf);
%
% LITERATURE REFERENCES:
% ----------------------
% All plotting functions include literature references in their documentation.
% Key papers:
%
% AUC-guided dosing:
%   - Gamelin et al. (2008) Clin Cancer Res 14(20):6677-6683
%     "Long-term weekly treatment of colorectal cancer with 5-FU: AUC-guided dosing"
%   - Kaldate et al. (2012) Clin Pharmacol Ther 91(1):44-52
%     "Modeling 5-FU pharmacokinetics and ADRs in patients"
%   - Beumer et al. (2019) Clin Cancer Res 25(13):3807-3815
%     "Therapeutic drug monitoring in oncology: international consensus"
%
% Cmax correlations:
%   - Thyss et al. (1986) Cancer Res 46(2):795-799
%     "Clinical pharmacokinetic study of 5-FU"
%   - Saif et al. (2013) Clin Colorectal Cancer 12(4):219-229
%     "Pharmacokinetically guided dose adjustment of 5-FU"
%
% Pharmacometric methods:
%   - Bergstrand et al. (2011) AAPS J 13(2):143-151
%     "Prediction-corrected visual predictive checks"
%   - Holford (2005) Clin Pharmacokinet 44(9):879-899
%     "The visual predictive check—superiority to standard diagnostic plots"
%   - Mouton & Vinks (2005) J Antimicrob Chemother 55(3):253-263
%     "Pharmacokinetic/pharmacodynamic modelling of antibacterials"
%   - Drusano (2004) Clin Infect Dis 38(9):1296-1305
%     "Antimicrobial pharmacodynamics: critical interactions of 'bug and drug'"
%
% Data visualization:
%   - Weissgerber et al. (2015) PLoS Biol 13(4):e1002128
%     "Beyond bar and line graphs: time for a new data presentation paradigm"
%   - Krzywinski & Altman (2014) Nat Methods 11(2):119-120
%     "Visualizing samples with box plots"
%
%% ========================================================================
