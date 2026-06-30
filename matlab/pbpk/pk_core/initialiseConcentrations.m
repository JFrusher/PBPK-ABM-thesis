function concentrations = initialiseConcentrations(nTimePoints)
% CONCENTRATION ARRAY INITIALISATION FUNCTION
% ================================================================================
%
% FUNCTION: initialiseConcentrations(nTimePoints)
%
% PURPOSE:  Pre-allocate concentration arrays for all PBPK compartments,
%           metabolites, and pharmacodynamic markers. All values initialized
%           to zero (baseline). Arrays are column vectors with length nTimePoints,
%           enabling efficient time-series storage during ODE integration.
%
% INPUT:    nTimePoints (integer) - Number of time steps in simulation
% OUTPUT:   struct 'concentrations' with 20 fields, each [nTimePoints × 1] vector
%
% ================================================================================
% COMPARTMENTAL STRUCTURE
% ================================================================================
%
% PARENT DRUG (5-FU) - µM units:
%   • C_central: Plasma and highly perfused organs (blood, heart, brain)
%   • C_peripheral: Muscle, adipose, skin (poorly perfused)
%   • C_liver: Hepatic tissue (primary metabolism site)
%   • C_kidney: Renal tissue (filtration and secretion)
%   • C_tumor: Malignant tissue (target compartment)
%   • C_muscle: Skeletal muscle (major distribution sink)
%   • C_fat: Adipose tissue (hydrophobic reservoir, low Kp = 0.3)
%
% ACTIVE METABOLITES - µM units (intracellular + circulating):
%   • C_FdUMP: Fluorodeoxyuridine monophosphate
%     → Inhibits thymidylate synthase (TS), blocks dTMP synthesis
%     → Primary mechanism of cytotoxicity
%
%   • C_FdUTP: Fluorodeoxyuridine triphosphate
%     → Incorporated into DNA during replication
%     → Causes replication stress and apoptosis
%
%   • C_FUTP: Fluorouridine triphosphate
%     → Incorporated into RNA
%     → Disrupts transcription and translation
%     → Longest intracellular retention (slow elimination)
%
% TUMOR-SPECIFIC METABOLITES - µM units (local concentration):
%   • C_tumor_FdUMP, C_tumor_FdUTP, C_tumor_FUTP
%   → Track metabolite accumulation in tumor tissue separately
%   → Tumor cells form metabolites 2.5× faster than normal tissue
%   → Enables PD modeling of tumor response vs. normal tissue toxicity
%
% CATABOLITES - µM units (5-FU degradation products):
%   • C_DHFU: Dihydrofluorouracil (first DPD-catalyzed product)
%   • C_FBAL: α-fluoro-β-alanine (final urinary excretion form)
%   → Track clearance pathway and metabolite fate
%
% CUMULATIVE EXCRETION - µmol total (not concentration):
%   • excreted_5FU: Unchanged 5-FU in urine (~10% of dose)
%   • excreted_FBAL: FBAL in urine (primary route, ~80-90% of dose)
%   • excreted_metabolites: Other minor pathways
%   • cumulative_input_5FU: Total IV input over time (µmol)
%
% ================================================================================
% PHARMACODYNAMIC BIOMARKERS
% ================================================================================
%
% These compartments link drug concentration to cellular effect and toxicity:
%
% 1. dNTP_POOL_FRACTION (dimensionless, range 0-1)
%    • Tracks depletion of deoxynucleotide triphosphate pool
%    • TS inhibition by FdUMP → dTTP depletion → dNTP imbalance
%    • Initialized to 1.0 (fully replenished)
%    • During therapy: decreases toward 0.5 or lower
%    • Models: C_FdUMP → ODE derivative → dNTP dynamics
%
% 2. TS_INHIBITION_FRACTION (dimensionless, range 0-1)
%    • Fraction of thymidylate synthase enzyme bound by FdUMP
%    • Initialized to 0 (no inhibition at baseline)
%    • During therapy: rises toward 0.8-0.95 at peak drug levels
%    • Driving force: C_FdUMP concentration and TS abundance
%    • Pharmacodynamic model: Emax function of C_FdUMP
%
% 3. DNA_DAMAGE_INDEX (accumulating score)
%    • Cumulative DNA damage from dNTP imbalance
%    • Caused by incorporation of dUMP into DNA (due to dTMP deficiency)
%    • Initialized to 0
%    • Increases monotonically with TS inhibition duration
%    • Links to apoptosis triggering: when index > threshold, cells die
%
% 4. S_PHASE_FRACTION (fixed parameter — not stored as a concentration array)
%    • Held constant at params.S_phase_fraction = 0.15; see Rupa et al. 2003
%    • Only S-phase cells are vulnerable to dNTP-imbalance-driven DNA damage
%
% ================================================================================
% ARRAY STRUCTURE AND ODE INTEGRATION
% ================================================================================
%
% During ode45 integration, the output Y matrix is [nTimePoints × nStates],
% where nStates = 20 (one for each concentration variable).
%
% The initialiseConcentrations output preallocates storage:
%
%   for timeStep = 1:nTimePoints
%     Y(timeStep, :) contains all 20 state variables at that time point
%     concentrations.C_central(timeStep) = Y(timeStep, 1)
%     concentrations.C_FdUMP(timeStep) = Y(timeStep, 9)
%     etc.
%
% This structure enables:
%   • Efficient memory allocation (no dynamic array growth)
%   • Post-processing visualization (plot(T, concentrations.C_central))
%   • Metabolite kinetics analysis (compare C_FdUMP vs C_FUTP)
%   • PD endpoint calculations (integral of TS_inhibition_fraction)
%
% ================================================================================
% UNIT CONVENTIONS
% ================================================================================
%
% All concentrations in µM (micromolar = µmol/L) EXCEPT:
%
%   • Cumulative excretion (excreted_*): µmol total excreted (time integral)
%   • dNTP_pool_fraction: dimensionless (fraction, 0-1)
%   • TS_inhibition_fraction: dimensionless (fraction, 0-1)
%   • DNA_damage_index: dimensionless accumulation (arbitrary units)
%   • S_phase_fraction: now a fixed parameter in params, not a concentration array
%
% Conversion: If needed to convert from mg to µmol:
%   C_µM = (C_mg/L × 1000) / MW_5FU
%   where MW_5FU = 130.08 g/mol
%
% Example: 5 mg/L = (5 × 1000) / 130.08 ≈ 38.4 µM
%
% ================================================================================
    % Parent drug (5-FU) in PBPK compartments
    concentrations.C_central = zeros(nTimePoints, 1);      % Central compartment (µM)
    concentrations.C_peripheral = zeros(nTimePoints, 1);   % Peripheral compartment (µM)
    concentrations.C_liver = zeros(nTimePoints, 1);        % Liver tissue (µM)
    concentrations.C_kidney = zeros(nTimePoints, 1);       % Kidney tissue (µM)
    concentrations.C_tumor = zeros(nTimePoints, 1);        % Tumor tissue (µM)
    concentrations.C_muscle = zeros(nTimePoints, 1);       % Muscle tissue (µM)
    concentrations.C_fat = zeros(nTimePoints, 1);          % Adipose tissue (µM)

    % Active metabolites (primarily intracellular, but tracked systemically)
    concentrations.C_FdUMP = zeros(nTimePoints, 1);        % FdUMP (TS inhibitor) (µM)
    concentrations.C_FdUTP = zeros(nTimePoints, 1);        % FdUTP (DNA incorporation) (µM)
    concentrations.C_FUTP = zeros(nTimePoints, 1);         % FUTP (RNA incorporation) (µM)

    % Tumor metabolites (local concentration in tumor)
    concentrations.C_tumor_FdUMP = zeros(nTimePoints, 1);  % Tumor FdUMP (µM)
    concentrations.C_tumor_FdUTP = zeros(nTimePoints, 1);  % Tumor FdUTP (µM)
    concentrations.C_tumor_FUTP = zeros(nTimePoints, 1);   % Tumor FUTP (µM)

    % Catabolites (degradation products)
    concentrations.C_DHFU = zeros(nTimePoints, 1);         % Dihydrofluorouracil (µM)
    concentrations.C_FBAL = zeros(nTimePoints, 1);         % α-fluoro-β-alanine (µM)

    % Cumulative excretion (µmol total excreted)
    concentrations.excreted_5FU = zeros(nTimePoints, 1);   % Unchanged 5-FU in urine
    concentrations.excreted_FBAL = zeros(nTimePoints, 1);  % FBAL in urine (main route)
    concentrations.excreted_metabolites = zeros(nTimePoints, 1); % Other metabolites
    concentrations.cumulative_input_5FU = zeros(nTimePoints, 1); % Total drug input (µmol)

    % Pharmacodynamic compartments tracking drug effect (NEW)
    concentrations.dNTP_pool_fraction = ones(nTimePoints, 1); % dNTP pool depletion (1.0 = normal, 0.5 = 50% depleted)
    concentrations.TS_inhibition_fraction = zeros(nTimePoints, 1); % Fraction of TS bound by FdUMP (0 = none, 1 = all)
    concentrations.DNA_damage_index = zeros(nTimePoints, 1); % Accumulation of damage from dNTP imbalance
    % S_phase_fraction is now a fixed parameter (params.S_phase_fraction), not a time-varying array
end
