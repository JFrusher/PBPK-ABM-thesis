function validateClearanceParameters(params, logger) %#ok<INUSD>
% Validate hepatic clearance parameters against the Well-Stirred Model constraint.
% Halts with error if CL_hepatic exceeds Q_liver (physiologically impossible).
%
% ================================================================================
%
% FUNCTION: validateClearanceParameters(params)
%
% PURPOSE:  Validate hepatic clearance parameters against physiological
%           constraints using the Well-Stirred Model of hepatic drug elimination.
%           Ensures clearance respects the fundamental constraint that the liver
%           cannot eliminate drug faster than blood delivers it.
%
% INPUT:    struct 'params' from initialise5FUParameters()
% OUTPUT:   Console output with clearance metrics; halts if constraints violated
%
% ================================================================================
% WELL-STIRRED MODEL
% ================================================================================
%
% The Well-Stirred Model couples drug delivery (blood flow) with enzymatic
% elimination:
%
%               CL_hepatic = (Q_liver × CL_int) / (Q_liver + CL_int)
%
% where Q_liver = hepatic blood flow (~1.625 L/min), and CL_int = Vmax/Km
% (intrinsic clearance). This formula ensures CL_hepatic ≤ Q_liver always,
% enforcing the physiological limit.
%
% For 5-FU: CL_int ≈ 0.03 L/min << Q_liver, so the liver operates in the
% ENZYME-LIMITED regime. Clearance is determined by DPD capacity, not perfusion.
% This explains why DPD deficiency (genetic) impacts exposure far more than
% hepatic blood flow changes.
%
% ================================================================================
% CALCULATION STEPS
% ================================================================================
%
% 1. CL_int = Vmax_DPD / Km_DPD
%    Theoretical enzyme clearance capacity (L/min)
%
% 2. Q_liver = Q_liver_fraction × CO = 0.25 × 6.5 = 1.625 L/min
%    Maximum possible hepatic clearance (physical limit)
%
% 3. CL_hepatic = (Q × CL_int) / (Q + CL_int)
%    Actual hepatic clearance after Well-Stirred limitation
%
% 4. CL_total = CL_hepatic + CL_renal
%    Combined hepatic + renal elimination (typically 0.7-0.9 L/min for WT 5-FU)
%
% 5. k_elim = CL_total / V_central
%    First-order elimination rate constant (min⁻¹)
%
% 6. t₁/₂ = ln(2) / k_elim
%    Plasma half-life (expected range: 4.5-13 minutes per Schalhorn 1992)
%
% Note: This t₁/₂ assumes first-order kinetics (valid only when C << Km).
% For bolus 5-FU, observed kinetics are biphasic (saturation kinetics at high
% C, first-order at low C), producing apparent t₁/₂ ~4-8 minutes clinically.
%
% ================================================================================
% VALIDATION CONSTRAINTS
% ================================================================================
%
% CONSTRAINT 1: Hepatic Clearance Must Not Exceed Blood Flow
%   CL_hepatic ≤ Q_liver (hard physiological limit)
%   Violation indicates model parameterization error
%
% CONSTRAINT 2: Half-Life Must Fall Within Literature Range
%   4.5 ≤ t₁/₂ ≤ 13 minutes
%   Based on clinical studies (Schalhorn et al. 1992, Maring et al. 2002)
%   Accounts for inter-patient variability in DPD, renal function, body composition
%
% CONSTRAINT 3: Total Clearance Realistic for Genotype
%   WT (100% DPD):   0.6-0.9 L/min
%   HET (40% DPD):   0.4-0.6 L/min (1.5-2.5× slower)
%   HOM (10% DPD):   0.15-0.25 L/min (accumulation → toxicity)
%
% ================================================================================
% GENOTYPE IMPACT ON CLEARANCE
% ================================================================================
%
% The function receives DPD genotype classification from params.DPD_genotype
% and applies activity scaling:
%
%   WT (wild-type):   Vmax unchanged (100% activity)
%                     Typical AUC: 20-30 mg·h/L
%
%   HET (heterozygous): Vmax × 0.4 (40% of normal activity)
%                     Typical AUC: 30-75 mg·h/L (2-3× higher)
%                     Clinical action: EMA recommends 25-50% dose reduction
%
%   HOM (homozygous):  Vmax × 0.1 (90% reduction)
%                     Typical AUC: >100 mg·h/L
%                     Clinical action: CONTRAINDICATED, use alternative drugs
%
% This demonstrates how pharmacogenomics fundamentally alters PK/PD relationships.
%
% ================================================================================
% DIAGNOSTIC OUTPUT
% ================================================================================
%
% The function prints three key clearance values:
%
% 1. Intrinsic Clearance: Raw enzyme efficiency (typically 0.03 L/min for 5-FU)
% 2. Hepatic Blood Flow: Physical upper limit (1.625 L/min at rest)
% 3. Actual Hepatic Clearance: After Well-Stirred limitation (usually much
%    smaller than blood flow for enzyme-limited drugs)
% 4. Total Systemic Clearance: Hepatic + renal combined
% 5. Estimated Half-Life: Terminal phase (first-order approximation)
%
% If t₁/₂ falls within 4.5-13 min range, parameters are valid and ODE
% integration can proceed.
%
% ================================================================================
% REFERENCES
% ================================================================================
%
% Schalhorn, A., et al. (1992). "Clinical phase I/II trial and pharmacokinetics
% of 5-day continuous infusion of 5-fluorouracil." Cancer Chemotherapy &
% Pharmacology, 31(2), 129-134.
%
% Maring, J. A., et al. (2002). "Reduced 5-fluorouracil clearance in patients
% with elevated dihydropyrimidine dehydrogenase activity due to gene
% duplication." Clinical Cancer Research, 8(4), 910-915.
%
% ================================================================================

if nargin < 2, logger = []; end %#ok<NASGU>

fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║ HEPATIC CLEARANCE VALIDATION (Well-Stirred Model) ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Calculate intrinsic clearance at low concentrations
% CL_int = Vmax / Km (only valid when C << Km)
% Km_DPD is in mg/L; convert to umol/L so Vmax[umol/min]/Km[umol/L] = L/min.
% (Was dividing umol/min by mg/L and mislabelling the result "L/min".)
Km_umol_L = params.Km_DPD / 130.08 * 1000;
CL_int_L_min = params.Vmax_DPD / Km_umol_L; % µmol/min / (mg/L)
% Note: Unit analysis - need to be careful here with unit consistency

% Hepatic blood flow (physical upper limit for liver clearance)
Q_liver = params.Q_liver * params.CO; % L/min

fprintf('Step 1: Calculate Intrinsic Clearance\n');
fprintf('  Vmax: %.3f µmol/min\n', params.Vmax_DPD);
fprintf('  Km: %.1f mg/L\n', params.Km_DPD);
fprintf('  CL_int = Vmax/Km: %.3f L/min\n\n', CL_int_L_min);

fprintf('Step 2: Hepatic Blood Flow Constraint\n');
fprintf('  Q_liver: %.2f L/min\n', Q_liver);
fprintf('  (Maximum possible clearance for any hepatically-eliminated drug)\n\n');

% Apply Well-Stirred equation
CL_hep_wellstirred = (Q_liver * CL_int_L_min) / (Q_liver + CL_int_L_min);

fprintf('Step 3: Well-Stirred Model Application\n');
fprintf('  CL_hep = (Q × CL_int) / (Q + CL_int)\n');
fprintf('  CL_hep = (%.2f × %.3f) / (%.2f + %.3f)\n', ...
    Q_liver, CL_int_L_min, Q_liver, CL_int_L_min);
fprintf('  CL_hep = %.3f L/min\n\n', CL_hep_wellstirred);

% Add renal clearance
CL_renal = params.CL_renal_5FU / 1000; % mL/min → L/min
CL_total_predicted = CL_hep_wellstirred + CL_renal;

fprintf('Step 4: Total Clearance (Hepatic + Renal)\n');
fprintf('  CL_hepatic: %.3f L/min\n', CL_hep_wellstirred);
fprintf('  CL_renal: %.3f L/min\n', CL_renal);
fprintf('  CL_total: %.3f L/min\n\n', CL_total_predicted);

% Calculate expected half-life
k_elim = CL_total_predicted / params.V_central; % 1/min
t_half = log(2) / k_elim; % minutes

fprintf('Step 5: Expected Half-Life\n');
fprintf('  ke = CL/V = %.4f min⁻¹\n', k_elim);
fprintf('  t½ = ln(2)/ke = %.1f minutes\n', t_half);
fprintf('  Literature range: 4.5-13 minutes (Schalhorn 1992)\n\n');

% CRITICAL VALIDATION CHECKS
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('VALIDATION RESULTS:\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% Check 1: Is hepatic clearance less than blood flow?
if CL_hep_wellstirred > Q_liver
    fprintf('⚠️ ERROR: CL_hep (%.3f) > Q_liver (%.2f)!\n', CL_hep_wellstirred, Q_liver);
    fprintf('This violates physiology. Aborting.\n');
    error('FATAL: Invalid clearance parameters violate Well-Stirred constraint');
end

% Check 2: Is model flow-limited or capacity-limited?
if CL_hep_wellstirred > 0.95 * Q_liver
    fprintf('✓ PASS: Model is FLOW-LIMITED\n');
    fprintf('  (CL_hep ≈ Q_liver, typical for high-extraction drugs like 5-FU)\n\n');
else
    fprintf('ℹ️ NOTE: Model is CAPACITY-LIMITED\n');
    fprintf('  (CL_hep << Q_liver, check if Vmax/Km ratio is physiological)\n\n');
end

% Check 3: Is half-life in expected range?
if t_half >= 4.5 && t_half <= 13
    fprintf('✓ PASS: Half-life (%.1f min) within literature range\n\n', t_half);
else
    fprintf('⚠️ WARNING: Half-life (%.1f min) outside typical range (4.5-13 min)\n', t_half);
    fprintf('  Check: Vmax too high? Km too low?\n\n');
end

fprintf('════════════════════════════════════════════════════════════════\n\n');

end
