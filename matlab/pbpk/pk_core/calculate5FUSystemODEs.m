function dCdt = calculate5FUSystemODEs(C, t, dosingRate, DPD_factor, params, time_current, dt_actual, logger)
% Calculate rates of change for all compartments and metabolites
% CORRECTED: Consistent unit handling throughout
%% =========================================================================================
%  FUNCTION: calculate5FUSystemODEs
%  =========================================================================================
%  CORE DIFFERENTIAL EQUATION SOLVER FOR 5-FU PBPK MODEL
%
%  DESCRIPTION:
%  Calculates the mass balance derivatives (dC/dt) for 5-Fluorouracil and its metabolites
%  across physiological compartments. This function represents the "physics engine" of the
%  simulation, enforcing mass conservation, saturable enzyme kinetics, and flow-limited
%  tissue transport.
%
%  MECHANISTIC OVERVIEW (THE PHYSICS):
%  1. TRANSPORT (Flow-Limited):
%     Drug moves between blood and tissues based on organ blood flow (Q).
%     Crucially, it respects the Partition Coefficient (Kp).
%     - Influx: Rate = Q * C_central
%     - Efflux: Rate = Q * (C_tissue / Kp)  <-- Ensures equilibrium at Ratio = Kp
%
%  2. METABOLISM (Hepatic DPD):
%     Primary elimination follows Michaelis-Menten kinetics (Saturable).
%     Rate = (Vmax * C) / (Km + C)
%     * Includes unit scaling to ensure Vmax matches physiological capacity (Âµmol/min).
%     * Modulated by Circadian Rhythm (DPD_factor).
%
%  3. TUMOR DYNAMICS (Anabolic):
%     Drug enters tumor via blood flow but is trapped via conversion to active
%     nucleotides (FdUMP, FdUTP, FUTP).
%     * Formation rates are scaled by Tumor Volume to ensure dimensional consistency.
%
%  DATA FLOWCHART:
%
%       [Input: Dosing]
%             â”‚
%             â–¼
%      â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”  Renal Cl   â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
%      â”‚   CENTRAL   â”‚ â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â–º â”‚  URINE   â”‚
%      â”‚ (Plasma/IV) â”‚             â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
%      â””â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”˜
%             â”‚
%       Q_liv â”‚ Kp_liv            Q_tum â”‚ Kp_tum
%             â–¼                         â–¼
%      â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”           â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
%      â”‚    LIVER    â”‚           â”‚    TUMOR    â”‚
%      â””â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”˜           â””â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”˜
%             â”‚                         â”‚
%        [DPD Enzyme]             [Anabolism]
%        (Saturable)             (Active Metab)
%             â”‚                         â”‚
%             â–¼                         â–¼
%       [Catabolites]        [FdUMP, FdUTP, FUTP]
%       (DHFU -> FBAL)       (Cytotoxic Accumulation)
%
%
%  CITATIONS & PARAMETER SOURCES:
%  [1] Diasio et al. (2001): DPD Kinetics (Km ~5 mg/L, Vmax variability).
%  [2] Schalhorn et al. (1992): Clinical PK profiles and renal clearance.
%  [3] Konings et al. (2010): Tumor partition coefficients (Kp ~0.61).
%  [4] Harris et al. (1990): Circadian rhythm of DPD activity.
%
%  INPUTS:
%    C            - Structure containing current concentrations (ÂµM)
%    t            - Current time step index
%    dosingRate   - Rate of drug input from IV (Âµmol/min)
%    DPD_factor   - Current circadian multiplier (0.6 - 1.4)
%    params       - Structure of physiological parameters (Q, V, Kp, Vmax, etc.)
%    time_current - Actual simulation time in minutes
%
%  OUTPUTS:
%    dCdt         - Structure of derivatives (Rate of change per minute)
%                   used by the Euler integration step.
%
%  =========================================================================================
% Extract current concentrations
C_central = C.C_central(t);
C_peripheral = C.C_peripheral(t);
C_liver = C.C_liver(t);
C_tumor = C.C_tumor(t);
C_FdUMP = C.C_FdUMP(t);
C_FdUTP = C.C_FdUTP(t);
C_FUTP = C.C_FUTP(t);
C_tumor_FdUMP = C.C_tumor_FdUMP(t);
C_tumor_FdUTP = C.C_tumor_FdUTP(t);
C_tumor_FUTP = C.C_tumor_FUTP(t);
C_DHFU = C.C_DHFU(t);
C_FBAL = C.C_FBAL(t);

% [Vmax Unit Scaling Check]
% Ensure Vmax is in Âµmol/min. If < 1.0, it's likely mmol/min -> convert.
if params.Vmax_DPD < 1.0
    real_Vmax = params.Vmax_DPD * 1000;
else
    real_Vmax = params.Vmax_DPD;
end
Vmax_current = real_Vmax * DPD_factor;


%% 5-FU METABOLISM - DPD-MEDIATED WITH SATURATION CHECK
% [FIX 1b: Use the already calculated Vmax_current from top of function]
% (We removed the redundant calculation here that was overwriting the fix)

%% 5-FU METABOLISM - DPD-MEDIATED WITH SATURATION & FLOW LIMITATION
MW_5FU = 130.08;
C_central_mg_L = C_central * MW_5FU / 1000; % ÂµM â†’ mg/L

% HEPATIC DPD METABOLISM
% Using the corrected Vmax_current derived at top of function
% **CRITICAL FIX**: DPD is located in hepatocytes, so use C_liver (tissue concentration),
% NOT C_central (plasma concentration). This ensures the liver compartment actually
% accumulates drug and shows realistic kinetics.
Vmax_current_mg_min = Vmax_current * MW_5FU / 1000;
C_liver_mg_L = C_liver * MW_5FU / 1000;  % Convert liver concentration to mg/L

% Michaelis-Menten kinetics (saturable) - uses LIVER concentration
rate_DPD_liver_mg_min = (Vmax_current_mg_min * C_liver_mg_L) / ...
    (params.Km_DPD + C_liver_mg_L);
rate_DPD_liver = rate_DPD_liver_mg_min / (MW_5FU / 1000); % Convert back to Âµmol/min

% SYSTEMIC (EXTRA-HEPATIC) DPD METABOLISM (20% of total)
rate_DPD_systemic_mg_min = (Vmax_current_mg_min * 0.15 * C_central_mg_L) / ...
                           (params.Km_DPD + C_central_mg_L);
rate_DPD_systemic = rate_DPD_systemic_mg_min / (MW_5FU / 1000);

% DIAGNOSTIC: Check saturation status at high concentrations
fraction_Vmax_liver = (C_liver_mg_L / (params.Km_DPD + C_liver_mg_L)) * 100;
fraction_Vmax_central = (C_central_mg_L / (params.Km_DPD + C_central_mg_L)) * 100;

% Print saturation warning at bolus peak (first 2 hours)
if mod(time_current, 60) == 0 && time_current <= 120
    if fraction_Vmax_liver > 80 || fraction_Vmax_central > 80
        fprintf('[t=%.0f min] DPD Saturation: Liver=%.0f%%, Central=%.0f%%\n', ...
            time_current, fraction_Vmax_liver, fraction_Vmax_central);
        fprintf('          (High saturation â†’ clearance drops â†’ AUC increases)\n');
    end
end

%% 5-FU DISTRIBUTION BETWEEN COMPARTMENTS (FIXED DIMENSIONALLY)
% Calculate Mass Flow Rates in Âµmol/min
% Flow = Q (L/min) * Concentration (Âµmol/L)

% Central â†” Peripheral
rate_central_to_periph = (params.k_cp * params.V_central) * C_central;
rate_periph_to_central = (params.k_pc * params.V_peripheral) * C_peripheral;

% Central â†” Liver (Blood Flow Limited)
Q_liver_L = params.Q_liver * params.CO;  % L/min
V_liver_L = params.V_liver * params.BW;  % L
rate_central_to_liver = Q_liver_L * C_central;

% [FIX 2b: Physics Calculation]
% Apply Kp to efflux. Kp=1.2 means Liver is 1.2x Blood at equilibrium.
% Rate = Flow * (Concentration / Partition_Coeff)
rate_liver_to_central = (Q_liver_L * C_liver) / params.Kp_liver;

% Central â†” Tumor (Blood Flow Limited with correct Kp application)
V_tumor_L = params.V_tumor * params.BW;  % Tumour volume in L
Q_tumor_L = params.Q_tumor * params.CO;  % L/min

rate_central_to_tumor = Q_tumor_L * C_central;

% [FIX 2c: Physics Calculation]
% Apply Kp to efflux. Kp=0.61 means Tumor is 0.61x Blood.
rate_tumor_to_central = (Q_tumor_L * C_tumor) / params.Kp_tumor;

%% ACTIVE METABOLITE FORMATION
% Formation rates (Âµmol/min) = rate constant * concentration * volume

% UMPS (uridine monophosphate synthase) catalyses formation of FdUMP
rate_form_FdUMP_sys = (params.Vmax_UMPS * C_central) / (params.Km_UMPS + C_central);

% Ribonucleotide reductase (RR) forms FUTP
rate_form_FUTP_sys = (params.Vmax_RR * C_central) / (params.Km_RR + C_central);

% dCTP deaminase forms FdUTP
rate_form_FdUTP_sys = (params.Vmax_CDA * C_central) / (params.Km_CDA + C_central);

% Adjust formation rates based on cell cycle
rate_form_FdUMP_sys = rate_form_FdUMP_sys * params.cycle_modulation_factor;
rate_form_FdUTP_sys = rate_form_FdUTP_sys * params.cycle_modulation_factor;
rate_form_FUTP_sys = rate_form_FUTP_sys * params.cycle_modulation_factor;

% Tumor metabolite formation
% [FIX 3: Dimensional Scaling]
% Multiply by V_tumor_L here to get Mass Rate (Âµmol/min)
% The ODE below divides by V_tumor_L, so units cancel to (ÂµM/min)
rate_form_FdUMP_tumor = params.k_form_FdUMP * params.tumor_metabolite_factor * C_tumor * V_tumor_L;
rate_form_FdUTP_tumor = params.k_form_FdUTP * params.tumor_metabolite_factor * C_tumor * V_tumor_L;
rate_form_FUTP_tumor = params.k_form_FUTP * params.tumor_metabolite_factor * C_tumor * V_tumor_L;

% Metabolite elimination
rate_elim_FdUMP = params.k_elim_FdUMP * C_FdUMP * params.V_central;
rate_elim_FdUTP = params.k_elim_FdUTP * C_FdUTP * params.V_central;
rate_elim_FUTP = params.k_elim_FUTP * C_FUTP * params.V_central;

% Tumor metabolite elimination
rate_elim_FdUMP_tumor = params.k_elim_FdUMP * C_tumor_FdUMP * V_tumor_L;
rate_elim_FdUTP_tumor = params.k_elim_FdUTP * C_tumor_FdUTP * V_tumor_L;
rate_elim_FUTP_tumor = params.k_elim_FUTP * C_tumor_FUTP * V_tumor_L;

sum_metabolites_sys = C_FdUMP + C_FdUTP + C_FUTP;
rate_excrete_metabolites = ((params.CL_renal_metabolites/1000) / params.V_central) * sum_metabolites_sys * params.V_central;

% Cap elimination to available mass per step to avoid negative pools
if dt_actual > 0
    avail_FdUMP = C_FdUMP * params.V_central;
    avail_FdUTP = C_FdUTP * params.V_central;
    avail_FUTP = C_FUTP * params.V_central;
    avail_FdUMP_tumor = C_tumor_FdUMP * V_tumor_L;
    avail_FdUTP_tumor = C_tumor_FdUTP * V_tumor_L;
    avail_FUTP_tumor = C_tumor_FUTP * V_tumor_L;

    rate_elim_FdUMP = min(rate_elim_FdUMP, avail_FdUMP / dt_actual);
    rate_elim_FdUTP = min(rate_elim_FdUTP, avail_FdUTP / dt_actual);
    rate_elim_FUTP = min(rate_elim_FUTP, avail_FUTP / dt_actual);
    rate_elim_FdUMP_tumor = min(rate_elim_FdUMP_tumor, avail_FdUMP_tumor / dt_actual);
    rate_elim_FdUTP_tumor = min(rate_elim_FdUTP_tumor, avail_FdUTP_tumor / dt_actual);
    rate_elim_FUTP_tumor = min(rate_elim_FUTP_tumor, avail_FUTP_tumor / dt_actual);

    avail_metabolites_sys = sum_metabolites_sys * params.V_central;
    rate_excrete_metabolites = min(rate_excrete_metabolites, avail_metabolites_sys / dt_actual);
end

%% CIRCADIAN MODULATION OF THYMIDYLATE SYNTHASE (NEW)
hourOfDay = mod(time_current, 1440) / 60;
TS_circadian = 1.0 + params.TS_acrophase * cos(2*pi*(hourOfDay - params.TS_peak_hour)/24);

% Apply circadian modulation to FdUMP formation (TS is direct target)
rate_form_FdUMP_sys = rate_form_FdUMP_sys * TS_circadian;

%% CATABOLITE FORMATION AND EXCRETION
rate_form_DHFU = rate_DPD_liver + rate_DPD_systemic;
rate_form_FBAL = params.k_DHFU_to_FBAL * C_DHFU * params.V_central;
rate_excrete_5FU = ((params.CL_renal_5FU/1000) / params.V_central) * C_central * params.V_central;
rate_excrete_FBAL = ((params.CL_renal_FBAL/1000) / params.V_central) * C_FBAL * params.V_central;
rate_clear_metabolites = rate_excrete_metabolites + rate_elim_FdUMP + rate_elim_FdUTP + rate_elim_FUTP + ...
    rate_elim_FdUMP_tumor + rate_elim_FdUTP_tumor + rate_elim_FUTP_tumor;

%% PHARMACODYNAMIC ODE CONTRIBUTIONS
% TS inhibition from tumor FdUMP drives dNTP depletion; TS_inhibition_fraction is also
% updated algebraically per step in the main loop after updateConcentrations.
TS_inhib_current = C_tumor_FdUMP / (params.IC50_TS_tumor + C_tumor_FdUMP);
dCdt.dNTP_pool_fraction = updateDNTPPool(C.dNTP_pool_fraction(t), TS_inhib_current, params);
dCdt.DNA_damage_index   = updateDNADamage(C.DNA_damage_index(t), C.dNTP_pool_fraction(t), params);


%% DIFFERENTIAL EQUATIONS
% Central compartment (5-FU)
dCdt.C_central = (dosingRate / params.V_central) + ...
    (rate_periph_to_central / params.V_central) - ...
    (rate_central_to_periph / params.V_central) + ...
    (rate_liver_to_central / params.V_central) - ...
    (rate_central_to_liver / params.V_central) + ...
    (rate_tumor_to_central / params.V_central) - ...
    (rate_central_to_tumor / params.V_central) - ...
    (rate_DPD_systemic / params.V_central) - ...
    (rate_form_FdUMP_sys / params.V_central) - ...
    (rate_form_FdUTP_sys / params.V_central) - ...
    (rate_form_FUTP_sys / params.V_central) - ...
    (rate_excrete_5FU / params.V_central);

% Peripheral compartment (5-FU)
dCdt.C_peripheral = (rate_central_to_periph / params.V_peripheral) - ...
    (rate_periph_to_central / params.V_peripheral);

% Liver compartment (5-FU)
dCdt.C_liver = (rate_central_to_liver / V_liver_L) - ...
    (rate_liver_to_central / V_liver_L) - ...
    (rate_DPD_liver / V_liver_L);

% Tumor compartment (5-FU)
% Dividing by V_tumor_L correctly converts the Mass Rate (Âµmol/min) into Concentration Rate (ÂµM/min)
dCdt.C_tumor = (rate_central_to_tumor / V_tumor_L) - ...
    (rate_tumor_to_central / V_tumor_L) - ...
    (rate_form_FdUMP_tumor / V_tumor_L) - ...
    (rate_form_FdUTP_tumor / V_tumor_L) - ...
    (rate_form_FUTP_tumor / V_tumor_L);

% Kidney, muscle, fat (simplified)
dCdt.C_kidney = 0;
dCdt.C_muscle = 0;
dCdt.C_fat = 0;

% Active metabolites - systemic
dCdt.C_FdUMP = (rate_form_FdUMP_sys / params.V_central) - ...
    (rate_elim_FdUMP / params.V_central);
dCdt.C_FdUTP = (rate_form_FdUTP_sys / params.V_central) - ...
    (rate_elim_FdUTP / params.V_central);
dCdt.C_FUTP = (rate_form_FUTP_sys / params.V_central) - ...
    (rate_elim_FUTP / params.V_central);

if sum_metabolites_sys > 0
    rate_excrete_FdUMP = rate_excrete_metabolites * (C_FdUMP / sum_metabolites_sys);
    rate_excrete_FdUTP = rate_excrete_metabolites * (C_FdUTP / sum_metabolites_sys);
    rate_excrete_FUTP = rate_excrete_metabolites * (C_FUTP / sum_metabolites_sys);
else
    rate_excrete_FdUMP = 0;
    rate_excrete_FdUTP = 0;
    rate_excrete_FUTP = 0;
end

dCdt.C_FdUMP = dCdt.C_FdUMP - (rate_excrete_FdUMP / params.V_central);
dCdt.C_FdUTP = dCdt.C_FdUTP - (rate_excrete_FdUTP / params.V_central);
dCdt.C_FUTP = dCdt.C_FUTP - (rate_excrete_FUTP / params.V_central);

% Active metabolites - tumor
dCdt.C_tumor_FdUMP = (rate_form_FdUMP_tumor - rate_elim_FdUMP_tumor) / V_tumor_L;
dCdt.C_tumor_FdUTP = (rate_form_FdUTP_tumor - rate_elim_FdUTP_tumor) / V_tumor_L;
dCdt.C_tumor_FUTP = (rate_form_FUTP_tumor - rate_elim_FUTP_tumor) / V_tumor_L;

% Catabolites
dCdt.C_DHFU = (rate_form_DHFU / params.V_central) - ...
    (rate_form_FBAL / params.V_central);
dCdt.C_FBAL = (rate_form_FBAL / params.V_central) - ...
    (rate_excrete_FBAL / params.V_central);

% Cumulative excretion
dCdt.excreted_5FU = rate_excrete_5FU;
dCdt.excreted_FBAL = rate_excrete_FBAL;
dCdt.excreted_metabolites = rate_clear_metabolites;
dCdt.cumulative_input_5FU = dosingRate;
end
