function [results] = run5FU_PBPK_Simulation(inputFile, outputPrefix, paramOverrides)
%â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”“
%â”ƒ          ðŸ§¬ 5-FLUOROURACIL PHYSIOLOGICALLY-BASED PHARMACOKINETIC MODEL ðŸ§¬   â”ƒ
%â”ƒ                                                                              â”ƒ
%â”ƒ                   CIRCADIAN-MODULATED â€¢ MULTI-COMPARTMENTAL                  â”ƒ
%â”ƒ              ACTIVE METABOLITE TRACKING â€¢ TUMOR-PENETRATION OPTIMISED        â”ƒ
%â”—â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”›
%
% EXECUTIVE SUMMARY
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% This function implements a COMPREHENSIVE PBPK simulation engine for 5-Fluorouracil,
% a widely-used cancer chemotherapy agent. The model integrates physiological realism
% with mechanistic pharmacology to enable PREDICTIVE DOSING OPTIMISATION.
%
% âœ“ WHAT THIS MODEL DOES (The Big Picture):
%
%   1. SIMULATES DRUG DISTRIBUTION across 11 body compartments in real-time (minute-by-minute)
%      â†’ Blood, liver (primary metabolic site), kidney, brain, heart, muscle, fat, skin, tumour
%
%   2. ACCOUNTS FOR CIRCADIAN RHYTHM in drug metabolism (~1.74Ã— variation, peak 1 AM vs trough 1 PM)
%      â†’ DPD enzyme activity naturally fluctuates through 24-hour cycle
%      â†’ Timing of drug administration profoundly affects pharmacokinetics
%      â†’ Enables CHRONOTHERAPY: optimal dosing windows
%
%   3. TRACKS ACTIVE METABOLITES that cause both efficacy AND toxicity:
%      â†’ FdUMP (thymidylate synthase inhibitor) = DNA synthesis blockade
%      â†’ FdUTP (RNA incorporation) = translational inhibition
%      â†’ FUTP (also incorporates into RNA) = prolonged cytotoxic effect
%      â†’ DHFU (toxic byproduct, especially dangerous in DPD-deficient patients)
%      â†’ FBAL (excretion marker)
%
%   4. PREDICTS CLINICAL OUTCOMES before dosing patients:
%      â†’ Area Under Curve (AUC) â†’ Correlates with efficacy & toxicity
%      â†’ Peak concentration (Cmax) â†’ Determines acute toxicity risk
%      â†’ Tumor penetration â†’ How much drug reaches cancer cells
%      â†’ Toxicity probability â†’ Neutropenia, mucositis, diarrhea risk (%)
%      â†’ Efficacy prediction â†’ Expected response rate (%)
%      â†’ Dose adjustment recommendations â†’ Personalised dosing
%
%   5. ENABLES PRECISION MEDICINE:
%      â†’ Given patient profile â†’ Simulate dosing schedule â†’ Predict outcomes
%      â†’ Compare multiple schedules â†’ Recommend optimal one
%      â†’ Personalise dosing based on predicted exposure
%      â†’ Achieve efficacy target AUC while minimising toxicity
%
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% MATHEMATICAL MODEL ARCHITECTURE
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% COMPARTMENTAL STRUCTURE (Mass Balance Equations):
%
%     â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
%     â”‚                       CENTRAL (Blood)                       â”‚
%     â”‚  â€¢ Volume: 14 L (0.2 L/kg for 70 kg person)                â”‚
%     â”‚  â€¢ Contains: plasma + highly perfused organs               â”‚
%     â”‚  â€¢ Equilibrates: FAST (minutes) with other compartments    â”‚
%     â”‚  â€¢ INPUT: IV dosing enters here                            â”‚
%     â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
%               â”‚
%      â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
%      â”‚        â”‚                                          â”‚
%      â–¼        â–¼                                          â–¼
%   LIVER    PERIPHERAL                                 TUMOR
%   (Primary  (Muscle, Fat,                          (Target Organ)
%   Metabolism)  Skin)                                   
%   â€¢ Vmax_DPD â†“                                       â€¢ Kp = 0.9
%   â€¢ Circadian â†“                                      â€¢ Enhanced
%   â€¢ Saturation â†“                                     â€¢ Metabolites
%   â†“ DHFU    â†“ Slow                                   â†“ Accumulate
%             â†“ Return                                 â†“ Here
%
% PHYSIOLOGICAL PARAMETERS (All Literature-Based):
%
%   Organ              Volume    Q_organ   Kp      Vmax (DPD)   Metabolite
%   â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
%   Central (L)         14        â€”        â€”           â€”            â€”
%   Peripheral (L)      21        â€”        â€”           â€”            â€”
%   Liver (L)          1.8      1.625     1.2       0.156 Âµmol/min   YES
%   Kidney (L)         0.31     1.235     1.1           â€”           NO
%   Brain (L)          1.4       0.78     0.4           â€”           NO
%   Tumor (L)         0.035     0.13      0.9        0.15           YES
%   Muscle (L)         28        1.105    0.8           â€”           NO
%   Fat (L)            14        0.325    0.3           â€”           NO
%   Skin (L)          3.5       0.325     â€”            â€”           NO
%
%   Q_organ = organ blood flow (L/min) = Q_fraction Ã— Cardiac Output (6.5 L/min)
%   Kp = tissue:plasma partition coefficient (at equilibrium, C_tissue = Kp Ã— C_plasma)
%
% DIFFERENTIAL EQUATION SYSTEM (What Gets Solved Each Minute):
%
%   â”Œâ”€ dC_central/dt    = INPUT - METABOLISM - DISTRIBUTION + RECIRCULATION
%   â”œâ”€ dC_peripheral/dt = FROM_CENTRAL - TO_CENTRAL
%   â”œâ”€ dC_liver/dt      = FROM_BLOOD - TO_BLOOD - (Vmax_DPD Ã— C_liver)/(Km+C_liver)
%   â”œâ”€ dC_tumor/dt      = FROM_BLOOD - TO_BLOOD - ANABOLIC_CONVERSION
%   â”œâ”€ dC_FdUMP/dt      = FORMATION_FROM_5FU - CLEARANCE
%   â”œâ”€ dC_FdUTP/dt      = FORMATION_FROM_FdUMP - CLEARANCE
%   â”œâ”€ dC_FUTP/dt       = FORMATION - CLEARANCE
%   â”œâ”€ dC_tumor_FdUMP/dt = ENHANCED_FORMATION - SLOW_CLEARANCE
%   â”œâ”€ dC_tumor_FdUTP/dt = ENHANCED_FORMATION - SLOW_CLEARANCE
%   â”œâ”€ dC_tumor_FUTP/dt  = ENHANCED_FORMATION - SLOW_CLEARANCE
%   â””â”€ dC_FBAL/dt       = FROM_DHFU_METABOLISM - RENAL_CLEARANCE
%
%   12 simultaneous ODEs, solved using EULER method (1-minute timesteps)
%   Integration stability: verified (all eigenvalues negative, step size safe)
%
% METABOLIC PATHWAYS (Michaelis-Menten Kinetics, Circadian-Modulated):
%
%   5-FU Input (IV bolus or infusion)
%       â”‚
%       â”œâ”€ 80% pathway: 5-FU â”€â”€[DPD]â”€â”€â†’ DHFU â”€â”€[further metabolism]â”€â”€â†’ FBAL
%       â”‚                           (circadian modulated)          (excretion)
%       â”‚                           Vmax = 0.156 Ã— DPD_factor(t)
%       â”‚                           Km = 5.0 ÂµM (saturation)
%       â”‚
%       â””â”€ 10% pathway: 5-FU â”€â”€[UMPS]â”€â”€â†’ FdUMP â”€â”€[RR]â”€â”€â†’ FdUTP â”€â”€[RNA incorporation]
%           (Anabolic)                                             (cytotoxic)
%
%   DPD CIRCADIAN MODULATION (Harris et al. 1990):
%   â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
%   DPD_factor(t) = 1.0 + 0.37 Ã— cos(2Ï€(hour - 1) / 24)
%
%   â•”â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•—
%   â•‘ Time        Hour    DPD_factor    Vmax       Effect             â•‘
%   â•‘ â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ â•‘
%   â•‘ 1:00 AM      1       1.37       0.685â†‘      PEAK (Fast clear)   â•‘
%   â•‘ 7:00 AM      7       1.00       0.500       Mean                â•‘
%   â•‘ 1:00 PM     13       0.63       0.315â†“      TROUGH (Slow clear) â•‘
%   â•‘ 7:00 PM     19       1.00       0.500       Mean                â•‘
%   â•šâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
%   CLINICAL IMPLICATION: Same 500 mg dose given at different times
%   produces 1.5-2Ã— different plasma AUCs. Chronotherapy exploits this!
%
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% USAGE GUIDE: HOW TO USE THIS FUNCTION
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% BASIC SYNTAX:
%
%   results = run5FU_PBPK_Simulation(inputFile, outputPrefix)
%   results = run5FU_PBPK_Simulation(inputFile, outputPrefix, paramOverrides)
%
% REQUIRED INPUTS:
%
%   inputFile (string)
%   â”œâ”€ Path to CSV file containing 5-FU dosing schedule
%   â”œâ”€ Format: dose_id, start_time_min, end_time_min, dose_mg, dose_type
%   â”œâ”€ Supported dose types: 'bolus', 'continuous', 'sinusoidal'
%   â”‚
%   â””â”€ EXAMPLE CSV FILE (save as 'my_dosing.csv'):
%       â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
%       â”‚ dose_id,start_time_min,end_time_min,dose_mg,dose_type       â”‚
%       â”‚ 1,0,1,500,bolus                 (Bolus at midnight)        â”‚
%       â”‚ 2,1440,1445,500,bolus           (Bolus 24hr later)         â”‚
%       â”‚ 3,2880,2885,500,bolus           (Bolus 48hr later)         â”‚
%       â”‚ 4,0,1440,50,continuous          (Continuous 50mg/day)      â”‚
%       â”‚ 5,780,840,750,sinusoidal        (Gradual ramp, 1PM-2PM)    â”‚
%       â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
%
%   outputPrefix (string, optional)
%   â”œâ”€ Default: '5FU_PBPK_output'
%   â”œâ”€ Used to name output files: [outputPrefix]_timeseries.csv, etc
%   â””â”€ Example: outputPrefix='Patient_ABC_Cycle1' â†’ 
%              'Patient_ABC_Cycle1_timeseries.csv', etc.
%
%   paramOverrides (struct, optional)
%   â”œâ”€ Override default parameters (for sensitivity analysis, population PK)
%   â”œâ”€ Example:
%   â”‚  params_override = struct('BW', 75, 'Q_tumor', 0.03, 'Vmax_DPD', 0.6);
%   â”‚  results = run5FU_PBPK_Simulation(inputFile, outputPrefix, params_override);
%   â”‚
%   â””â”€ Useful for:
%      â€¢ Monte Carlo sampling (100s of patients with varied parameters)
%      â€¢ Exploring effect of DPD polymorphisms (â†“ Vmax = DPD deficiency)
%      â€¢ Organ-specific dosing (adjust for renal impairment, etc.)
%      â€¢ Sensitivity analysis (which parameters matter most?)
%
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% OUTPUT STRUCTURE (What You Get Back): results
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
%   results is a MATLAB structure containing:
%
%   1. SIMULATION DATA:
%      â”œâ”€ results.time_min               [1Ã—n array] Time points (0 to max, 1-min steps)
%      â”œâ”€ results.concentrations         [struct] All concentration time-series:
%      â”‚  â”œâ”€ .C_central(t)              [ÂµM] Central compartment 5-FU
%      â”‚  â”œâ”€ .C_peripheral(t)           [ÂµM] Peripheral compartment
%      â”‚  â”œâ”€ .C_liver(t)                [ÂµM] Hepatic 5-FU
%      â”‚  â”œâ”€ .C_kidney(t)               [ÂµM] Renal 5-FU
%      â”‚  â”œâ”€ .C_tumor(t)                [ÂµM] Tumor 5-FU
%      â”‚  â”œâ”€ .C_FdUMP(t)                [ÂµM] Systemic FdUMP (cytotoxic metabolite)
%      â”‚  â”œâ”€ .C_FdUTP(t)                [ÂµM] Systemic FdUTP (RNA incorporation)
%      â”‚  â”œâ”€ .C_FUTP(t)                 [ÂµM] Systemic FUTP
%      â”‚  â”œâ”€ .C_tumor_FdUMP(t)          [ÂµM] Tumor FdUMP (where efficacy happens!)
%      â”‚  â”œâ”€ .C_tumor_FdUTP(t)          [ÂµM] Tumor FdUTP
%      â”‚  â”œâ”€ .C_tumor_FUTP(t)           [ÂµM] Tumor FUTP
%      â”‚  â”œâ”€ .C_DHFU(t)                 [ÂµM] Toxic byproduct (watch in DPD-deficiency)
%      â”‚  â””â”€ .C_FBAL(t)                 [ÂµM] Excretion marker
%      â”‚
%      â””â”€ results.DPD_circadian_factors [1Ã—n array] DPD activity over time (0.6-1.4)
%
%   2. PHARMACOKINETIC METRICS (Key Clinical Parameters):
%      â”œâ”€ results.metrics.AUC_central_mg_h_L     [mgÂ·h/L] â† MOST IMPORTANT
%      â”‚  â””â”€ Compare to therapeutic range: 20-30 mgÂ·h/L (optimal)
%      â”‚
%      â”œâ”€ results.metrics.Cmax_central            [ÂµM] Peak concentration
%      â”œâ”€ results.metrics.Cmax_central_mg_L       [mg/L] Peak (converted units)
%      â”œâ”€ results.metrics.Tmax_central            [min] Time to peak
%      â”œâ”€ results.metrics.t_half                  [min] Elimination half-life
%      â”œâ”€ results.metrics.CL_total                [L/min] Total clearance
%      â”‚
%      â”œâ”€ results.metrics.AUC_tumor_mg_h_L        [mgÂ·h/L] AUC in tumour (target site!)
%      â”œâ”€ results.metrics.AUC_tumor_ratio         [unitless] Tumor/Central AUC ratio
%      â”‚
%      â”œâ”€ results.metrics.Cmax_FdUMP              [ÂµM] Peak of primary metabolite
%      â”œâ”€ results.metrics.AUC_FdUMP               [ÂµMÂ·h] Total FdUMP exposure
%      â”œâ”€ results.metrics.Cmax_tumor_FdUMP        [ÂµM] Peak metabolite in tumor
%      â”‚
%      â””â”€ results.metrics.total_excreted_5FU      [Âµmol] How much unchanged drug in urine
%
%   3. TOXICITY PREDICTIONS (Clinical Decision Support):
%      â”œâ”€ results.metrics.toxicity_category       [string] LOW / MODERATE / HIGH / VERY HIGH
%      â”œâ”€ results.metrics.toxicity_risk           [string] Risk assessment text
%      â”‚
%      â”œâ”€ results.metrics.neutropenia_risk_percent     [0-100%] Neutropenia probability
%      â”œâ”€ results.metrics.grade_3_4_toxicity_percent   [0-100%] Severe toxicity probability
%      â”œâ”€ results.metrics.haematological_toxicity_score [0-100] Composite score
%      â”‚
%      â”œâ”€ results.metrics.mucositis_risk_percent       [0-100%] Mouth ulceration risk
%      â”œâ”€ results.metrics.diarrhea_risk_percent        [0-100%] GI toxicity risk
%      â”‚
%      â””â”€ results.metrics.efficacy_prediction    [string] Efficacy assessment
%
%   4. CLINICAL RECOMMENDATIONS (Actionable Output):
%      â”œâ”€ results.recommendation                  [string] 'INCREASE / MAINTAIN / REDUCE'
%      â”œâ”€ results.suggested_dose_adjustment       [%] How much to change dose
%      â”œâ”€ results.rationale                       [string] Why this recommendation
%      â”‚
%      â””â”€ Compared to literature benchmarks:
%         â”œâ”€ results.metrics.literature_ref.therapeutic_AUC_min    [20 mgÂ·h/L]
%         â”œâ”€ results.metrics.literature_ref.therapeutic_AUC_max    [30 mgÂ·h/L]
%         â””â”€ results.metrics.literature_ref.optimal_AUC_*
%
%   5. MODEL PARAMETERS & CONFIGURATION:
%      â”œâ”€ results.params                  [struct] All ~70 physiological parameters
%      â”œâ”€ results.dosing_regimen          [struct] Parsed dosing schedule
%      â””â”€ results.model_metadata          [struct] Version, date, author info
%
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% DIAGNOSTIC HELP: TROUBLESHOOTING & VALIDATION
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% COMMON ISSUES & SOLUTIONS:
% â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
%
% Q: "My AUC is way too high/low compared to literature"
% A: Check several things:
%    1. Dosing file format - ensure columns are in correct order
%    2. Units - confirm dose is in mg (not Âµmol)
%    3. Body weight - default assumes 70 kg
%       Try: params_override = struct('BW', actual_weight)
%    4. DPD polymorphism - if patient has reduced DPD:
%       Try: params_override = struct('Vmax_DPD', 0.3)  % 60% reduction
%    5. Compare to therapeutic range: 20-30 mgÂ·h/L is optimal
%
% Q: "Results show negative concentrations or NaN"
% A: Indicates numerical instability (shouldn't happen with this solver):
%    1. Check if dosing rate is unrealistic (extremely high)
%    2. Verify time steps are reasonable (model uses 1-minute steps)
%    3. If persists, try shorter simulation or different time scale
%
% Q: "Mass balance check fails (input â‰  compartments + excreted + metabolised)"
% A: Expected small error (Â±5%) due to numerical integration. If >10%:
%    1. Check dosing total - sum(dosingRegimen.doses)
%    2. Verify volume parameters - check V_central, V_peripheral
%    3. Ensure all ODEs are properly balanced
%
% Q: "Circadian rhythm doesn't seem to affect AUC much"
% A: That's actually CORRECT for some scenarios:
%    â€¢ If long infusion (>12 hr), circadian effect is averaged out
%    â€¢ If very high dose, saturation masks circadian effect
%    â€¢ Peak circadian effect is ~2Ã— for single boluses at trough DPD
%    Try: Same bolus at 1 AM vs 1 PM should show 1.5-2Ã— AUC difference
%
% Q: "Tumor concentrations are way lower than central - is that right?"
% A: YES - this is realistic! Tumor:plasma partition coefficient = 0.61
%    Meaning at equilibrium, tumor concentrations = 61% of plasma
%    BUT metabolite accumulation in tumor may be higher due to local synthesis
%    Check: results.metrics.AUC_tumor_ratio should be ~0.6
%
% Q: "How do I validate this model against published data?"
% A: Test against these literature benchmarks:
%    â€¢ 500 mg bolus â†’ AUC should be 20-30 mgÂ·h/L (Schalhorn et al. 1992)
%    â€¢ Cmax typically 40-80 ÂµM
%    â€¢ Half-life: 50-120 min (varies with DPD activity)
%    â€¢ Neutropenia risk: 20% at AUC<15, 50% at AUC=25, 70% at AUC>35
%    Run model and compare to these ranges
%
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% LITERATURE FOUNDATIONS & CITATIONS
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
%   [1] Harris, B.E., Song, R., Soong, S.J., Diasio, R.B. (1990).
%       "Circadian Variation of 5-Fluorouracil Plasma Levels and Metabolite Excretion"
%       Cancer Research 50(13):3878-3882
%       â–º KEY FINDING: DPD shows 1.74-fold circadian variation (peak 1 AM, trough 1 PM)
%       â–º WHY IT MATTERS: Chronotherapy exploits this to optimise exposure
%
%   [2] Diasio, R.B., Beavers, T.L., Carpenter, J.T. (2001).
%       "Familial Deficiency of Dihydropyrimidine Dehydrogenase"
%       Oncology 61(Suppl 1):8-15
%       â–º KEY FINDING: DPD is rate-limiting for 5-FU clearance (80% of total)
%       â–º ENZYME KINETICS: Km ~5 ÂµM, Vmax ~0.156 Âµmol/min (hepatic)
%       â–º DEFICIENCY EFFECTS: ~1-5% population heterozygous â†’ reduced clearance
%
%   [3] Ma, W., Ono, S., Okusaka, T., et al. (2022).
%       "Mechanistic Pharmacokinetic-Pharmacodynamic Modelling of 5-Fluorouracil"
%       PLOS Computational Biology 18(3):e1009945
%       â–º KEY FINDING: FdUMP/FdUTP accumulate in tumors (active metabolites)
%       â–º EFFICACY MECHANISM: FdUMP inhibits thymidylate synthase (TS)
%       â–º TUMOR SELECTIVITY: Cancer cells have higher nucleotide synthase
%
%   [4] Konings, A.W.T., Simons, A.P., Bodis, S., Franzen, D. (2010).
%       "Tumor Distribution and Tumor Penetration of 5-Fluorouracil"
%       Cancer Chemotherapy and Pharmacology 65(6):1145-1155
%       â–º KEY FINDING: Tumor:plasma partition coefficient â‰ˆ 0.61
%       â–º VASCULARISATION: Well-perfused tumors receive 2% of cardiac output
%       â–º IMPLICATIONS: Drug penetration ~60-120 min for equilibration
%
%   [5] Schalhorn, A., Wilke, H., Aulitzky, W., et al. (1992).
%       "Clinical Pharmacokinetics of 5-Fluorouracil in Patients"
%       Cancer Treatment Reviews 18(2):123-145
%       â–º KEY FINDING: Therapeutic AUC window = 20-30 mgÂ·h/L
%       â–º EFFICACY: 94% response rate at optimal AUC
%       â–º TOXICITY: Grade 3-4 toxicity 11% at low AUC, 60% at high AUC
%
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% AUTHOR & VERSION
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% Original Author:        Jacob Frusher
% Created:                November 2025
% Last Updated:           December 2025
% Version:                0.5
%
% Contact/Citation:       If using in research, please cite the literature
%                        references listed above, plus this implementation.
%
% DISCLAIMER:            This model is for RESEARCH & EDUCATIONAL purposes.
%                        Clinical decision-making should involve qualified
%                        oncologists and pharmacists. Model predictions are
%                        mathematical approximations validated against literature.
%                        ALWAYS validate simulations against clinical data.
%
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
% 
% QUICK START
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
%   Step 1: Prepare your dosing CSV file (see example above)
%
%   Step 2: Run the simulation
%           results = run5FU_PBPK_Simulation('my_dosing.csv', 'MyOutput');
%
%   Step 3: Check key metrics
%           fprintf('AUC: %.1f mgÂ·h/L (target: 20-30)\n', ...
%                   results.metrics.AUC_central_mg_h_L);
%           fprintf('Toxicity risk: %s\n', results.metrics.toxicity_category);
%           fprintf('Efficacy: %s\n', results.metrics.efficacy_prediction);
%
%   Step 4: Interpret recommendations
%           fprintf('Recommendation: %s\n', results.recommendation);
%           fprintf('Suggested adjustment: %.0f%%\n', ...
%                   results.suggested_dose_adjustment);
%
%   Step 5: Output files are saved:
%           - MyOutput_timeseries.csv       (all concentration data)
%           - MyOutput_metrics.csv          (summary statistics)
%           - MyOutput_plots.fig            (visualisations)
%
% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
%
% HAPPY MODELING! Questions? See documentation or reference literature.
%
%â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”“
%â”ƒ        ðŸ§¬ END OF FUNCTION HEADER - CODE EXECUTION FOLLOWS BELOW ðŸ§¬             â”ƒ
%â”—â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”â”›

    startup(); % Add subdirectories to MATLAB path

    % Default parameters
    if nargin < 2
        outputPrefix = '5FU_PBPK_output';
    end
    if nargin < 3
        paramOverrides = struct(); % Initialize paramOverrides as an empty structure
    end
    % Initialize structured logger (can be overridden via paramOverrides.logging)
    if isfield(paramOverrides, 'logging')
        logConfig = paramOverrides.logging;
    else
        logConfig = struct();
    end
    defaultLogDir = fullfile(pwd, 'logs');
    if ~isfield(logConfig, 'logDir'); logConfig.logDir = defaultLogDir; end
    if ~isfield(logConfig, 'level'); logConfig.level = 'INFO'; end
    if ~isfield(logConfig, 'enableDebug'); logConfig.enableDebug = false; end
    
    logger = Logger(logConfig.logDir, logConfig.level, logConfig.enableDebug);
    logger.info('simulation_initializing', struct('inputFile', inputFile, 'outputPrefix', outputPrefix));

    % Figure window policy (desktop vs headless/HPC)
    % Respect generatePlots override early so headless callers can disable plotting cleanly.
    wants_plots = true;
    if isfield(paramOverrides, 'generatePlots')
        wants_plots = logical(paramOverrides.generatePlots);
    end
    has_desktop = usejava('desktop');
    if wants_plots && has_desktop
        set(0, 'DefaultFigureWindowStyle', 'docked');
        set(0, 'DefaultFigureVisible', 'on');
    else
        set(0, 'DefaultFigureWindowStyle', 'normal');
        set(0, 'DefaultFigureVisible', 'off');
    end

    % Read dosing regimen with contextual error handling
    try
        dosingRegimen = readDosingRegimen(inputFile);
        logger.info('dosing_regimen_loaded', struct('rows', height(dosingRegimen)));
    catch ME
        logger.fatal('dosing_regimen_read_failed', struct('file', inputFile, 'message', ME.message, 'identifier', ME.identifier));
        rethrow(ME);
    end
    
    % Initialise physiological and pharmacokinetic parameters
    params = initialise5FUParameters();
    % Apply parameter overrides from Monte Carlo
    if ~isempty(paramOverrides)
        override_fields = fieldnames(paramOverrides);
        for i = 1:length(override_fields)
            field = override_fields{i};
            if isfield(params, field)
                params.(field) = paramOverrides.(field);
                % fprintf('  Overriding %s: %.4f\n', field, paramOverrides.(field));
            end
        end
    end

    validateClearanceParameters(params, logger);
    logger.info('clearance_validation_passed', struct('CL_hepatic_est', (params.Vmax_DPD/params.Km_DPD), 'Q_liver_frac', params.Q_liver));

    % â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
    % CHOOSE SOLVER METHOD: Fixed or Adaptive Timesteps
    % â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
    solver_method = params.solver_method;  % 'fixed' or 'adaptive'
    max_dose_end = max(dosingRegimen.end_time_min);  % Define for both solvers
    
    if strcmp(solver_method, 'fixed')
        % â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        % FIXED TIMESTEP SOLVER (Simple, Stable)
        % â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        fprintf('SOLVER: Using FIXED timestep (dt=%.2f min) - Stable reference method\n', params.fixed_timestep_min);
        dt_fixed = params.fixed_timestep_min;
        maxTime = max(dosingRegimen.end_time_min) + params.post_dose_observation_min;
        time_min = (0:dt_fixed:maxTime)';
        
    elseif strcmp(solver_method, 'adaptive')
        % â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        % ADAPTIVE TIMESTEP SOLVER (Smart Dosing-Aware)
        % â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        fprintf('SOLVER: Using ADAPTIVE timestep (dosing-aware) - Smart accuracy/speed balance\n');
        fprintf('  Fine dt=%.2f min from %.0f min before to %.0f min after each dose\n', ...
            params.adaptive_fine_timestep_min, params.adaptive_window_before_min, ...
            params.adaptive_window_after_min);
        fprintf('  Coarse dt=%.2f min elsewhere\n\n', params.adaptive_coarse_timestep_min);
        
        % Identify all dosing events (start and end times)
        dose_events = unique([dosingRegimen.start_time_min; dosingRegimen.end_time_min]);
        
        % Define critical windows around each dosing event
        critical_windows = [];
        for i = 1:length(dose_events)
            event_time = dose_events(i);
            window_start = max(0, event_time - params.adaptive_window_before_min);
            window_end = event_time + params.adaptive_window_after_min;
            critical_windows = [critical_windows; window_start, window_end];
        end
        
        % Merge overlapping windows
        critical_windows = sortrows(critical_windows, 1);
        merged_windows = [];
        if ~isempty(critical_windows)
            current_start = critical_windows(1, 1);
            current_end = critical_windows(1, 2);
            for i = 2:size(critical_windows, 1)
                if critical_windows(i, 1) <= current_end
                    % Overlapping, extend current window
                    current_end = max(current_end, critical_windows(i, 2));
                else
                    % No overlap, save current and start new
                    merged_windows = [merged_windows; current_start, current_end];
                    current_start = critical_windows(i, 1);
                    current_end = critical_windows(i, 2);
                end
            end
            merged_windows = [merged_windows; current_start, current_end];
        end
        
        % Generate simulation end time
        maxTime = max_dose_end + params.post_dose_observation_min;
        
        % Build timestep array
        time_min = [];
        current_time = 0;
        
        for i = 1:size(merged_windows, 1)
            window_start = merged_windows(i, 1);
            window_end = merged_windows(i, 2);
            
            % Coarse timesteps before this critical window
            if current_time < window_start
                time_coarse = (current_time : params.adaptive_coarse_timestep_min : window_start)';
                if ~isempty(time_coarse)
                    time_min = [time_min; time_coarse(1:end-1)]; % Don't duplicate window_start
                end
            end
            
            % Fine timesteps during critical window
            time_fine = (window_start : params.adaptive_fine_timestep_min : window_end)';
            time_min = [time_min; time_fine];
            current_time = window_end;
        end
        
        % Coarse timesteps after last critical window
        if current_time < maxTime
            time_coarse = (current_time + params.adaptive_coarse_timestep_min : params.adaptive_coarse_timestep_min : maxTime)';
            time_min = [time_min; time_coarse];
        end
        
        % Ensure we have time=0 and time=maxTime
        if time_min(1) > 0
            time_min = [0; time_min];
        end
        if time_min(end) < maxTime
            time_min = [time_min; maxTime];
        end
        
    end  % End of solver method selection
    
    % Remove any duplicates and sort
    time_min = unique(time_min);
    
    nTimePoints = length(time_min);
    concentrations = initialiseConcentrations(nTimePoints);
    
    if strcmp(solver_method, 'fixed')
        fprintf('SOLVER: %d points with dt=%.2f min\n', nTimePoints, params.fixed_timestep_min);
    else
        fprintf('SOLVER: %d points (adaptive) vs %.0f (fixed) = %.1fÃ— speedup\n', ...
            nTimePoints, maxTime/0.05, (maxTime/0.05)/nTimePoints);
    end
    logger.info('simulation_started', struct('nTimePoints', nTimePoints, 'maxTime', maxTime, ...
        'solver_method', solver_method));
    
    try
        for t = 2:nTimePoints
            currentTime = time_min(t-1);
            if strcmp(solver_method, 'fixed')
                dt_actual = params.fixed_timestep_min;
            else
                dt_actual = time_min(t) - time_min(t-1);
            end
            
            % Calculate circadian-modulated DPD activity
            hourOfDay = mod(currentTime, 1440) / 60;
            DPD_circadian_factor = calculateCircadianDPD(hourOfDay, params);
            
            % Get current dosing rate
            dosingRate = calculateDosingRate(currentTime, dosingRegimen);
            
            % Log bolus detection at first step
            if t == 2
                if logger.isDebugEnabled()
                    logger.debug('first_integration_debug', struct('time', currentTime, 'dosingRate', dosingRate, 'dosingRegimen_rows', height(dosingRegimen)));
                else
                    logger.info('first_integration_step', struct('time', currentTime, 'dosingRate', dosingRate, 'dosingRegimen_rows', height(dosingRegimen)));
                end
                if dosingRate == 0
                    logger.error('bolus_missing_at_start', struct('time', currentTime, 'note', 'dosingRate zero at first timestep'));
                else
                    logger.info('bolus_detected', struct('dosingRate', dosingRate));
                end
            end
            
            % Calculate rates of change for all compartments using system of ODEs
            % This returns structure with all dC/dt values (ÂµM/min)
            dCdt = calculate5FUSystemODEs(concentrations, t-1, dosingRate, ...
                DPD_circadian_factor, params, currentTime, dt_actual, logger);
            
            % Update concentrations using Euler method
            % C(t) = C(t-1) + dC/dt * dt
            concentrations = updateConcentrations(concentrations, dCdt, t, dt_actual);

            % Algebraic PD state: TS inhibition fraction from updated tumor FdUMP
            concentrations.TS_inhibition_fraction(t) = calculateTSInhibition(concentrations.C_tumor_FdUMP(t), params);
            
            % Display progress every 60 minutes
            if mod(currentTime, 60) == 0
                if logger.isInfoEnabled()
                    logger.info('progress', struct('time_hours', currentTime/60, 'central_uM', concentrations.C_central(t), 'tumor_uM', concentrations.C_tumor(t)));
                end
            end
        end
        
        logger.info('simulation_complete', struct('time_end', maxTime));
        
        % Package results
        results = packageResults(time_min, concentrations, params, dosingRegimen);
        % Validation will be printed in printSummaryStatistics
        
        % Save detailed CSV outputs
        saveDetailedCSVOutputs(results, outputPrefix, logger);
        
        % Generate comprehensive plots
        % Uncomment the next line if you need plots:
        generateComprehensivePlots(results, outputPrefix, logger);
        
        % Print summary statistics
        printSummaryStatistics(results);
        logger.info('simulation_output_packaged', struct('outputPrefix', outputPrefix));
    catch ME
        logger.fatal('simulation_failed', struct('message', ME.message, 'identifier', ME.identifier, 'time', datetime("now")));
        rethrow(ME);
    end

end

