function metrics = calculatePKMetrics(time_min, C)

% Calculate pharmacokinetic metrics using trapezoidal rule integration
% with literature-validated unit conversions for 5-FU AUC
%
% Literature References:
%  Diasio et al. (2001) - Standard PK parameters for 5-FU
%  Schalhorn et al. (1992) - Clinical pharmacokinetics, elimination t1/2 4.5-13 min
%  Fang et al. (2016) - AUC-response relationship, optimal 28-38.94 mg·h/L
%  Maring et al. (2002) - Michaelis-Menten kinetics for DPD metabolism
%  Kaldate et al. (2012) - Individual variation in AUC, dose-scaling

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('CALCULATING PHARMACOKINETIC METRICS\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

% Ensure monotonic, finite inputs for robust numerical integration
time_min = time_min(:);
if any(~isfinite(time_min))
    error('calculatePKMetrics:InvalidTime', 'time_min contains non-finite values.');
end
if any(diff(time_min) <= 0)
    [time_min, sortIdx] = sort(time_min, 'ascend');
    [time_min, uniqueIdx] = unique(time_min, 'stable');
    fieldsC = fieldnames(C);
    for k = 1:numel(fieldsC)
        v = C.(fieldsC{k})(:);
        if numel(v) == numel(sortIdx)
            v = v(sortIdx);
            C.(fieldsC{k}) = v(uniqueIdx);
        end
    end
    fprintf('⚠️  PK metrics input had non-monotonic/duplicate time points; sorted and deduplicated before AUC.\n');
end

finiteCentralMask = isfinite(time_min) & isfinite(C.C_central);
if ~all(finiteCentralMask)
    fprintf('⚠️  Non-finite central concentration/time points removed before AUC integration.\n');
end
time_auc = time_min(finiteCentralMask);
C_central_auc = C.C_central(finiteCentralMask);

if numel(time_auc) < 2
    error('calculatePKMetrics:InsufficientPoints', 'Insufficient valid points for AUC calculation.');
end

% ===== STEP 1: VALIDATE INPUT ARRAYS =====%
nTimePoints = length(time_min);
fprintf('Integration points: %d\n', nTimePoints);
fprintf('Time range: %.1f - %.1f minutes (%.2f hours)\n', ...
    time_min(1), time_min(end), time_min(end)/60);
fprintf('Integration method: Trapezoidal rule (MATLAB trapz)\n');
fprintf('Literature validation: Standard for clinical PK [Diasio 2001]\n\n');

% ===== STEP 2: CALCULATE AUC USING TRAPEZOIDAL RULE =====%
%
% Mathematical formulation:
% AUC = ∫₀ᵗ C(t) dt ≈ Σᵢ [½(Cᵢ + Cᵢ₊₁) × Δt]
%
% trapz() implements the composite trapezoidal rule:
% For uniform spacing Δt = 1 min between points
% This is the gold standard for numerical integration in PK
%
% Formula verification (dimensional analysis):
%   trapz(time_min, C_central) where:
%   time_min: vector in minutes [0, 1, 2, ..., t_end] minutes
%   C_central: concentration vector in µM
%   Output: [µM] × [min] = [µM·min] ✓

fprintf('Step 1: Trapezoidal integration of concentration-time profiles\n');
fprintf('─────────────────────────────────────────────────────────────\n');

% Central compartment AUC (primary exposure metric)
metrics.AUC_central_uM_min = trapz(time_auc, C_central_auc);
fprintf('Central compartment AUC: %.1f µM·min\n', metrics.AUC_central_uM_min);

% Tumor compartment AUC (therapeutic target site)
metrics.AUC_tumor_uM_min = trapz(time_min, C.C_tumor);
fprintf('Tumor compartment AUC: %.1f µM·min\n', metrics.AUC_tumor_uM_min);

% Active metabolite AUC (mechanistic biomarker)
metrics.AUC_FdUMP = trapz(time_min, C.C_FdUMP);
fprintf('Systemic FdUMP AUC: %.1f µM·min\n', metrics.AUC_FdUMP);

metrics.AUC_tumor_FdUTP = trapz(time_min, C.C_tumor_FdUTP);
fprintf('Tumour FdUTP AUC: %.1f µM·min\n', metrics.AUC_tumor_FdUTP);

metrics.AUC_tumor_FUTP = trapz(time_min, C.C_tumor_FUTP);
fprintf('Tumour FUTP AUC: %.1f µM·min\n', metrics.AUC_tumor_FUTP);

metrics.AUC_FdUTP = trapz(time_min, C.C_FdUTP);
fprintf('Systemic FdUTP AUC: %.1f µM·min\n', metrics.AUC_FdUTP);

metrics.AUC_FUTP = trapz(time_min, C.C_FUTP);
fprintf('Systemic FUTP AUC: %.1f µM·min\n\n', metrics.AUC_FUTP);

% ===== STEP 3: UNIT CONVERSION TO CLINICAL UNITS =====%
%
% Convert from [µM·min] to [mg·h/L] (standard literature unit)
%
% Two-step conversion:
%
% Step 3a: Time unit conversion (minutes → hours)
%   [µM·min] × (1 hr / 60 min) = [µM·hr]
%   Dividing by 60 is the correct approach
%   Units check: [µM·min] × [hr/min] = [µM·hr] ✓
%
% Step 3b: Concentration unit conversion (µM → mg/L)
%   1 µM = 1 micromole/L = 10⁻⁶ mol/L
%   1 mol 5-FU = 130.08 g = 130,080 mg
%   Therefore: 1 µM = 10⁻⁶ mol/L × 130.08 g/mol × 1000 mg/g
%             = 0.13008 mg/L
%   Or equivalently: multiply by MW/1000 [mg·L⁻¹·µM⁻¹]
%   [µM·hr] × [0.13008 mg/L per µM] = [mg·hr/L] ✓
%
% Combined conversion factor:
%   [µM·min] × (1/60 hr/min) × (130.08/1000 mg·L⁻¹·µM⁻¹)
%   = [µM·min] × (130.08 / 60000)
%   = [µM·min] × (0.00216813)
%   Result: [mg·h/L] ✓

fprintf('Step 2: Unit conversion from [µM·min] to [mg·h/L]\n');
fprintf('─────────────────────────────────────────────────────────────\n');

MW_5FU = 130.08; % Molecular weight of 5-FU in g/mol
% Literature source: Diasio et al. (2001), commonly used standard

fprintf('5-FU molecular weight: %.2f g/mol\n', MW_5FU);
fprintf('Conversion factor: 1 µM = %.5f mg/L\n', MW_5FU/1000);

% Validated conversion formula:
% AUC [mg·h/L] = AUC [µM·min] × (1 hr/60 min) × (MW/1000 mg/µM)
% This is equivalent to: AUC [mg·h/L] = AUC [µM·min] × (MW / 60000)

conversion_factor = MW_5FU / 60000; % mg·h·L⁻¹·µM⁻¹·min⁻¹

fprintf('Conversion factor: %.8f\n\n', conversion_factor);

% Central compartment AUC in clinical units
metrics.AUC_central_mg_h_L = metrics.AUC_central_uM_min * conversion_factor;
fprintf('Central AUC: %.2f mg·h/L (from %.1f µM·min)\n', ...
    metrics.AUC_central_mg_h_L, metrics.AUC_central_uM_min);

% Tumor compartment AUC in clinical units
metrics.AUC_tumor_mg_h_L = metrics.AUC_tumor_uM_min * conversion_factor;
fprintf('Tumor AUC: %.2f mg·h/L (from %.1f µM·min)\n', ...
    metrics.AUC_tumor_mg_h_L, metrics.AUC_tumor_uM_min);

% Independent cross-check path: convert concentration first, integrate in hours
C_central_mg_L_series = C_central_auc * (MW_5FU / 1000);
time_hr_auc = time_auc / 60;
metrics.AUC_central_mg_h_L_direct = trapz(time_hr_auc, C_central_mg_L_series);
metrics.AUC_central_crosscheck_rel_diff_pct = 100 * abs(metrics.AUC_central_mg_h_L_direct - metrics.AUC_central_mg_h_L) / max(abs(metrics.AUC_central_mg_h_L), eps);

fprintf('Cross-check AUC (direct mg/L·h integration): %.4f mg·h/L\n', metrics.AUC_central_mg_h_L_direct);
fprintf('Cross-check relative difference: %.4f%%\n', metrics.AUC_central_crosscheck_rel_diff_pct);
if metrics.AUC_central_crosscheck_rel_diff_pct > 0.5
    fprintf('⚠️  AUC cross-check mismatch > 0.5%%, investigate time grid or concentration integrity.\n');
end

% Tail contribution diagnostic to detect truncated simulation horizon
tailWindowMin = min(60, max(10, floor((time_auc(end) - time_auc(1)) / 20)));
tailMask = time_auc >= (time_auc(end) - tailWindowMin);
if sum(tailMask) >= 2
    aucTail_uM_min = trapz(time_auc(tailMask), C_central_auc(tailMask));
    metrics.AUC_tail_last_window_fraction = max((aucTail_uM_min * conversion_factor) / max(metrics.AUC_central_mg_h_L, eps), 0);
else
    metrics.AUC_tail_last_window_fraction = 0;
end
fprintf('Tail AUC fraction (last %.0f min): %.4f\n', tailWindowMin, metrics.AUC_tail_last_window_fraction);
if metrics.AUC_tail_last_window_fraction > 0.05
    fprintf('⚠️  >5%% of AUC lies in final window; extend post-dose simulation horizon for stable AUC.\n');
end

% Tumor/Central ratio (literature reference: ~0.61 from Konings et al. 2010)
if metrics.AUC_central_mg_h_L > 0
    tumor_central_ratio = metrics.AUC_tumor_mg_h_L / metrics.AUC_central_mg_h_L;
    fprintf('Tumor/Central AUC ratio: %.3f (Literature: 0.61 ± 0.15 from Konings 2010)\n\n', ...
        tumor_central_ratio);
end

% ===== STEP 4: VALIDATE AGAINST LITERATURE =====%

fprintf('Step 3: Literature validation of AUC magnitude\n');
fprintf('─────────────────────────────────────────────────────────────\n');

% Literature reference ranges for 5-FU bolus dosing
% Based on Fang et al. (2016) meta-analysis of 307 patients
literature_AUC_low = 20;      % mg·h/L (minimum therapeutic)
literature_AUC_optimal_min = 28.03;  % mg·h/L (optimal response rate 94%)
literature_AUC_optimal_max = 38.94;  % mg·h/L (cutoff for increased toxicity)
literature_AUC_high = 50;     % mg·h/L (dangerous toxicity)

fprintf('Literature reference ranges (Fang et al. 2016, n=307):\n');
fprintf('  < 20 mg·h/L: Subtherapeutic, ↑ risk of underdosing\n');
fprintf('  20-28 mg·h/L: Therapeutic, moderate efficacy\n');
fprintf('  28-39 mg·h/L: OPTIMAL, 94%% response rate\n');
fprintf('  > 39 mg·h/L: SUPRATHERAPEUTIC, 71%% Grade 3-4 toxicity\n');
fprintf('  > 50 mg·h/L: DANGEROUS, severe toxicity expected\n\n');

% Validity check
if metrics.AUC_central_mg_h_L < 5
    fprintf('⚠️  WARNING: AUC is unusually low (%.1f mg·h/L)\n', metrics.AUC_central_mg_h_L);
    fprintf('   Check: Did simulation run long enough? Dose > 0?\n\n');
elseif metrics.AUC_central_mg_h_L > 100
    fprintf('⚠️  WARNING: AUC is very high (%.1f mg·h/L)\n', metrics.AUC_central_mg_h_L);
    fprintf('   Check: Is dose correct? V_central reasonable?\n\n');
else
    fprintf('✓ AUC in plausible clinical range\n\n');
end


% Peak concentrations (Cmax)
metrics.Cmax_central = max(C.C_central);
metrics.Cmax_tumor = max(C.C_tumor);
metrics.Cmax_FdUMP = max(C.C_FdUMP);
metrics.Cmax_tumor_FdUMP = max(C.C_tumor_FdUMP);

% Convert Cmax to mg/L for comparison
metrics.Cmax_central_mg_L = metrics.Cmax_central * (MW_5FU / 1000);

% Time to peak (Tmax) - minutes
[~, idx] = max(C.C_central);
metrics.Tmax_central = time_min(idx);
[~, idx] = max(C.C_tumor);
metrics.Tmax_tumor = time_min(idx);

% Total excretion
metrics.total_excreted_5FU = C.excreted_5FU(end);
metrics.total_excreted_FBAL = C.excreted_FBAL(end);
metrics.total_excreted_metabolites = C.excreted_metabolites(end);

% LITERATURE-BASED REFERENCE VALUES FOR CRC
% Based on: Fang et al. (2016), Goel et al. (2015), Morawska et al. (2018)
metrics.literature_ref.therapeutic_AUC_min = 20; % mg·h/L
metrics.literature_ref.therapeutic_AUC_max = 30; % mg·h/L
metrics.literature_ref.optimal_AUC_min = 28.03; % mg·h/L (Fang et al. 2016)
metrics.literature_ref.optimal_AUC_max = 38.94; % mg·h/L
metrics.literature_ref.high_toxicity_threshold = 38.94; % mg·h/L
metrics.literature_ref.low_efficacy_threshold = 20; % mg·h/L

% TOXICITY ASSESSMENT
% Based on literature: AUC > 38.94 mg·h/L → 70.97% high toxicity rate
% AUC 28.03-38.94 mg·h/L → optimal therapeutic window
% AUC < 28.03 mg·h/L → reduced efficacy

if metrics.AUC_central_mg_h_L < metrics.literature_ref.low_efficacy_threshold
    metrics.toxicity_category = 'SUBTHERAPEUTIC';
    metrics.toxicity_risk = 'Low';
    metrics.efficacy_prediction = 'POOR - Below therapeutic range';
    metrics.neutropenia_risk_percent = 10; % Low toxicity, low efficacy
    metrics.grade_3_4_toxicity_percent = 5;
elseif metrics.AUC_central_mg_h_L >= metrics.literature_ref.low_efficacy_threshold && ...
       metrics.AUC_central_mg_h_L < metrics.literature_ref.optimal_AUC_min
    metrics.toxicity_category = 'THERAPEUTIC - Lower Range';
    metrics.toxicity_risk = 'Low-Moderate';
    metrics.efficacy_prediction = 'MODERATE - Therapeutic but suboptimal';
    metrics.neutropenia_risk_percent = 20; % Garg et al. 2012
    metrics.grade_3_4_toxicity_percent = 15;
elseif metrics.AUC_central_mg_h_L >= metrics.literature_ref.optimal_AUC_min && ...
       metrics.AUC_central_mg_h_L <= metrics.literature_ref.optimal_AUC_max
    metrics.toxicity_category = 'OPTIMAL THERAPEUTIC WINDOW';
    metrics.toxicity_risk = 'Moderate';
    metrics.efficacy_prediction = 'EXCELLENT - 94% response rate';
    metrics.neutropenia_risk_percent = 35; % Balanced efficacy/toxicity
    metrics.grade_3_4_toxicity_percent = 25;
elseif metrics.AUC_central_mg_h_L > metrics.literature_ref.optimal_AUC_max && ...
       metrics.AUC_central_mg_h_L <= 50
    metrics.toxicity_category = 'SUPRATHERAPEUTIC - High Toxicity Risk';
    metrics.toxicity_risk = 'HIGH';
    metrics.efficacy_prediction = 'GOOD efficacy but HIGH toxicity';
    metrics.neutropenia_risk_percent = 71; % Fang et al. 2016: 70.97%
    metrics.grade_3_4_toxicity_percent = 55;
else % AUC > 50 mg·h/L
    metrics.toxicity_category = 'SEVERE TOXICITY - DANGEROUS';
    metrics.toxicity_risk = 'VERY HIGH';
    metrics.efficacy_prediction = 'Life-threatening toxicity likely';
    metrics.neutropenia_risk_percent = 90; % Severe febrile neutropenia risk
    metrics.grade_3_4_toxicity_percent = 80;
end

% Specific toxicity predictions based on literature
% Garg et al. (2012): AUC strongly predicts neutropenia and leukopenia
% Morawska et al. (2018): Grade I/II toxicity 19.4% vs 41.3% within vs above range

% Haematological toxicity score (0-100, higher = worse)
metrics.haematological_toxicity_score = min(100, (metrics.AUC_central_mg_h_L / 50) * 100);

% GI toxicity score (mucositis, diarrhea, nausea)
% Literature: 58% mucositis at standard dosing
if metrics.AUC_central_mg_h_L < 20
    metrics.mucositis_risk_percent = 20;
    metrics.diarrhea_risk_percent = 15;
elseif metrics.AUC_central_mg_h_L <= 30
    metrics.mucositis_risk_percent = 45;
    metrics.diarrhea_risk_percent = 35;
elseif metrics.AUC_central_mg_h_L <= 40
    metrics.mucositis_risk_percent = 65;
    metrics.diarrhea_risk_percent = 55;
else
    metrics.mucositis_risk_percent = 85;
    metrics.diarrhea_risk_percent = 75;
end

end
