function rate = calculateDosingRate(currentTime, dosingRegimen, allowCustomEval)
% DOSING RATE CALCULATION FUNCTION
% ================================================================================
%
% FUNCTION: calculateDosingRate(currentTime, dosingRegimen)
%
% PURPOSE:  Calculate the instantaneous 5-FU infusion rate (µmol/min) at any
%           given time point, accounting for multiple concurrent dosing regimens
%           (bolus, continuous infusion, sinusoidal patterns). Handles unit
%           conversion from mg/min to µmol/min and robust string parsing of
%           dosing type identifiers.
%
% INPUT:    currentTime (scalar, minutes from t=0)
%           dosingRegimen (table with columns: start_time_min, end_time_min,
%           dosing_type, dose_mg or dose_amount, infusion_rate, mean_rate,
%           amplitude, frequency_per_min)
%
% OUTPUT:   rate (scalar, µmol/min) - Total infusion rate from all active
%           dosing regimens at currentTime. Returns 0 if no active dosing.
%
% ================================================================================
% DOSING REGIMEN ARCHITECTURE
% ================================================================================
%
% This function supports three distinct dosing patterns, each with different
% temporal dynamics and clinical applications:
%
% 1. BOLUS (Intravenous Push)
%    • Single rapid injection over short duration (typically <5 minutes)
%    • Characteristic: High peak concentration, rapid decline
%    • Example: 400 mg over 2 minutes → rate = 400/2 = 200 mg/min
%    • Clinical use: Standard 5-FU bolus regimens (e.g., Mayo Clinic protocol)
%
% 2. CONTINUOUS INFUSION
%    • Steady-state drug delivery over extended period (hours to days)
%    • Characteristic: Flat concentration profile during infusion
%    • Example: 15 mg/min constant rate for 8 hours
%    • Clinical use: De Gramont regimen, protracted 5-FU schedules
%
% 3. SINUSOIDAL INFUSION
%    • Time-varying rate following periodic function
%    • Characteristic: Oscillating concentration, exploiting circadian DPD variation
%    • Model: rate(t) = mean_rate + amplitude × sin(2πf × t)
%    • Clinical use: Chronotherapy protocols leveraging circadian biology
%
% The function iterates through ALL regimen rows and SUMS rates from concurrent
% regimens, enabling complex polypharmacy scenarios.
%
% ================================================================================
% UNIT CONVERSION AND MOLECULAR WEIGHT
% ================================================================================
%
% 5-FU molecular weight: MW = 130.08 g/mol
%
% Conversion formula: rate_µmol_min = (rate_mg_min × 1000) / MW
%
% Example:
%   • Bolus: 400 mg over 2 minutes
%   • rate_mg_min = 400 mg / 2 min = 200 mg/min
%   • rate_µmol_min = (200 × 1000) / 130.08 = 1538 µmol/min
%
% Inverse conversion (for reference):
%   • rate_mg_min = (rate_µmol_min × MW) / 1000
%
% ================================================================================
% BOLUS HANDLING (CRITICAL FIX)
% ================================================================================
%
% PROBLEM: Zero-length boluses (start_time = end_time) cause division by zero
% or undefined behavior in concentration calculations.
%
% SOLUTION: Enforce minimum duration of 0.1 minutes (6 seconds)
%
% Code logic:
%   duration = max(endTime_raw - startTime, 0.1)
%   effectiveEndTime = startTime + duration
%
% If start=60, end=60:
%   → duration = max(0, 0.1) = 0.1 min
%   → effectiveEndTime = 60.1 min
%   → Bolus rate = dose_mg / 0.1 (high but finite)
%
% This prevents singularities while maintaining physical realism for rapid boluses.
%
% ================================================================================
% COLUMN NAME ROBUSTNESS
% ================================================================================
%
% The function handles two equivalent column name conventions for dose amount:
%
%   if column 'dose_mg' exists:
%     → Use dosingRegimen.dose_mg(i)
%   else if column 'dose_amount' exists:
%     → Use dosingRegimen.dose_amount(i)
%   else:
%     → Set dose_mg = NaN (invalid, skipped)
%
% This flexibility accommodates different table generation standards or data
% sources. The validation check (if ~isnan(dose_mg) && dose_mg > 0) ensures
% only positive, defined doses are processed.
%
% ================================================================================
% STRING PARSING FOR DOSING TYPE
% ================================================================================
%
% Robust handling of dosing_type column (may be cell array or string array):
%
%   if iscell(dosingRegimen.dosing_type):
%     → Extract cell content: dosingRegimen.dosing_type{i}
%   else:
%     → Convert to string: string(dosingRegimen.dosing_type(i))
%
%   Apply: lower(strtrim(...))
%     • strtrim: removes leading/trailing whitespace
%     • lower: converts to lowercase for case-insensitive matching
%
% Allows input variants: 'Bolus', 'BOLUS', ' bolus ', all map to 'bolus'
%
% ================================================================================
% REGIMEN MATCHING LOGIC
% ================================================================================
%
% For each dosing type, the function checks if currentTime falls within the
% active window:
%
% BOLUS:
%   Condition: currentTime >= startTime AND currentTime < effectiveEndTime
%   Action: Add bolus rate = dose_mg / duration
%   (Note: Use < (strict less than) to prevent double-counting at boundaries)
%
% CONTINUOUS/CONSTANT:
%   Condition: currentTime >= startTime AND currentTime <= endTime_raw
%   Action: Add constant infusion_rate
%   (Note: Use <= to include endpoint for continuous regimens)
%
% SINUSOIDAL:
%   Condition: currentTime >= startTime AND currentTime <= endTime_raw
%   Action: Evaluate rate(t) = mean_rate + amplitude × sin(2πft)
%   Safety: max(rate_mg_per_min, 0) enforces non-negative dosing
%   (Prevents negative rates from sine function oscillation)
%
% ================================================================================
% SINUSOIDAL CHRONOTHERAPY MODEL
% ================================================================================
%
% Circadian dosing exploits the ~1.5-1.7 fold variation in DPD activity across
% 24 hours [Harris et al. 1990]. Delivering higher doses when DPD activity is
% high (daytime) minimizes exposure when activity is low (low DPD = high exposure
% = toxicity risk).
%
% Mathematical form:
%   rate(t) = mean_rate + amplitude × sin(2πf × t)
%
% Parameters (example):
%   • mean_rate = 10 mg/min
%   • amplitude = 3 mg/min
%   • frequency = 1/(1440) min⁻¹ (one cycle per 24 hours = 1440 min)
%
%   → At t=0 (midnight): rate = 10 + 3×sin(0) = 10 mg/min
%   → At t=360 (6 AM): rate = 10 + 3×sin(π/2) = 13 mg/min (peak)
%   → At t=720 (noon): rate = 10 + 3×sin(π) = 10 mg/min
%   → At t=1080 (6 PM): rate = 10 + 3×sin(3π/2) = 7 mg/min (trough)
%
% max(rate, 0) ensures rate never goes negative (physical constraint).
%
% ================================================================================
% MULTIPLE CONCURRENT REGIMENS
% ================================================================================
%
% The loop structure iterates through ALL rows of dosingRegimen table and SUMS
% the active rates:
%
%   rate_total = 0
%   for each row i:
%     if currentTime in active window for regimen i:
%       rate_total += rate_i
%
% This enables complex scenarios:
%   • Two boluses at same time (additive effect)
%   • Continuous infusion + bolus overlay (typical in clinical protocols)
%   • Multiple sinusoidal regimens with different phases (experimental design)
%
% Example (De Gramont protocol variant):
%   Row 1: Bolus 400 mg at t=0-2 min → rate₁ = 200 mg/min (0-2 min only)
%   Row 2: Continuous 15 mg/min at t=0-120 min → rate₂ = 15 mg/min (throughout)
%
%   At t=1 min: rate_total = 200 + 15 = 215 mg/min
%   At t=5 min: rate_total = 0 + 15 = 15 mg/min
%   At t=150 min: rate_total = 0 + 0 = 0 mg/min
%
% ================================================================================
% CLINICAL APPLICATION EXAMPLES
% ================================================================================
%
% BOLUS REGIMEN (Mayo Clinic):
%   start: 0 min, end: 3 min, type: 'bolus', dose: 400 mg
%   → Produces spike at t=0-3, rate ≈ 133 mg/min = 1023 µmol/min
%   → Rapid rise and fall in plasma concentration
%
% CONTINUOUS INFUSION (5-day schedule):
%   start: 0 min, end: 7200 min (5 days), type: 'continuous',
%   infusion_rate: 15 mg/min
%   → Constant exposure over extended period
%   → Steady-state reached within 2-3 half-lives (~30 min)
%
% CHRONOTHERAPY (Lévi protocol):
%   start: 0, end: 1440 (24 hours), type: 'sinusoidal',
%   mean_rate: 12 mg/min, amplitude: 4 mg/min, frequency: 1/1440 min⁻¹
%   → Variable infusion matched to circadian DPD activity
%   → Higher rates during high-DPD period, lower during low-DPD period
%   → May reduce toxicity while maintaining efficacy
%
% ================================================================================
% RETURN VALUE AND DOWNSTREAM USE
% ================================================================================
%
% The returned rate (µmol/min) is used as the dosing term in the ODE system:
%
%   d[5-FU]/dt = rate_in(t) - rate_out_hepatic(t) - rate_out_renal(t) + ...
%                │────────────┬────────────────────────────────────────│
%                Input from       Output via metabolism & clearance
%                calculateDosingRate()
%
% At each ODE time step, this function is called with the current time to obtain
% the instantaneous infusion rate, ensuring realistic dosing dynamics throughout
% the simulation.
%
% ================================================================================
% REFERENCES
% ================================================================================
%
% Harris, B. E., Song, R., & Soong, S. J., et al. (1990). "Circadian variation
% in plasma 5-fluorouracil levels during prolonged infusion." Journal of Clinical
% Oncology, 8(7), 1192-1197.
%
% Lévi, F., Zidani, R., & Misset, J. L. (1997). "Randomised multicentre trial of
% chronotherapy with oxaliplatin, fluorouracil, and folinic acid in metastatic
% colorectal cancer." The Lancet, 350(9075), 681-686.
%
% ================================================================================
    rate = 0;
    MW_5FU = 130.08;
    if nargin < 3 || isempty(allowCustomEval), allowCustomEval = false; end

    for i = 1:height(dosingRegimen)
        startTime = dosingRegimen.start_time_min(i);
        endTime_raw = dosingRegimen.end_time_min(i);
        effectiveEnd = endTime_raw;
        if ismember('effective_end_time_min', dosingRegimen.Properties.VariableNames)
            effectiveEnd = dosingRegimen.effective_end_time_min(i);
        end

        % Robust string handling
        if iscell(dosingRegimen.dosing_type)
            dosingType = lower(strtrim(string(dosingRegimen.dosing_type{i})));
        else
            dosingType = lower(strtrim(string(dosingRegimen.dosing_type(i))));
        end

        dosingTypeStr = string(dosingType);
        is_bolus = any(contains(dosingTypeStr, "bolus"));
        is_constant = any(contains(dosingTypeStr, ["constant", "continuous", "infusion", "step"]));
        is_sinusoidal = any(contains(dosingTypeStr, ["sinusoidal", "chrono", "circadian"]));
        is_custom = any(contains(dosingTypeStr, ["custom", "function", "piecewise"]));

        % Generic direct rate column support
        if ismember('rate_mg_per_min', dosingRegimen.Properties.VariableNames)
            if currentTime >= startTime && currentTime < effectiveEnd
                directRate = dosingRegimen.rate_mg_per_min(i);
                if ~isnan(directRate) && isfinite(directRate) && directRate > 0
                    rate = rate + (directRate * 1000) / MW_5FU;
                end
            end
        end

        if is_bolus
            % FIX: Calculate duration FIRST to handle 0-length boluses
            % If start=60, end=60 -> duration=0.1, effectiveEnd=60.1
            duration = max(effectiveEnd - startTime, 0.1);
            effectiveEndTime = startTime + duration;

            % Check against EFFECTIVE window
            if currentTime >= startTime && currentTime < effectiveEndTime
                % Get dose in mg (with column detection)
                dose_mg = NaN;
                if ismember('dose_amount', dosingRegimen.Properties.VariableNames)
                    dose_mg = dosingRegimen.dose_amount(i);
                elseif ismember('dose_mg', dosingRegimen.Properties.VariableNames)
                    dose_mg = dosingRegimen.dose_mg(i);
                end

                % Validate and Add
                if ~isnan(dose_mg) && dose_mg > 0
                    % Rate = Amount / Duration
                    dose_rate_mg_per_min = dose_mg / duration;
                    % Convert to µmol/min
                    rate = rate + (dose_rate_mg_per_min * 1000) / MW_5FU;
                end
            end
        end

        if is_constant
            if currentTime >= startTime && currentTime < effectiveEnd
                rate_mg_per_min = dosingRegimen.infusion_rate(i);
                if isnan(rate_mg_per_min)
                    dose_mg = NaN;
                    if ismember('dose_amount', dosingRegimen.Properties.VariableNames)
                        dose_mg = dosingRegimen.dose_amount(i);
                    elseif ismember('dose_mg', dosingRegimen.Properties.VariableNames)
                        dose_mg = dosingRegimen.dose_mg(i);
                    end
                    duration = effectiveEnd - startTime;
                    if ~isnan(dose_mg) && dose_mg > 0 && duration > 0
                        rate_mg_per_min = dose_mg / duration;
                    end
                end
                if isnan(rate_mg_per_min) && ismember('mean_rate', dosingRegimen.Properties.VariableNames)
                    mr = dosingRegimen.mean_rate(i);
                    if isfinite(mr) && ~isnan(mr)
                        rate_mg_per_min = mr;
                    end
                end
                if ~isnan(rate_mg_per_min)
                    rate = rate + (rate_mg_per_min / MW_5FU) * 1000;
                end
            end
        end

        if is_sinusoidal
            if currentTime >= startTime && currentTime < effectiveEnd
                mean_rate = dosingRegimen.mean_rate(i);
                amplitude = dosingRegimen.amplitude(i);
                frequency = dosingRegimen.frequency_per_min(i);

                if isnan(mean_rate) && ismember('infusion_rate', dosingRegimen.Properties.VariableNames)
                    mean_rate = dosingRegimen.infusion_rate(i);
                end
                if isnan(amplitude)
                    amplitude = 0;
                end
                if isnan(frequency) || frequency <= 0
                    period = max(effectiveEnd - startTime, 1.0);
                    frequency = 1 / period;
                end

                rate_mg_per_min = mean_rate + amplitude * sin(2 * pi * frequency * currentTime);
                rate_mg_per_min = max(rate_mg_per_min, 0);

                rate = rate + (rate_mg_per_min / MW_5FU) * 1000;
            end
        end

        if is_custom
            if currentTime >= startTime && currentTime < effectiveEnd
                customRate = NaN;
                if ismember('custom_function', dosingRegimen.Properties.VariableNames)
                    exprRaw = dosingRegimen.custom_function(i);
                    expr = strtrim(string(exprRaw));
                    if strlength(expr) > 0 && allowCustomEval
                        try
                            % SECURITY: executes arbitrary MATLAB from the input CSV's
                            % custom_function column. Gated behind allowCustomEval
                            % (params.allow_custom_function_eval, default false). Only
                            % enable for self-authored, trusted dosing files.
                            fn = str2func(['@(t)' char(expr)]);
                            customRate = fn(currentTime);
                        catch
                            customRate = NaN;
                        end
                    end
                end

                if isnan(customRate)
                    if ismember('infusion_rate', dosingRegimen.Properties.VariableNames)
                        customRate = dosingRegimen.infusion_rate(i);
                    end
                end

                if ~isnan(customRate) && isfinite(customRate)
                    customRate = max(customRate, 0);
                    rate = rate + (customRate / MW_5FU) * 1000;
                end
            end
        end
    end
end
