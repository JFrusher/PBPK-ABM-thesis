function rate = calculateDosingRate(currentTime, dosingRegimen, params, logger)
% calculateDosingRate  Compute instantaneous dosing rate (µmol/min)
% Supports bolus, continuous, sinusoidal dosing. Robust to column names.
if nargin < 3, params = []; end
if nargin < 4, logger = []; end
rate = 0;
MW_5FU = 130.08;

if ~ismember('dosing_type', dosingRegimen.Properties.VariableNames)
    typeCandidates = {'dose_type','doseType','type'};
    for tc = 1:length(typeCandidates)
        if ismember(typeCandidates{tc}, dosingRegimen.Properties.VariableNames)
            dosingRegimen.dosing_type = dosingRegimen.(typeCandidates{tc});
            break;
        end
    end
end

if ~ismember('effective_end_time_min', dosingRegimen.Properties.VariableNames)
    dosingRegimen.effective_end_time_min = dosingRegimen.end_time_min;
end

for i = 1:height(dosingRegimen)
    startTime = dosingRegimen.start_time_min(i);
    endTime_raw = dosingRegimen.end_time_min(i);
    effectiveEnd = dosingRegimen.effective_end_time_min(i);

    % Detect dosing type column
    typeCandidates = {'dosing_type','dose_type','doseType','type'};
    colFound = '';
    for tc = 1:length(typeCandidates)
        if ismember(typeCandidates{tc}, dosingRegimen.Properties.VariableNames)
            colFound = typeCandidates{tc}; break; end
    end
    if isempty(colFound)
        if ~isempty(logger), logger.warnSampled('dosing_type_missing', struct('row', i), 5); end
        dosingType = 'unknown';
    else
        rawVal = dosingRegimen.(colFound)(i);
        if iscell(rawVal)
            dosingType = lower(strtrim(rawVal{1}));
        else
            dosingType = lower(strtrim(string(rawVal)));
        end
    end

    dosingTypeStr = string(dosingType);
    is_bolus = any(contains(dosingTypeStr, "bolus"));
    is_constant = any(contains(dosingTypeStr, ["constant","continuous","infusion","step"]));
    is_sinusoidal = any(contains(dosingTypeStr, ["sinusoidal","chrono","circadian"]));
    is_custom = any(contains(dosingTypeStr, ["custom","function","piecewise"]));

    if ismember('rate_mg_per_min', dosingRegimen.Properties.VariableNames)
        if currentTime >= startTime && currentTime < effectiveEnd
            directRate = dosingRegimen.rate_mg_per_min(i);
            if ~isnan(directRate) && isfinite(directRate) && directRate > 0
                rate = rate + (directRate * 1000) / MW_5FU;
            end
        end
    end

    if is_bolus
            duration = max(endTime_raw - startTime, 0.1);
            effectiveEnd = startTime + duration;
            if currentTime >= startTime && currentTime < effectiveEnd
                % get dose in mg
                dose_mg = NaN;
                if ismember('dose_amount', dosingRegimen.Properties.VariableNames)
                    dose_mg = dosingRegimen.dose_amount(i);
                elseif ismember('dose_mg', dosingRegimen.Properties.VariableNames)
                    dose_mg = dosingRegimen.dose_mg(i);
                elseif ismember('dose_mg_per_m2', dosingRegimen.Properties.VariableNames) && ~isempty(params) && isfield(params,'BSA')
                    dose_mg = dosingRegimen.dose_mg_per_m2(i) * params.BSA;
                end
                if ~isnan(dose_mg) && dose_mg > 0
                    dose_rate_mg_per_min = dose_mg / duration;
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
                duration = max(effectiveEnd - startTime, 0);
                if ~isnan(dose_mg) && dose_mg > 0 && duration > 0
                    rate_mg_per_min = dose_mg / duration;
                end
            end
            if isnan(rate_mg_per_min) && ismember('mean_rate', dosingRegimen.Properties.VariableNames)
                mr = dosingRegimen.mean_rate(i);
                if isfinite(mr)
                    rate_mg_per_min = mr;
                end
            end
            if ~isnan(rate_mg_per_min)
                rate = rate + (rate_mg_per_min * 1000) / MW_5FU;
            end
        end
    end

    if is_sinusoidal
        if currentTime >= startTime && currentTime < effectiveEnd
            mean_rate = dosingRegimen.mean_rate(i);
            amp = dosingRegimen.amplitude(i);
            freq = dosingRegimen.frequency_per_min(i);
            if isnan(mean_rate) && ismember('infusion_rate', dosingRegimen.Properties.VariableNames)
                mean_rate = dosingRegimen.infusion_rate(i);
            end
            if isnan(amp), amp = 0; end
            if isnan(freq) || freq <= 0
                period = max(effectiveEnd - startTime, 1.0);
                freq = 1 / period;
            end
            instantaneous = mean_rate + amp * sin(2*pi*freq*currentTime);
            instantaneous = max(instantaneous, 0);
            if isfinite(instantaneous)
                rate = rate + (instantaneous * 1000) / MW_5FU;
            end
        end
    end

    if is_custom
        if currentTime >= startTime && currentTime < effectiveEnd
            customRate = NaN;
            if ismember('custom_function', dosingRegimen.Properties.VariableNames)
                expr = strtrim(string(dosingRegimen.custom_function(i)));
                if strlength(expr) > 0
                    try
                        fn = str2func(['@(t)' char(expr)]);
                        customRate = fn(currentTime);
                    catch
                        customRate = NaN;
                    end
                end
            end
            if isnan(customRate) && ismember('infusion_rate', dosingRegimen.Properties.VariableNames)
                customRate = dosingRegimen.infusion_rate(i);
            end
            if ~isnan(customRate) && isfinite(customRate)
                customRate = max(customRate, 0);
                rate = rate + (customRate * 1000) / MW_5FU;
            end
        end
    end
end
end