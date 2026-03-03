function rate = calculateDosingRate(currentTime, dosingRegimen, params, logger)
% calculateDosingRate  Compute instantaneous dosing rate (µmol/min)
% Supports bolus, continuous, sinusoidal dosing. Robust to column names.
if nargin < 3, params = []; end
if nargin < 4, logger = []; end
rate = 0;
MW_5FU = 130.08;

for i = 1:height(dosingRegimen)
    startTime = dosingRegimen.start_time_min(i);
    endTime_raw = dosingRegimen.end_time_min(i);
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

    switch dosingType
        case 'bolus'
            duration = max(endTime_raw - startTime, 0.1);
            effectiveEnd = startTime + duration;
            if currentTime >= startTime && currentTime < effectiveEnd
                % get dose in mg
                dose_mg = NaN;
                if ismember('dose_mg', dosingRegimen.Properties.VariableNames)
                    dose_mg = dosingRegimen.dose_mg(i);
                elseif ismember('dose_amount', dosingRegimen.Properties.VariableNames)
                    dose_mg = dosingRegimen.dose_amount(i);
                elseif ismember('dose_mg_per_m2', dosingRegimen.Properties.VariableNames) && ~isempty(params) && isfield(params,'BSA')
                    dose_mg = dosingRegimen.dose_mg_per_m2(i) * params.BSA;
                end
                if ~isnan(dose_mg) && dose_mg > 0
                    dose_rate_mg_per_min = dose_mg / duration;
                    rate = rate + (dose_rate_mg_per_min * 1000) / MW_5FU;
                end
            end
        case {'continuous','constant'}
            if currentTime >= startTime && currentTime <= endTime_raw
                rate_mg_per_min = dosingRegimen.infusion_rate(i);
                if ~isnan(rate_mg_per_min)
                    rate = rate + (rate_mg_per_min * 1000) / MW_5FU;
                end
            end
        case 'sinusoidal'
            if currentTime >= startTime && currentTime <= endTime_raw
                mean_rate = dosingRegimen.mean_rate(i);
                amp = dosingRegimen.amplitude(i);
                freq = dosingRegimen.frequency_per_min(i);
                instantaneous = mean_rate + amp * sin(2*pi*freq*currentTime);
                rate = rate + (instantaneous * 1000) / MW_5FU;
            end
        otherwise
            % unknown -> ignore
    end
end
end