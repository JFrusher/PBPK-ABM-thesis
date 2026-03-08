function dosing = readDosingRegimen(csvFile)
% readDosingRegimen  Read & validate a dosing CSV into a standardized table
% Required columns: start_time_min, end_time_min, dosing_type (or dose_type)

if ~isfile(csvFile)
    error('readDosingRegimen:NotFound','File not found: %s', csvFile);
end
T = readtable(csvFile);
required = {'start_time_min','end_time_min'};
if ~all(ismember(required, T.Properties.VariableNames))
    error('MATLAB:unassignedOutputs','Missing required columns in dosing CSV');
end
% Dosing type column: accept several names
typeCandidates = {'dosing_type','dose_type','doseType','type'};
foundType = '';
for i=1:length(typeCandidates)
    if ismember(typeCandidates{i}, T.Properties.VariableNames)
        foundType = typeCandidates{i}; break;
    end
end
if isempty(foundType)
    error('MATLAB:unassignedOutputs','Missing dosing type column in dosing CSV');
end

if ~strcmp(foundType, 'dosing_type')
    T.dosing_type = T.(foundType);
end

if ismember('dose_mg', T.Properties.VariableNames) && ~ismember('dose_amount', T.Properties.VariableNames)
    T.dose_amount = T.dose_mg;
end

optionalNumeric = {'dose_amount','infusion_rate','mean_rate','amplitude','frequency_per_min','rate_mg_per_min'};
for i = 1:numel(optionalNumeric)
    c = optionalNumeric{i};
    if ~ismember(c, T.Properties.VariableNames)
        T.(c) = nan(height(T), 1);
    end
end
if ~ismember('custom_function', T.Properties.VariableNames)
    T.custom_function = repmat({''}, height(T), 1);
end

numericCols = [{'start_time_min','end_time_min'}, optionalNumeric];
for i = 1:numel(numericCols)
    c = numericCols{i};
    col = T.(c);
    if iscell(col) || isstring(col) || iscategorical(col)
        col = str2double(string(col));
    end
    T.(c) = col;
end

T.dosing_type = lower(strtrim(string(T.dosing_type)));
[~, ord] = sort(T.start_time_min, 'ascend');
T = T(ord,:);

n = height(T);
effective_end = T.end_time_min;
for i = 1:n
    s = T.start_time_min(i);
    e = T.end_time_min(i);
    if ~isfinite(s)
        error('readDosingRegimen:InvalidTime', 'Non-finite start_time_min at row %d.', i);
    end
    if ~isfinite(e)
        e = s;
    end
    if e < s
        tmp = s; s = e; e = tmp;
        T.start_time_min(i) = s;
        T.end_time_min(i) = e;
    end

    type_i = string(T.dosing_type(i));
    is_bolus = any(contains(type_i, "bolus"));
    is_nonbolus = any(contains(type_i, ["continuous","constant","infusion","sinusoidal","chrono","custom","step"]));
    if e <= s && is_nonbolus && ~is_bolus
        nextStart = NaN;
        if i < n
            futureStarts = T.start_time_min((i+1):end);
            futureStarts = futureStarts(futureStarts > s);
            if ~isempty(futureStarts)
                nextStart = min(futureStarts);
            end
        end
        if isfinite(nextStart)
            e = nextStart;
        else
            e = s + 1.0;
        end
    end
    effective_end(i) = max(e, s);
end

T.effective_end_time_min = effective_end;
T.effective_duration_min = max(T.effective_end_time_min - T.start_time_min, 0);

dosing = T;
end