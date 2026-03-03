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
% Normalize dosing type to cellstr
T.(foundType) = cellstr(string(T.(foundType)));
% Standardize output columns to expected names
% If dose_amount exists, prefer it for compatibility
if ismember('dose_amount', T.Properties.VariableNames) && ~ismember('dose_mg', T.Properties.VariableNames)
    T.dose_mg = T.dose_amount;
end
% Return table
dosing = T;
end