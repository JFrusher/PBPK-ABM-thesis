function results = packageResults(time_min, concentrations, params, dosingRegimen)
% Package all results into organized structure

results.time_min = time_min;
results.time_hr = time_min / 60;

% Concentrations
results.concentrations = concentrations;

% Parameters
results.params = params;

% Dosing regimen
results.dosingRegimen = dosingRegimen;

% Calculate AUC and other PK metrics
results.metrics = calculatePKMetrics(time_min, concentrations);
results.mass_balance = calculateMassBalance(time_min, concentrations, params);

rec = determineClinicalRecommendation(results.metrics);
results.recommendation = rec.recommendation;
results.suggested_dose_adjustment = rec.suggested_dose_adjustment;
results.rationale = rec.rationale;
end
