function rec = determineClinicalRecommendation(metrics)
% Return recommendation, suggested_dose_adjustment (%), and rationale string
% based on simulated AUC vs literature thresholds.
    AUC = metrics.AUC_central_mg_h_L;
    ref = metrics.literature_ref;

    if AUC < ref.low_efficacy_threshold
        target_AUC = 25;
        dose_adj = ((target_AUC / AUC) - 1) * 100;
        rec.recommendation = 'INCREASE';
        rec.suggested_dose_adjustment = dose_adj;
        rec.rationale = sprintf('AUC %.1f mg·h/L below therapeutic range; +%.0f%% targets 25 mg·h/L', AUC, dose_adj);

    elseif AUC >= ref.optimal_AUC_min && AUC <= ref.optimal_AUC_max
        rec.recommendation = 'MAINTAIN';
        rec.suggested_dose_adjustment = 0;
        rec.rationale = sprintf('AUC %.1f mg·h/L in optimal window (%.0f–%.0f mg·h/L)', AUC, ref.optimal_AUC_min, ref.optimal_AUC_max);

    elseif AUC > ref.optimal_AUC_max
        target_AUC = 30;
        dose_adj = ((target_AUC / AUC) - 1) * 100;
        rec.recommendation = 'REDUCE';
        rec.suggested_dose_adjustment = dose_adj;
        rec.rationale = sprintf('AUC %.1f mg·h/L exceeds safe range; %.0f%% reduction targets 30 mg·h/L', AUC, abs(dose_adj));

    else
        rec.recommendation = 'MONITOR';
        rec.suggested_dose_adjustment = 0;
        rec.rationale = sprintf('AUC %.1f mg·h/L; continue PK-guided dose adjustment', AUC);
    end
end
