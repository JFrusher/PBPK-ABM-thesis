function report_filename = generate_comprehensive_stats_report(output_dir, Cmax, AUC, ...
    n_samples, analysis_date_str, input_csv, param_sensitivities, param_names)
% Exports a detailed statistics report suitable for a paper appendix.

    if nargin < 8
        param_names = [];
    end
    if nargin < 7
        param_sensitivities = [];
    end

    report_filename = fullfile(output_dir, 'MC_Statistical_Analysis_Report.txt');
    fid = fopen(report_filename, 'w');
    if fid == -1
        error('generate_comprehensive_stats_report:FileError', ...
            'Could not open %s for writing', report_filename);
    end

    write_header = @(title) fprintf(fid, '\n%s\n%s\n', title, repmat('=', 1, length(title)));
    write_subheader = @(title) fprintf(fid, '\n%s\n%s\n', title, repmat('-', 1, length(title)));

    fprintf(fid, '5-FU PHARMACOKINETICS: MONTE CARLO SENSITIVITY ANALYSIS\n');
    fprintf(fid, 'COMPREHENSIVE STATISTICAL REPORT FOR PUBLICATION\n\n');

    fprintf(fid, 'ANALYSIS METADATA:\n');
    fprintf(fid, '  Input Dosing Scenario: %s\n', input_csv);
    fprintf(fid, '  Monte Carlo Samples:   %d simulations\n', n_samples);
    fprintf(fid, '  Analysis Date:         %s\n', analysis_date_str);
    fprintf(fid, '  Report Generated:      %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    write_header('CMAX (PEAK CONCENTRATION) ANALYSIS');
    fprintf(fid, 'Units: uM (micromolar)\n');
    write_subheader('Core Statistics');

    mean_cmax = mean(Cmax);
    median_cmax = median(Cmax);
    std_cmax = std(Cmax);
    cv_cmax = std_cmax / mean_cmax * 100;

    fprintf(fid, '\nCentral Tendency:\n');
    fprintf(fid, '  Mean:               %.3f +/- %.3f uM\n', mean_cmax, std_cmax);
    fprintf(fid, '  Median:             %.3f uM\n', median_cmax);
    fprintf(fid, '  Coefficient of Variation: %.1f%%\n\n', cv_cmax);

    fprintf(fid, 'Distribution Span:\n');
    fprintf(fid, '  Minimum:            %.3f uM\n', min(Cmax));
    fprintf(fid, '  Maximum:            %.3f uM\n', max(Cmax));
    fprintf(fid, '  Range:              %.3f uM\n', max(Cmax) - min(Cmax));
    fprintf(fid, '  Median Absolute Deviation: %.3f uM\n\n', mad(Cmax, 1));

    p5_cmax = prctile(Cmax, 5);
    p10_cmax = prctile(Cmax, 10);
    p25_cmax = prctile(Cmax, 25);
    p50_cmax = prctile(Cmax, 50);
    p75_cmax = prctile(Cmax, 75);
    p90_cmax = prctile(Cmax, 90);
    p95_cmax = prctile(Cmax, 95);
    iqr_cmax = p75_cmax - p25_cmax;

    write_subheader('Percentile Distribution');
    fprintf(fid, '  5th percentile:     %.3f uM (5%% of samples below)\n', p5_cmax);
    fprintf(fid, '  10th percentile:    %.3f uM\n', p10_cmax);
    fprintf(fid, '  25th percentile:    %.3f uM (Q1)\n', p25_cmax);
    fprintf(fid, '  50th percentile:    %.3f uM (Q2/Median)\n', p50_cmax);
    fprintf(fid, '  75th percentile:    %.3f uM (Q3)\n', p75_cmax);
    fprintf(fid, '  90th percentile:    %.3f uM\n', p90_cmax);
    fprintf(fid, '  95th percentile:    %.3f uM (95%% of samples below)\n\n', p95_cmax);

    fprintf(fid, 'Interquartile Range (IQR = Q3-Q1):\n');
    fprintf(fid, '  IQR:                %.3f uM\n', iqr_cmax);
    fprintf(fid, '  Interpretation:     50%% of values fall within %.3f-%.3f uM\n\n', p25_cmax, p75_cmax);

    write_subheader('Literature Benchmark Comparison');
    try
        bench_saif = get5FU_literature_benchmarks('Cmax', 'saif2013');
        fprintf(fid, '\nSaif et al. (2013) Clinical Cutoffs:\n');
        fprintf(fid, '  Safe range (<%.0f uM):           %.1f%% of samples\n', ...
            bench_saif.safe_max, 100 * sum(Cmax < bench_saif.safe_max) / length(Cmax));
        fprintf(fid, '  Moderate toxicity risk (%.0f-%.0f uM): %.1f%% of samples\n', ...
            bench_saif.safe_max, bench_saif.moderate_threshold, ...
            100 * sum(Cmax >= bench_saif.safe_max & Cmax < bench_saif.moderate_threshold) / length(Cmax));
        fprintf(fid, '  High toxicity risk (>%.0f uM):  %.1f%% of samples\n', ...
            bench_saif.high_threshold, ...
            100 * sum(Cmax >= bench_saif.high_threshold) / length(Cmax));
    catch
        fprintf(fid, 'Literature benchmarks unavailable\n');
    end

    write_header('AUC (SYSTEMIC EXPOSURE) ANALYSIS');
    fprintf(fid, 'Units: mg*h/L (milligrams*hours per liter)\n');
    write_subheader('Core Statistics');

    mean_auc = mean(AUC);
    median_auc = median(AUC);
    std_auc = std(AUC);
    cv_auc = std_auc / mean_auc * 100;

    fprintf(fid, '\nCentral Tendency:\n');
    fprintf(fid, '  Mean:               %.2f +/- %.2f mg*h/L\n', mean_auc, std_auc);
    fprintf(fid, '  Median:             %.2f mg*h/L\n', median_auc);
    fprintf(fid, '  Coefficient of Variation: %.1f%%\n\n', cv_auc);

    fprintf(fid, 'Distribution Span:\n');
    fprintf(fid, '  Minimum:            %.2f mg*h/L\n', min(AUC));
    fprintf(fid, '  Maximum:            %.2f mg*h/L\n', max(AUC));
    fprintf(fid, '  Range:              %.2f mg*h/L\n', max(AUC) - min(AUC));
    fprintf(fid, '  Median Absolute Deviation: %.2f mg*h/L\n\n', mad(AUC, 1));

    p5_auc = prctile(AUC, 5);
    p10_auc = prctile(AUC, 10);
    p25_auc = prctile(AUC, 25);
    p50_auc = prctile(AUC, 50);
    p75_auc = prctile(AUC, 75);
    p90_auc = prctile(AUC, 90);
    p95_auc = prctile(AUC, 95);
    iqr_auc = p75_auc - p25_auc;

    write_subheader('Percentile Distribution');
    fprintf(fid, '  5th percentile:     %.2f mg*h/L (5%% of samples below)\n', p5_auc);
    fprintf(fid, '  10th percentile:    %.2f mg*h/L\n', p10_auc);
    fprintf(fid, '  25th percentile:    %.2f mg*h/L (Q1)\n', p25_auc);
    fprintf(fid, '  50th percentile:    %.2f mg*h/L (Q2/Median)\n', p50_auc);
    fprintf(fid, '  75th percentile:    %.2f mg*h/L (Q3)\n', p75_auc);
    fprintf(fid, '  90th percentile:    %.2f mg*h/L\n', p90_auc);
    fprintf(fid, '  95th percentile:    %.2f mg*h/L (95%% of samples below)\n\n', p95_auc);

    fprintf(fid, 'Interquartile Range (IQR = Q3-Q1):\n');
    fprintf(fid, '  IQR:                %.2f mg*h/L\n', iqr_auc);
    fprintf(fid, '  Interpretation:     50%% of values fall within %.2f-%.2f mg*h/L\n\n', p25_auc, p75_auc);

    write_subheader('Therapeutic Window Analysis');
    try
        bench_gamelin = get5FU_literature_benchmarks('AUC', 'gamelin2008');
        pct_subtherapeutic = 100 * sum(AUC < bench_gamelin.subtherapeutic) / length(AUC);
        pct_optimal = 100 * sum(AUC >= bench_gamelin.optimal_low & AUC <= bench_gamelin.optimal_high) / length(AUC);
        pct_above = 100 * sum(AUC > bench_gamelin.optimal_high & AUC <= bench_gamelin.toxic) / length(AUC);
        pct_toxic = 100 * sum(AUC > bench_gamelin.toxic) / length(AUC);

        fprintf(fid, '\nGamelin et al. (2008) Clinical Targets:\n');
        fprintf(fid, '  Subtherapeutic (<%.0f mg*h/L):        %.1f%% (insufficient efficacy)\n', ...
            bench_gamelin.subtherapeutic, pct_subtherapeutic);
        fprintf(fid, '  OPTIMAL WINDOW (%.0f-%.0f mg*h/L):    %.1f%% (target)\n', ...
            bench_gamelin.optimal_low, bench_gamelin.optimal_high, pct_optimal);
        fprintf(fid, '  Above optimal (%.0f-%.0f mg*h/L):     %.1f%% (increased toxicity)\n', ...
            bench_gamelin.optimal_high, bench_gamelin.toxic, pct_above);
        fprintf(fid, '  Toxic range (>%.0f mg*h/L):          %.1f%% (high toxicity risk)\n\n', ...
            bench_gamelin.toxic, pct_toxic);

        fprintf(fid, 'Clinical Assessment:\n');
        if pct_optimal >= 60
            fprintf(fid, '  OK: %.1f%% of patients achieve optimal target (>=60%% criterion)\n', pct_optimal);
        elseif pct_optimal >= 40
            fprintf(fid, '  WARN: %.1f%% of patients achieve optimal target (40-60%% range)\n', pct_optimal);
        else
            fprintf(fid, '  FAIL: Only %.1f%% achieve target (<40%% criterion)\n', pct_optimal);
        end

        if pct_toxic <= 10
            fprintf(fid, '  OK: %.1f%% at high toxicity risk (<=10%% criterion)\n', pct_toxic);
        else
            fprintf(fid, '  WARN: %.1f%% patients exceed toxicity threshold\n', pct_toxic);
        end
    catch
        fprintf(fid, 'Therapeutic window analysis unavailable\n');
    end

    if ~isempty(param_sensitivities) && ~isempty(param_names)
        write_header('PARAMETER SENSITIVITY ANALYSIS');
        write_subheader('Spearman Rank Correlation with Cmax');

        [sorted_rho, order] = sort(abs(param_sensitivities(1,:)), 'descend');
        fprintf(fid, '\nParameters ranked by influence on Cmax:\n');
        for i = 1:min(10, length(order))
            fprintf(fid, '  %2d. %-30s   |rho| = %.3f\n', i, param_names{order(i)}, sorted_rho(i));
        end

        fprintf(fid, '\nInterpretation:\n');
        fprintf(fid, '  |rho| > 0.7:  Strong influence on Cmax variability\n');
        fprintf(fid, '  0.4 < |rho| < 0.7:  Moderate influence\n');
        fprintf(fid, '  |rho| < 0.4:  Weak influence\n');
    end

    write_header('STATISTICAL NOMENCLATURE AND DEFINITIONS');
    fprintf(fid, '\nPercentile:\n');
    fprintf(fid, '  The N-th percentile is a value below which N%% of observations fall.\n');
    fprintf(fid, '  Example: 25th percentile (Q1) means 25%% of samples fall below this value.\n\n');

    fprintf(fid, 'Interquartile Range (IQR):\n');
    fprintf(fid, '  IQR = Q3 - Q1 = 75th percentile - 25th percentile\n');
    fprintf(fid, '  Represents the range containing the middle 50%% of data.\n');
    fprintf(fid, '  More robust to outliers than mean +/- SD.\n\n');

    fprintf(fid, 'Coefficient of Variation (CV):\n');
    fprintf(fid, '  CV = (Standard Deviation / Mean) x 100%%\n');
    fprintf(fid, '  Normalized measure of spread; allows comparison across different scales.\n');
    fprintf(fid, '  CV < 15%%: Low variability\n');
    fprintf(fid, '  CV 15-30%%: Moderate variability\n');
    fprintf(fid, '  CV > 30%%: High variability\n\n');

    fprintf(fid, 'Median Absolute Deviation (MAD):\n');
    fprintf(fid, '  Robust measure of variability; less affected by outliers than SD.\n\n');

    write_header('RECOMMENDATIONS FOR MANUSCRIPT');
    fprintf(fid, '\n1. RESULTS SECTION:\n');
    fprintf(fid, '   "Monte Carlo analysis (N=%d samples) revealed a mean Cmax of %.2f +/- %.2f uM\n', ...
        n_samples, mean_cmax, std_cmax);
    fprintf(fid, '    (median %.2f uM, IQR %.2f-%.2f) and mean AUC of %.2f +/- %.2f mg*h/L\n', ...
        median_cmax, p25_cmax, p75_cmax, mean_auc, std_auc);
    fprintf(fid, '    (median %.2f mg*h/L, IQR %.2f-%.2f)."\n\n', median_auc, p25_auc, p75_auc);

    fprintf(fid, '2. FIGURE LEGENDS:\n');
    fprintf(fid, '   "Distribution plots show histogram of %d simulations with percentile overlays.\n', n_samples);
    fprintf(fid, '    Shaded zones indicate clinical efficacy (green) and toxicity (red) ranges."\n\n');

    fprintf(fid, '3. TABLE CAPTION:\n');
    fprintf(fid, '   "Descriptive statistics of pharmacokinetic parameters from %d Monte Carlo\n', n_samples);
    fprintf(fid, '    simulations. Values presented as mean +/- SD with median and IQR in parentheses."\n\n');

    fprintf(fid, '\n%s\n', repmat('=', 1, 80));
    fprintf(fid, 'Report generated by MC_5FU_PK_sensitivity.m\n');
    fprintf(fid, 'Contact: Pharmacokinetics and Personalized Medicine Laboratory\n');
    fprintf(fid, '%s\n', repmat('=', 1, 80));

    fclose(fid);
    fprintf('\nStatistical report saved: %s\n', report_filename);
end
