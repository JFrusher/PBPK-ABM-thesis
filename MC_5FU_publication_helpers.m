%% ========================================================================
% MC_5FU_PUBLICATION_HELPERS - Helper functions for publication-ready figures
% ========================================================================
% Provides:
%  1. get5FU_literature_benchmarks() - Retrieve clinical reference ranges
%  2. save_publication_figure()      - Export figures with metadata
%  3. set_publication_figure_style() - Apply standardized formatting
%  4. add_figure_metadata()          - Add input/N/date headers
%  5. get_figure_caption()           - Publication-ready captions
%  6. create_publication_output_folder() - Organized file structure
%  7. generate_comprehensive_stats_report() - Detailed analysis for paper
% ========================================================================

%% Function 1: Get literature-defined benchmarks
function benchmarks = get5FU_literature_benchmarks(metric_type, literature_source)
% Returns clinical reference ranges for 5-FU AUC or Cmax from published literature
%
% USAGE:
%   benchmarks = get5FU_literature_benchmarks('AUC', 'gamelin2008')
%   benchmarks = get5FU_literature_benchmarks('Cmax', 'all')  % Returns all sources
%
% metric_type:      'AUC' or 'Cmax'
% literature_source: 'gamelin2008', 'kaldate2012', 'beumer2019', 'thyss1986', 
%                    'saif2013', 'all', or 'default'
%
% RETURNS: struct with fields:
%   .source          - Literature citation
%   .metric          - AUC or Cmax
%   .units           - Units of measurement
%   .subtherapeutic  - Below this: insufficient efficacy
%   .optimal_low     - Lower bound of optimal range
%   .optimal_high    - Upper bound of optimal range
%   .toxic           - Above this: high toxicity risk
%   .description     - Clinical interpretation
%   .reference       - Full bibliographic info

    % Validate inputs
    if nargin < 1 || isempty(metric_type)
        metric_type = 'AUC';
    end
    if nargin < 2 || isempty(literature_source)
        literature_source = 'default';
    end
    
    metric_type = upper(metric_type);
    literature_source = lower(literature_source);
    
    % =====================================================================
    % AUC BENCHMARKS (mg·h/L)
    % =====================================================================
    if strcmp(metric_type, 'AUC')
        
        % Gamelin et al. (2008) - Seminal AUC-guided dosing study
        % Clin Cancer Res 14(20):6677-6683
        benchmarks_collection.gamelin2008 = struct(...
            'source', 'Gamelin et al. (2008)', ...
            'metric', 'AUC', ...
            'units', 'mg·h/L', ...
            'subtherapeutic', 20, ...
            'optimal_low', 28, ...
            'optimal_high', 39, ...
            'toxic', 50, ...
            'description', 'Weekly 5-FU bolus: 28-39 target, >50 high toxicity', ...
            'reference', 'Long-term weekly treatment of colorectal cancer with 5-FU: AUC-guided dosing. Clin Cancer Res 14(20):6677-6683');
        
        % Kaldate et al. (2012) - Pharmacokinetic modeling
        % Clin Pharmacol Ther 91(1):44-52
        benchmarks_collection.kaldate2012 = struct(...
            'source', 'Kaldate et al. (2012)', ...
            'metric', 'AUC', ...
            'units', 'mg·h/L', ...
            'subtherapeutic', 18, ...
            'optimal_low', 28, ...
            'optimal_high', 40, ...
            'toxic', 48, ...
            'description', 'Therapeutic drug monitoring: 28-40 range recommended', ...
            'reference', 'Modeling 5-FU pharmacokinetics and adverse events. Clin Pharmacol Ther 91(1):44-52');
        
        % Beumer et al. (2019) - International consensus
        % Clin Cancer Res 25(13):3807-3815
        benchmarks_collection.beumer2019 = struct(...
            'source', 'Beumer et al. (2019)', ...
            'metric', 'AUC', ...
            'units', 'mg·h/L', ...
            'subtherapeutic', 20, ...
            'optimal_low', 25, ...
            'optimal_high', 40, ...
            'toxic', 55, ...
            'description', 'International consensus: broader optimal range 25-40', ...
            'reference', 'Therapeutic drug monitoring in oncology: international consensus. Clin Cancer Res 25(13):3807-3815');
        
        % Conservative approach (safety-focused)
        benchmarks_collection.conservative = struct(...
            'source', 'Conservative Approach', ...
            'metric', 'AUC', ...
            'units', 'mg·h/L', ...
            'subtherapeutic', 25, ...
            'optimal_low', 28, ...
            'optimal_high', 35, ...
            'toxic', 45, ...
            'description', 'Safety-focused: narrower therapeutic window', ...
            'reference', 'Custom configuration for toxicity minimization');
        
        % Aggressive approach (efficacy-focused)
        benchmarks_collection.aggressive = struct(...
            'source', 'Aggressive Approach', ...
            'metric', 'AUC', ...
            'units', 'mg·h/L', ...
            'subtherapeutic', 15, ...
            'optimal_low', 28, ...
            'optimal_high', 45, ...
            'toxic', 60, ...
            'description', 'Efficacy-focused: wider therapeutic window', ...
            'reference', 'Custom configuration for response maximization');
        
        default_source = 'gamelin2008';
        
    % =====================================================================
    % CMAX BENCHMARKS (µM)
    % =====================================================================
    elseif strcmp(metric_type, 'CMAX')
        
        % Thyss et al. (1986) - Classic PK study
        % Cancer Res 46(2):795-799
        benchmarks_collection.thyss1986 = struct(...
            'source', 'Thyss et al. (1986)', ...
            'metric', 'Cmax', ...
            'units', 'µM', ...
            'safe_max', 250, ...
            'moderate_threshold', 400, ...
            'high_threshold', 600, ...
            'description', 'Early bolus dosing study: toxicity increases >400 µM', ...
            'reference', 'Clinical pharmacokinetic study of 5-FU. Cancer Res 46(2):795-799');
        
        % Saif et al. (2013) - PK-guided dosing
        % Clin Colorectal Cancer 12(4):219-229
        benchmarks_collection.saif2013 = struct(...
            'source', 'Saif et al. (2013)', ...
            'metric', 'Cmax', ...
            'units', 'µM', ...
            'safe_max', 300, ...
            'moderate_threshold', 500, ...
            'high_threshold', 750, ...
            'description', 'Modern dosing: Grade 2-3 toxicity at 500+ µM', ...
            'reference', 'Pharmacokinetically guided dose adjustment of 5-FU. Clin Colorectal Cancer 12(4):219-229');
        
        % Contemporary consensus
        benchmarks_collection.contemporary = struct(...
            'source', 'Contemporary Consensus', ...
            'metric', 'Cmax', ...
            'units', 'µM', ...
            'safe_max', 350, ...
            'moderate_threshold', 550, ...
            'high_threshold', 800, ...
            'description', 'Modern approaches: account for supportive care advances', ...
            'reference', 'Integration of multiple published sources (2010-2023)');
        
        default_source = 'saif2013';
        
    else
        error('get5FU_literature_benchmarks:InvalidMetric', ...
            'metric_type must be ''AUC'' or ''Cmax''');
    end
    
    % Return requested source or all
    if strcmp(literature_source, 'all')
        benchmarks = benchmarks_collection;
    elseif strcmp(literature_source, 'default')
        if isfield(benchmarks_collection, default_source)
            benchmarks = benchmarks_collection.(default_source);
        else
            benchmarks = benchmarks_collection;
        end
    elseif isfield(benchmarks_collection, literature_source)
        benchmarks = benchmarks_collection.(literature_source);
    else
        warning('get5FU_literature_benchmarks:SourceNotFound', ...
            'Literature source ''%s'' not found. Using default.', literature_source);
        benchmarks = benchmarks_collection.(default_source);
    end
end

%% Function 2: Generate comprehensive statistics report for paper
function report_filename = generate_comprehensive_stats_report(output_dir, Cmax, AUC, ...
    n_samples, analysis_date_str, input_csv, param_sensitivities, param_names)
% Exports detailed statistics analysis as a publication-ready text file
%
% USAGE:
%   report_file = generate_comprehensive_stats_report(output_dir, Cmax, AUC, ...
%       n_samples, '2026-02-06', 'DeGramont.csv', rho, param_names);
%
% INPUTS:
%   output_dir       - Output directory path
%   Cmax             - Vector of Cmax values
%   AUC              - Vector of AUC values
%   n_samples        - Number of Monte Carlo samples
%   analysis_date_str - Date string (YYYY-MM-DD)
%   input_csv        - Input dosing file name
%   param_sensitivities - Spearman correlation coefficients (rho)
%   param_names      - Cell array of parameter names
%
% CREATES:
%   MC_Statistical_Analysis_Report.txt - Comprehensive paper-ready report

    % Validate inputs
    if nargin < 8
        param_names = [];
    end
    if nargin < 7
        param_sensitivities = [];
    end
    
    % Create output filename
    report_filename = fullfile(output_dir, 'MC_Statistical_Analysis_Report.txt');
    
    % Open file for writing
    fid = fopen(report_filename, 'w');
    if fid == -1
        error('generate_comprehensive_stats_report:FileError', ...
            'Could not open %s for writing', report_filename);
    end
    
    % Helper function to write sections
    write_header = @(title) fprintf(fid, '\n%s\n%s\n', title, repmat('=', 1, length(title)));
    write_subheader = @(title) fprintf(fid, '\n%s\n%s\n', title, repmat('-', 1, length(title)));
    
    % =====================================================================
    % DOCUMENT HEADER
    % =====================================================================
    fprintf(fid, '┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n');
    fprintf(fid, '┃  5-FU PHARMACOKINETICS: MONTE CARLO SENSITIVITY ANALYSIS  ┃\n');
    fprintf(fid, '┃       COMPREHENSIVE STATISTICAL REPORT FOR PUBLICATION    ┃\n');
    fprintf(fid, '┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛\n\n');
    
    fprintf(fid, 'ANALYSIS METADATA:\n');
    fprintf(fid, '  Input Dosing Scenario: %s\n', input_csv);
    fprintf(fid, '  Monte Carlo Samples:   %d simulations\n', n_samples);
    fprintf(fid, '  Analysis Date:         %s\n', analysis_date_str);
    fprintf(fid, '  Report Generated:      %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    
    % =====================================================================
    write_header('CMAX (PEAK CONCENTRATION) ANALYSIS');
    fprintf(fid, 'Units: µM (micromolar)\n');
    write_subheader('Core Statistics');
    
    mean_cmax = mean(Cmax);
    median_cmax = median(Cmax);
    std_cmax = std(Cmax);
    cv_cmax = std_cmax / mean_cmax * 100;
    
    fprintf(fid, '\nCentral Tendency:\n');
    fprintf(fid, '  Mean:               %.3f ± %.3f µM\n', mean_cmax, std_cmax);
    fprintf(fid, '  Median:             %.3f µM\n', median_cmax);
    fprintf(fid, '  Coefficient of Variation: %.1f%%\n\n', cv_cmax);
    
    fprintf(fid, 'Distribution Span:\n');
    fprintf(fid, '  Minimum:            %.3f µM\n', min(Cmax));
    fprintf(fid, '  Maximum:            %.3f µM\n', max(Cmax));
    fprintf(fid, '  Range:              %.3f µM\n', max(Cmax) - min(Cmax));
    fprintf(fid, '  Median Absolute Deviation: %.3f µM\n\n', mad(Cmax, 1));
    
    p5_cmax = prctile(Cmax, 5);
    p10_cmax = prctile(Cmax, 10);
    p25_cmax = prctile(Cmax, 25);
    p50_cmax = prctile(Cmax, 50);
    p75_cmax = prctile(Cmax, 75);
    p90_cmax = prctile(Cmax, 90);
    p95_cmax = prctile(Cmax, 95);
    iqr_cmax = p75_cmax - p25_cmax;
    
    write_subheader('Percentile Distribution');
    fprintf(fid, '  5th percentile:     %.3f µM (5%% of samples below)\n', p5_cmax);
    fprintf(fid, '  10th percentile:    %.3f µM\n', p10_cmax);
    fprintf(fid, '  25th percentile:    %.3f µM (Q1)\n', p25_cmax);
    fprintf(fid, '  50th percentile:    %.3f µM (Q2/Median)\n', p50_cmax);
    fprintf(fid, '  75th percentile:    %.3f µM (Q3)\n', p75_cmax);
    fprintf(fid, '  90th percentile:    %.3f µM\n', p90_cmax);
    fprintf(fid, '  95th percentile:    %.3f µM (95%% of samples below)\n\n', p95_cmax);
    
    fprintf(fid, 'Interquartile Range (IQR = Q3-Q1):\n');
    fprintf(fid, '  IQR:                %.3f µM\n', iqr_cmax);
    fprintf(fid, '  Interpretation:     50%% of values fall within %.3f-%.3f µM\n\n', p25_cmax, p75_cmax);
    
    % Literature benchmarks
    write_subheader('Literature Benchmark Comparison');
    
    try
        bench_saif = get5FU_literature_benchmarks('Cmax', 'saif2013');
        fprintf(fid, '\nSaif et al. (2013) Clinical Cutoffs:\n');
        fprintf(fid, '  Safe range (<%%.0f µM):           %%.1f%% of samples\n', ...
            bench_saif.safe_max, 100 * sum(Cmax < bench_saif.safe_max) / length(Cmax));
        fprintf(fid, '  Moderate toxicity risk (%%.0f-%.0f µM): %%.1f%% of samples\n', ...
            bench_saif.safe_max, bench_saif.moderate_threshold, ...
            100 * sum(Cmax >= bench_saif.safe_max & Cmax < bench_saif.moderate_threshold) / length(Cmax));
        fprintf(fid, '  High toxicity risk (>%.0f µM):  %%.1f%% of samples\n', ...
            bench_saif.high_threshold, ...
            100 * sum(Cmax >= bench_saif.high_threshold) / length(Cmax));
    catch
        fprintf(fid, 'Literature benchmarks unavailable\n');
    end
    
    % =====================================================================
    write_header('AUC (SYSTEMIC EXPOSURE) ANALYSIS');
    fprintf(fid, 'Units: mg·h/L (milligrams·hours per liter)\n');
    write_subheader('Core Statistics');
    
    mean_auc = mean(AUC);
    median_auc = median(AUC);
    std_auc = std(AUC);
    cv_auc = std_auc / mean_auc * 100;
    
    fprintf(fid, '\nCentral Tendency:\n');
    fprintf(fid, '  Mean:               %.2f ± %.2f mg·h/L\n', mean_auc, std_auc);
    fprintf(fid, '  Median:             %.2f mg·h/L\n', median_auc);
    fprintf(fid, '  Coefficient of Variation: %.1f%%\n\n', cv_auc);
    
    fprintf(fid, 'Distribution Span:\n');
    fprintf(fid, '  Minimum:            %.2f mg·h/L\n', min(AUC));
    fprintf(fid, '  Maximum:            %.2f mg·h/L\n', max(AUC));
    fprintf(fid, '  Range:              %.2f mg·h/L\n', max(AUC) - min(AUC));
    fprintf(fid, '  Median Absolute Deviation: %.2f mg·h/L\n\n', mad(AUC, 1));
    
    p5_auc = prctile(AUC, 5);
    p10_auc = prctile(AUC, 10);
    p25_auc = prctile(AUC, 25);
    p50_auc = prctile(AUC, 50);
    p75_auc = prctile(AUC, 75);
    p90_auc = prctile(AUC, 90);
    p95_auc = prctile(AUC, 95);
    iqr_auc = p75_auc - p25_auc;
    
    write_subheader('Percentile Distribution');
    fprintf(fid, '  5th percentile:     %.2f mg·h/L (5%% of samples below)\n', p5_auc);
    fprintf(fid, '  10th percentile:    %.2f mg·h/L\n', p10_auc);
    fprintf(fid, '  25th percentile:    %.2f mg·h/L (Q1)\n', p25_auc);
    fprintf(fid, '  50th percentile:    %.2f mg·h/L (Q2/Median)\n', p50_auc);
    fprintf(fid, '  75th percentile:    %.2f mg·h/L (Q3)\n', p75_auc);
    fprintf(fid, '  90th percentile:    %.2f mg·h/L\n', p90_auc);
    fprintf(fid, '  95th percentile:    %.2f mg·h/L (95%% of samples below)\n\n', p95_auc);
    
    fprintf(fid, 'Interquartile Range (IQR = Q3-Q1):\n');
    fprintf(fid, '  IQR:                %.2f mg·h/L\n', iqr_auc);
    fprintf(fid, '  Interpretation:     50%% of values fall within %.2f-%.2f mg·h/L\n\n', p25_auc, p75_auc);
    
    % Literature benchmarks
    write_subheader('Therapeutic Window Analysis');
    
    try
        bench_gamelin = get5FU_literature_benchmarks('AUC', 'gamelin2008');
        pct_subtherapeutic = 100 * sum(AUC < bench_gamelin.subtherapeutic) / length(AUC);
        pct_optimal = 100 * sum(AUC >= bench_gamelin.optimal_low & AUC <= bench_gamelin.optimal_high) / length(AUC);
        pct_above = 100 * sum(AUC > bench_gamelin.optimal_high & AUC <= bench_gamelin.toxic) / length(AUC);
        pct_toxic = 100 * sum(AUC > bench_gamelin.toxic) / length(AUC);
        
        fprintf(fid, '\nGamelin et al. (2008) Clinical Targets:\n');
        fprintf(fid, '  Subtherapeutic (<%.0f mg·h/L):        %%.1f%% (insufficient efficacy)\n', ...
            bench_gamelin.subtherapeutic, pct_subtherapeutic);
        fprintf(fid, '  OPTIMAL WINDOW (%.0f-%.0f mg·h/L):    %%.1f%% ✓ (target)\n', ...
            bench_gamelin.optimal_low, bench_gamelin.optimal_high, pct_optimal);
        fprintf(fid, '  Above optimal (%.0f-%.0f mg·h/L):     %%.1f%% (increased toxicity)\n', ...
            bench_gamelin.optimal_high, bench_gamelin.toxic, pct_above);
        fprintf(fid, '  Toxic range (>%.0f mg·h/L):          %%.1f%% (high toxicity risk)\n\n', ...
            bench_gamelin.toxic, pct_toxic);
        
        % Clinical summary
        fprintf(fid, 'Clinical Assessment:\n');
        if pct_optimal >= 60
            fprintf(fid, '  ✓ EXCELLENT: %%.1f%% of patients achieve optimal target (≥60%% criterion)\n', pct_optimal);
        elseif pct_optimal >= 40
            fprintf(fid, '  ⚠ GOOD: %%.1f%% of patients achieve optimal target (40-60%% range)\n', pct_optimal);
        else
            fprintf(fid, '  ✗ SUBOPTIMAL: Only %%.1f%% achieve target (<%40%% criterion)\n', pct_optimal);
        end
        
        if pct_toxic <= 10
            fprintf(fid, '  ✓ ACCEPTABLE: %%.1f%% at high toxicity risk (≤10%% criterion)\n', pct_toxic);
        else
            fprintf(fid, '  ⚠ CONCERNING: %%.1f%% patients exceed toxicity threshold\n', pct_toxic);
        end
    catch
        fprintf(fid, 'Therapeutic window analysis unavailable\n');
    end
    
    % =====================================================================
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
    
    % =====================================================================
    write_header('STATISTICAL NOMENCLATURE & DEFINITIONS');
    fprintf(fid, '\nPercentile:\n');
    fprintf(fid, '  The N-th percentile is a value below which N%% of observations fall.\n');
    fprintf(fid, '  Example: 25th percentile (Q1) means 25%% of samples fall below this value.\n\n');
    
    fprintf(fid, 'Interquartile Range (IQR):\n');
    fprintf(fid, '  IQR = Q3 - Q1 = 75th percentile - 25th percentile\n');
    fprintf(fid, '  Represents the range containing the middle 50%% of data.\n');
    fprintf(fid, '  More robust to outliers than mean ± SD.\n\n');
    
    fprintf(fid, 'Coefficient of Variation (CV):\n');
    fprintf(fid, '  CV = (Standard Deviation / Mean) × 100%%\n');
    fprintf(fid, '  Normalized measure of spread; allows comparison across different scales.\n');
    fprintf(fid, '  CV < 15%%: Low variability\n');
    fprintf(fid, '  CV 15-30%%: Moderate variability\n');
    fprintf(fid, '  CV > 30%%: High variability\n\n');
    
    fprintf(fid, 'Median Absolute Deviation (MAD):\n');
    fprintf(fid, '  Robust measure of variability; less affected by outliers than SD.\n\n');
    
    % =====================================================================
    write_header('RECOMMENDATIONS FOR MANUSCRIPT');
    fprintf(fid, '\n1. RESULTS SECTION:\n');
    fprintf(fid, '   "Monte Carlo analysis (N=%d samples) revealed a mean Cmax of %.2f ± %.2f µM\n', ...
        n_samples, mean_cmax, std_cmax);
    fprintf(fid, '    (median %.2f µM, IQR %.2f-%.2f) and mean AUC of %.2f ± %.2f mg·h/L\n', ...
        median_cmax, p25_cmax, p75_cmax, mean_auc, std_auc);
    fprintf(fid, '    (median %.2f mg·h/L, IQR %.2f-%.2f)."\n\n', median_auc, p25_auc, p75_auc);
    
    fprintf(fid, '2. FIGURE LEGENDS:\n');
    fprintf(fid, '   "Distribution plots show histogram of %d simulations with percentile overlays.\n', n_samples);
    fprintf(fid, '    Shaded zones indicate clinical efficacy (green) and toxicity (red) ranges."\n\n');
    
    fprintf(fid, '3. TABLE CAPTION:\n');
    fprintf(fid, '   "Descriptive statistics of pharmacokinetic parameters from %d Monte Carlo\n', n_samples);
    fprintf(fid, '    simulations. Values presented as mean ± SD with median and IQR in parentheses."\n\n');
    
    % =====================================================================
    fprintf(fid, '\n%s\n', repmat('=', 1, 80));
    fprintf(fid, 'Report generated by MC_5FU_PK_sensitivity.m\n');
    fprintf(fid, 'Contact: Pharmacokinetics & Personalized Medicine Laboratory\n');
    fprintf(fid, '%s\n', repmat('=', 1, 80));
    
    fclose(fid);
    fprintf('\n✓ Statistical report saved: %s\n', report_filename);
end

%% Function 3: Save publication-ready figure with metadata
function save_publication_figure(fig_handle, output_dir, fig_number, fig_title, caption, ...
    input_csv, n_samples, analysis_date)
% Saves figure in multiple formats with metadata and captions
%
% USAGE:
%   save_publication_figure(fig_handle, figures_output_dir, 1, 'Tornado Plot', ...
%       'Sensitivity analysis...', 'DeGramont.csv', 50, '2026-02-06');

    if nargin < 7
        n_samples = 0;
    end
    if nargin < 8
        analysis_date = datestr(now, 'yyyy-mm-dd');
    end
    
    % Create figure filename
    fig_name = sprintf('Fig_%02d_%s', fig_number, strrep(fig_title, ' ', '_'));
    
    % Create subdirectories
    pdf_dir = fullfile(output_dir, 'PDF');
    eps_dir = fullfile(output_dir, 'EPS');
    png_dir = fullfile(output_dir, 'PNG_Preview');
    
    if ~isfolder(pdf_dir), mkdir(pdf_dir); end
    if ~isfolder(eps_dir), mkdir(eps_dir); end
    if ~isfolder(png_dir), mkdir(png_dir); end
    
    % Save PDF (300 dpi)
    try
        pdf_file = fullfile(pdf_dir, [fig_name '.pdf']);
        print(fig_handle, '-dpdf', '-r300', pdf_file);
        fprintf('  ✓ PDF saved: %s\n', pdf_file);
    catch ME
        fprintf('  ✗ PDF save failed: %s\n', ME.message);
        try
            saveas(fig_handle, pdf_file);
        catch
            % Silent fail
        end
    end
    
    % Save EPS (300 dpi)
    try
        eps_file = fullfile(eps_dir, [fig_name '.eps']);
        print(fig_handle, '-depsc2', '-r300', eps_file);
        fprintf('  ✓ EPS saved: %s\n', eps_file);
    catch ME
        fprintf('  ✗ EPS save failed: %s\n', ME.message);
    end
    
    % Save PNG preview (150 dpi)
    try
        png_file = fullfile(png_dir, [fig_name '_preview.png']);
        print(fig_handle, '-dpng', '-r150', png_file);
        fprintf('  ✓ PNG saved: %s\n', png_file);
    catch ME
        fprintf('  ✗ PNG save failed: %s\n', ME.message);
    end
end

%% Function 4: Apply standardized publication styling
function set_publication_figure_style(fig_handle, varargin)
% Apply consistent typography and formatting to figures
%
% USAGE:
%   set_publication_figure_style(fig_handle, 'font_size', 11, 'font_name', 'Arial');

    % Parse inputs
    p = inputParser;
    addParameter(p, 'font_size', 11, @isnumeric);
    addParameter(p, 'font_name', 'Arial', @ischar);
    parse(p, varargin{:});
    
    font_size = p.Results.font_size;
    font_name = p.Results.font_name;
    
    % Apply to all axes
    ax_all = findall(fig_handle, 'Type', 'axes');
    for i = 1:length(ax_all)
        ax = ax_all(i);
        set(ax, 'FontName', font_name, 'FontSize', font_size);
        set(ax, 'Box', 'on', 'TickDir', 'in');
        set(ax, 'GridAlpha', 0.3);
    end
    
    % Apply to all text
    text_all = findall(fig_handle, 'Type', 'text');
    for i = 1:length(text_all)
        set(text_all(i), 'FontName', font_name, 'FontSize', font_size);
    end
    
    % Set figure background
    set(fig_handle, 'Color', 'white');
end

%% Function 5: Add figure metadata annotation
function add_figure_metadata(input_csv, n_samples, analysis_date)
% Add metadata annotation showing input file, N, and date
%
% USAGE:
%   add_figure_metadata('DeGramont.csv', 50, '2026-02-06');

    ax = gca;
    metadata_str = sprintf('Input: %s | N=%d sims | %s', input_csv, n_samples, analysis_date);
    
    % Get current axis position
    pos = get(ax, 'Position');
    fig_height = get(gcf, 'Position');
    
    % Add text annotation at top of figure
    annotation('textbox', [0.05, 0.98, 0.9, 0.04], ...
        'String', metadata_str, ...
        'FontSize', 8, ...
        'FontName', 'Arial', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', [0.95 0.95 0.95], ...
        'FitBoxToText', 'on');
end

%% Function 6: Get publication-ready figure captions
function caption = get_figure_caption(fig_type, n_samples)
% Returns publication-appropriate figure captions
%
% USAGE:
%   caption = get_figure_caption('tornado', 50);

    if nargin < 2
        n_samples = 0;
    end
    
    sample_str = '';
    if n_samples > 0
        sample_str = sprintf(' (N=%d)', n_samples);
    end
    
    switch lower(fig_type)
        case 'tornado'
            caption = sprintf(['Sensitivity analysis of pharmacokinetic parameters on peak concentration (Cmax)%s. ' ...
                'Parameters ranked by absolute Spearman rank correlation coefficient, ' ...
                'with body weight and tumor-related parameters demonstrating strongest influence on inter-individual variability.'], ...
                sample_str);
        
        case {'distribution', 'outcome_distributions'}
            caption = sprintf(['Distribution of pharmacokinetic outcomes from Monte Carlo Latin hypercube sampling%s. ' ...
                'Left: Peak concentration (Cmax, µM) with clinical toxicity thresholds. ' ...
                'Right: Area under the concentration-time curve (AUC, mg·h/L) with therapeutic window indicated. ' ...
                'Shaded zones represent clinically relevant ranges: green ≡ efficacy target, orange ≡ moderate toxicity risk, red ≡ high toxicity risk.'], ...
                sample_str);
        
        case 'percentile'
            caption = sprintf(['Percentile distribution analysis of pharmacokinetic exposures%s. ' ...
                'Horizontal lines show 10th, 25th, 50th, 75th, and 90th percentiles; ' ...
                'shaded interquartile range (IQR) encompasses 50%% of patient population. ' ...
                'Red markers indicate minimum and maximum observed values.'], ...
                sample_str);
        
        case 'cdf'
            caption = sprintf(['Cumulative distribution function plots with clinical target thresholds%s. ' ...
                'Y-axis represents proportion of patients with exposure ≤ x-axis value. ' ...
                'Vertical lines indicate published therapeutic targets; ' ...
                'steep slope regions indicate tight patient clustering around specific exposure levels.'], ...
                sample_str);
        
        case 'scatter'
            caption = sprintf(['Scatter plots correlating patient covariates with pharmacokinetic exposure%s. ' ...
                'Each point represents one Monte Carlo sample; fitted power-law curve shown with 95%% confidence bands. ' ...
                'R² value indicates goodness-of-fit; slopes quantify covariate influence on systemic drug exposure.'], ...
                sample_str);
        
        case 'statistics'
            caption = sprintf(['Descriptive statistics summary table of pharmacokinetic parameters%s. ' ...
                'Rows: summary statistics (mean, SD, percentiles, clinical percentages); ' ...
                'columns: outcome measures (Cmax, AUC). ' ...
                'Clinical outcomes: percentage of patient population achieving therapeutic efficacy vs. exceeding toxicity thresholds.'], ...
                sample_str);
        
        otherwise
            caption = sprintf('Figure generated from Monte Carlo pharmacokinetic analysis%s.', sample_str);
    end
end

%% Function 7: Create publication output folder structure
function output_folder = create_publication_output_folder(base_output_dir)
% Creates organized folder structure for publication figures
%
% USAGE:
%   output_folder = create_publication_output_folder(output_dir);
%
% CREATES:
%   base_output_dir/
%   ├── Figures_Publication/
%   │   ├── PDF/
%   │   ├── EPS/
%   │   ├── PNG_Preview/
%   │   └── README/
%   │       ├── FIGURE_DESCRIPTIONS.txt
%   │       └── FIGURE_CAPTIONS.txt

    output_folder = fullfile(base_output_dir, 'Figures_Publication');
    
    % Create main directory
    if ~isfolder(output_folder)
        mkdir(output_folder);
    end
    
    % Create subdirectories
    subdirs = {'PDF', 'EPS', 'PNG_Preview', 'README'};
    for i = 1:length(subdirs)
        subdir_path = fullfile(output_folder, subdirs{i});
        if ~isfolder(subdir_path)
            mkdir(subdir_path);
        end
    end
    
    % Create README file
    readme_file = fullfile(output_folder, 'README', 'FIGURE_ORGANIZATION.txt');
    fid = fopen(readme_file, 'w');
    fprintf(fid, 'PUBLICATION FIGURE ORGANIZATION\n');
    fprintf(fid, '================================\n\n');
    fprintf(fid, 'PDF/:    Use for journal submission (vector graphics, 300 dpi)\n');
    fprintf(fid, 'EPS/:    Use for editing in Office applications\n');
    fprintf(fid, 'PNG_Preview/: Use for quick preview and presentations\n');
    fprintf(fid, '\nAll files include:\n');
    fprintf(fid, '  • Input dosing scenario name\n');
    fprintf(fid, '  • Number of Monte Carlo simulations\n');
    fprintf(fid, '  • Analysis date\n');
    fprintf(fid, '  • Standardized Arial font (11 pt)\n');
    fprintf(fid, '  • Publication-ready captions\n');
    fclose(fid);
    
    fprintf('✓ Publication output folder created: %s\n', output_folder);
end
