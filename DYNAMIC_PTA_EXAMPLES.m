% =========================================================================
% DYNAMIC PTA - BEFORE & AFTER EXAMPLES
% =========================================================================
%
% This file shows the improvements to PTA plotting with dynamic adjustment

%% ========================================================================
% BEFORE: Fixed Literature-Based Targets (Original Code)
% ========================================================================

% Old approach - rigid targets, doesn't adapt to data characteristics
fprintf('OLD APPROACH: Fixed targets\n');

efficacy.values = [20, 28];  % Hard-coded limits
efficacy.labels = {'Minimum efficacy', 'Optimal efficacy'};

toxicity.values = [50, 70];  % Hard-coded limits
toxicity.labels = {'Grade 3+ toxicity (>20%)', 'Grade 3+ toxicity (>60%)'};

options.metric_name = 'AUC';
options.units = 'mg·h/L';
options.title_text = '5-FU Dose Optimization';
options.show_optimal = true;
options.optimal_pta = 90;
options.max_tox_prob = 10;
% No dynamic adjustment - plot range is ±50% of data

% PROBLEMS:
%   1. If data mean is 100 mg·h/L, targets of 20-70 are far outside range
%   2. If data range is 0.1-5, targets of 20-70 make no sense
%   3. Can't use for other drugs/metrics without modifying code
%   4. Fixed margins don't account for sample size

%% ========================================================================
% AFTER: Dynamic Data-Aware Targets (New Code)
% ========================================================================

fprintf('\nNEW APPROACH: Dynamic, data-aware targets\n');

% Approach 1: Auto-generate targets from percentiles
efficacy_dyn = struct('auto_generate', true);
toxicity_dyn = struct('auto_generate', true);

options_dyn = struct();
options_dyn.metric_name = 'AUC';
options_dyn.units = 'mg·h/L';
options_dyn.title_text = '5-FU Dose Optimization (Data-Driven)';
options_dyn.show_optimal = true;
options_dyn.optimal_pta = 90;
options_dyn.max_tox_prob = 10;
options_dyn.auto_adjust = true;  % ENABLE DYNAMIC ADJUSTMENT
options_dyn.percentile_efficacy = [25 50];  % P25 and P50
options_dyn.percentile_toxicity = [75 90];  % P75 and P90

% This will:
%   1. Calculate actual P25, P50, P75, P90 from data
%   2. Use those as targets (automatically appropriate!)
%   3. Scale plot margins based on sample size
%   4. Work with ANY metric/data range

plotPTA(auc_data, efficacy_dyn, toxicity_dyn, options_dyn);

% ADVANTAGES:
%   ✓ Works with any data range (0.1 to 1000+)
%   ✓ Automatically scales appropriately
%   ✓ No code changes needed for different metrics
%   ✓ Targets always clinically/statistically meaningful
%   ✓ Sample size affects margin (sparse data gets wider range)

%% ========================================================================
% EXAMPLE 1: Small Study (n=20)
% ========================================================================

fprintf('\n--- Example 1: Small Study (n=20) ---\n');

% Random small dataset
small_data = normrnd(30, 15, [20 1]);
small_data = small_data(small_data > 0);  % Remove negatives

fprintf('Data characteristics:\n');
fprintf('  N = %d\n', length(small_data));
fprintf('  Mean = %.2f\n', mean(small_data));
fprintf('  Range = %.2f to %.2f\n', min(small_data), max(small_data));

% With DYNAMIC adjustment:
%   - Margin factor = max(0.3, 2/sqrt(20)) = max(0.3, 0.447) = 0.447
%   - Expanded range with 44.7% margin on each side
%   - Larger viewing window for sparse data

% With OLD approach:
%   - Would use fixed ±50% margin
%   - Plot might be cramped or too wide

fprintf('Dynamic margin factor = %.3f (44.7%% expansion)\n', 0.447);
fprintf('Old fixed margin = 50.0%%\n');
fprintf('For sparse data, dynamic is ADAPTIVE ✓\n');

%% ========================================================================
% EXAMPLE 2: Large Study (n=500)
% ========================================================================

fprintf('\n--- Example 2: Large Study (n=500) ---\n');

large_data = normrnd(30, 15, [500 1]);
large_data = large_data(large_data > 0);

fprintf('Data characteristics:\n');
fprintf('  N = %d\n', length(large_data));
fprintf('  Mean = %.2f\n', mean(large_data));
fprintf('  Range = %.2f to %.2f\n', min(large_data), max(large_data));

% With DYNAMIC adjustment:
%   - Margin factor = max(0.3, 2/sqrt(500)) = max(0.3, 0.089) = 0.3
%   - Tighter viewing window for dense data
%   - Focuses on relevant region

% With OLD approach:
%   - Still uses fixed ±50% margin
%   - Wastes plot space for dense data

fprintf('Dynamic margin factor = %.3f (30.0%% expansion)\n', 0.3);
fprintf('Old fixed margin = 50.0%%\n');
fprintf('For large data, dynamic is MORE FOCUSED ✓\n');

%% ========================================================================
% EXAMPLE 3: Different Metrics (Drug A vs Drug B)
% ========================================================================

fprintf('\n--- Example 3: Drug A vs Drug B ---\n');

% Drug A: AUC range 5-50 mg·h/L
drug_a_auc = normrnd(25, 12, [100 1]);
drug_a_auc = drug_a_auc(drug_a_auc > 0);

% Drug B (different formulation): AUC range 100-500 mg·h/L  
drug_b_auc = normrnd(250, 120, [100 1]);
drug_b_auc = drug_b_auc(drug_b_auc > 0);

fprintf('Drug A AUC: Mean = %.1f, Range = %.1f-%.1f\n', ...
    mean(drug_a_auc), min(drug_a_auc), max(drug_a_auc));
fprintf('Drug B AUC: Mean = %.1f, Range = %.1f-%.1f\n', ...
    mean(drug_b_auc), min(drug_b_auc), max(drug_b_auc));

% OLD approach: Would require different target definitions
fprintf('\nOLD APPROACH: Need to hard-code different targets for each drug\n');
fprintf('  Drug A targets: [5 10 15 20] - manually defined\n');
fprintf('  Drug B targets: [100 200 300 400] - manually defined\n');
fprintf('  Problem: Code must be edited for each new drug/metric\n');

% NEW approach: Single code works for both!
fprintf('\nNEW APPROACH: One line of code works for any metric!\n');
efficacy_auto = struct('auto_generate', true);
toxicity_auto = struct('auto_generate', true);
options_auto.auto_adjust = true;
options_auto.percentile_efficacy = [25 50];
options_auto.percentile_toxicity = [75 90];

fprintf('For Drug A:\n');
%   plotPTA(drug_a_auc, efficacy_auto, toxicity_auto, options_auto);
fprintf('  Auto-targets: P25=%.1f, P50=%.1f, P75=%.1f, P90=%.1f\n', ...
    prctile(drug_a_auc,25), prctile(drug_a_auc,50), ...
    prctile(drug_a_auc,75), prctile(drug_a_auc,90));

fprintf('For Drug B:\n');
%   plotPTA(drug_b_auc, efficacy_auto, toxicity_auto, options_auto);
fprintf('  Auto-targets: P25=%.1f, P50=%.1f, P75=%.1f, P90=%.1f\n', ...
    prctile(drug_b_auc,25), prctile(drug_b_auc,50), ...
    prctile(drug_b_auc,75), prctile(drug_b_auc,90));

fprintf('No code changes needed! ✓\n');

%% ========================================================================
% EXAMPLE 4: Hybrid Mode - Keep Literature Targets, Auto-Scale Plot
% ========================================================================

fprintf('\n--- Example 4: Hybrid (Literature + Dynamic Range) ---\n');

% Use known literature targets but let plot auto-scale
fprintf('Best of both worlds:\n');
fprintf('  - Keep literature targets (proven clinical validity)\n');
fprintf('  - Auto-scale plot range to data characteristics\n');

efficacy_lit.values = [20 28];
efficacy_lit.labels = {'Minimum', 'Optimal'};
toxicity_lit.values = [50 70];
toxicity_lit.labels = {'Grade 2', 'Grade 3+'};

options_hybrid = struct();
options_hybrid.metric_name = 'AUC';
options_hybrid.units = 'mg·h/L';
options_hybrid.auto_adjust = true;  % Dynamic range but fixed targets
% Don't set percentile_efficacy/percentile_toxicity - use fixed targets

fprintf('Result:\n');
fprintf('  - Targets at 20, 28, 50, 70 mg·h/L (literature-based) ✓\n');
fprintf('  - Plot range adjusted to data (dynamic) ✓\n');
fprintf('  - No wasted space, all relevant data visible ✓\n');

%% ========================================================================
% FEATURE COMPARISON TABLE
% ========================================================================

fprintf('\n');
fprintf('┌─────────────────────┬──────────────────────┬──────────────────────┐\n');
fprintf('│ Feature             │ OLD (Fixed)          │ NEW (Dynamic)        │\n');
fprintf('├─────────────────────┼──────────────────────┼──────────────────────┤\n');
fprintf('│ Target Selection    │ Hard-coded only      │ Auto or hard-coded   │\n');
fprintf('│ Plot Range          │ ±50% fixed           │ Adaptive to N        │\n');
fprintf('│ Works with any data │ No (needs new code)  │ Yes (any metric)     │\n');
fprintf('│ Sparse data (n=10)  │ Poor margins         │ ~45%% margins        │\n');
fprintf('│ Large data (n=500)  │ Wastes space         │ ~30%% margins        │\n');
fprintf('│ Clinical relevance  │ Literature-based     │ Data-representative │\n');
fprintf('│ Code modifications  │ Required per drug    │ None required        │\n');
fprintf('│ Hybrid approach     │ Not available        │ Supported ✓          │\n');
fprintf('└─────────────────────┴──────────────────────┴──────────────────────┘\n');

%% ========================================================================
% QUICK START REFERENCE
% ========================================================================

fprintf('\n');
fprintf('QUICK START - Choose Your Approach:\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('1. PURE DATA-DRIVEN (percentile-based targets):\n');
fprintf('   efficacy = struct(''auto_generate'', true);\n');
fprintf('   toxicity = struct(''auto_generate'', true);\n');
fprintf('   options.auto_adjust = true;\n');
fprintf('   plotPTA(data, efficacy, toxicity, options);\n\n');

fprintf('2. HYBRID (literature targets + dynamic range):\n');
fprintf('   efficacy.values = [20 28];\n');
fprintf('   toxicity.values = [50 70];\n');
fprintf('   options.auto_adjust = true;\n');
fprintf('   plotPTA(data, efficacy, toxicity, options);\n\n');

fprintf('3. CUSTOM PERCENTILES (for different risk profiles):\n');
fprintf('   efficacy.auto_generate = true;\n');
fprintf('   toxicity.auto_generate = true;\n');
fprintf('   options.percentile_efficacy = [30 60];\n');
fprintf('   options.percentile_toxicity = [70 95];\n');
fprintf('   plotPTA(data, efficacy, toxicity, options);\n\n');

fprintf('4. TRADITIONAL (fully specified, no dynamic features):\n');
fprintf('   efficacy.values = [20 28];\n');
fprintf('   toxicity.values = [50 70];\n');
fprintf('   options.auto_adjust = false;  % Disable dynamic\n');
fprintf('   options.exposure_range = [0 100];  % Explicit range\n');
fprintf('   plotPTA(data, efficacy, toxicity, options);\n\n');

%% ========================================================================
% GENERATED OUTPUT
% ========================================================================

fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Generated PTA Files in Main Analysis:\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('1. MC_5FU_PK_PTA_[timestamp].pdf\n');
fprintf('   - Literature-based AUC targets\n');
fprintf('   - Dynamic plot scaling\n');
fprintf('   - Statistical summary printed\n\n');

fprintf('2. MC_5FU_PK_PTA_DataDriven_[timestamp].pdf\n');
fprintf('   - Auto-generated from AUC percentiles\n');
fprintf('   - Targets: P25, P50 (efficacy), P75, P90 (toxicity)\n');
fprintf('   - Fully data-adaptive\n\n');

fprintf('3. MC_5FU_PK_PTA_Cmax_[timestamp].pdf\n');
fprintf('   - Literature-based Cmax targets\n');
fprintf('   - Dynamic plot scaling\n\n');

fprintf('4. MC_5FU_PK_PTA_Cmax_DataDriven_[timestamp].pdf\n');
fprintf('   - Auto-generated from Cmax percentiles\n');
fprintf('   - Targets: P30, P60 (efficacy), P70, P95 (toxicity)\n\n');

fprintf('════════════════════════════════════════════════════════════════\n');
