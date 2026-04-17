% QUICK TEST TEMPLATE: Switch between Fixed and Adaptive Solvers
% Use this to test both methods side-by-side

clear all; close all;

% Get default params
p = struct();
p.solver_method = 'fixed';        % or 'adaptive'
p.fixed_timestep_min = 0.1;       % 0.1 minute steps
p.enable_ode_diagnostics = true;  % Enable detailed logging

fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║         5-FU PBPK SOLVER METHOD COMPARISON TEST               ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

fprintf('Testing solver_method = ''%s''\n', p.solver_method);
fprintf('Input file: BM500.csv\n\n');

% ═══════════════════════════════════════════════════════════════════════
% TEST 1: FIXED TIMESTEP (Stable Reference)
% ═══════════════════════════════════════════════════════════════════════

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('TEST 1: FIXED TIMESTEP METHOD (Stable)\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

tic;
try
    results_fixed = run5FU_PBPK_Simulation('BM500.csv', 'TEST_FIXED');
    elapsed_fixed = toc;
    
    fprintf('\n✓ FIXED method completed in %.1f seconds\n', elapsed_fixed);
    fprintf('  Simulated AUC: %.2f mg·h/L\n', results_fixed.metrics.AUC_central_mg_h_L);
    fprintf('  Time points: %d\n', length(results_fixed.time_min));
    
catch ME
    fprintf('\n✗ FIXED method FAILED:\n');
    fprintf('  %s\n', ME.message);
    results_fixed = [];
    elapsed_fixed = NaN;
end

% ═══════════════════════════════════════════════════════════════════════
% TEST 2: ADAPTIVE TIMESTEP (Fast but Verify)
% ═══════════════════════════════════════════════════════════════════════

fprintf('\n\n═══════════════════════════════════════════════════════════════\n');
fprintf('TEST 2: ADAPTIVE TIMESTEP METHOD (Fast)\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

p.solver_method = 'adaptive';  % Switch to adaptive

tic;
try
    results_adaptive = run5FU_PBPK_Simulation('BM500.csv', 'TEST_ADAPTIVE');
    elapsed_adaptive = toc;
    
    fprintf('\n✓ ADAPTIVE method completed in %.1f seconds\n', elapsed_adaptive);
    fprintf('  Simulated AUC: %.2f mg·h/L\n', results_adaptive.metrics.AUC_central_mg_h_L);
    fprintf('  Time points: %d\n', length(results_adaptive.time_min));
    
catch ME
    fprintf('\n✗ ADAPTIVE method FAILED:\n');
    fprintf('  %s\n', ME.message);
    results_adaptive = [];
    elapsed_adaptive = NaN;
end

% ═══════════════════════════════════════════════════════════════════════
% COMPARISON
% ═══════════════════════════════════════════════════════════════════════

if ~isempty(results_fixed) && ~isempty(results_adaptive)
    fprintf('\n\n╔═══════════════════════════════════════════════════════════════╗\n');
    fprintf('║                    COMPARISON RESULTS                          ║\n');
    fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');
    
    fprintf('COMPUTATIONAL PERFORMANCE:\n');
    fprintf('  Fixed method:     %.1f seconds (%d points)\n', elapsed_fixed, length(results_fixed.time_min));
    fprintf('  Adaptive method:  %.1f seconds (%d points)\n', elapsed_adaptive, length(results_adaptive.time_min));
    fprintf('  Speedup: %.1f×\n\n', elapsed_fixed / elapsed_adaptive);
    
    fprintf('PHARMACOKINETIC AGREEMENT:\n');
    AUC_fixed = results_fixed.metrics.AUC_central_mg_h_L;
    AUC_adaptive = results_adaptive.metrics.AUC_central_mg_h_L;
    AUC_error_pct = abs(AUC_adaptive - AUC_fixed) / AUC_fixed * 100;
    
    fprintf('  Fixed AUC:        %.2f mg·h/L\n', AUC_fixed);
    fprintf('  Adaptive AUC:     %.2f mg·h/L\n', AUC_adaptive);
    fprintf('  Difference:       %.2f%% (should be < 5%%)\n\n', AUC_error_pct);
    
    if AUC_error_pct < 5
        fprintf('  ✓ EXCELLENT: Methods agree within 5%% (acceptable)\n');
    elseif AUC_error_pct < 10
        fprintf('  ⚠️ WARNING: Methods differ by %.1f%% (investigate)\n', AUC_error_pct);
    else
        fprintf('  ✗ ERROR: Methods differ by %.1f%% (major problem!)\n', AUC_error_pct);
    end
    
    fprintf('\n  Cmax Fixed:       %.2f µM\n', results_fixed.metrics.Cmax_central);
    fprintf('  Cmax Adaptive:    %.2f µM\n', results_adaptive.metrics.Cmax_central);
    Cmax_error_pct = abs(results_adaptive.metrics.Cmax_central - results_fixed.metrics.Cmax_central) / results_fixed.metrics.Cmax_central * 100;
    fprintf('  Difference:       %.2f%%\n\n', Cmax_error_pct);
    
    % Compare shapes of curves
    fprintf('PLOT COMPARISON:\n');
    fprintf('  Fixed plots:      TEST_FIXED_5FU_compartments.png\n');
    fprintf('  Adaptive plots:   TEST_ADAPTIVE_5FU_compartments.png\n');
    fprintf('  → Open both and compare for artifacts, drops, or spikes\n\n');
    
    fprintf('DIAGNOSIS:\n');
    if AUC_error_pct < 5 && Cmax_error_pct < 5
        fprintf('  ✓✓✓ ADAPTIVE METHOD IS SAFE ✓✓✓\n');
        fprintf('  Both methods agree. Adaptive can be used for faster runs.\n');
    else
        fprintf('  ⚠️ ADAPTIVE METHOD NEEDS INVESTIGATION\n');
        fprintf('  Error exceeds 5%%. Check:\n');
        fprintf('  1. Bolus phase captures in first 5 minutes\n');
        fprintf('  2. No negative slopes (concentration drops) after peak\n');
        fprintf('  3. Tumor < Central at all times\n');
    end
    
else
    fprintf('\nOne or both methods failed. Check error messages above.\n');
end

fprintf('\n\n');
