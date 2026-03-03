% Test: 500 mg bolus IV
% Expected: AUC ≈ 20-30 mg·h/L (Schalhorn benchmark)

dosing_csv = 'test_500mg_bolus.csv';
% Content: start_time_min,end_time_min,dosing_type,dose_amount,infusion_rate
%          0,1,bolus,500,0

results = run5FU_PBPK_Simulation(dosing_csv, 'TEST_500mg');

fprintf('\\nTEST RESULT:\\n');
fprintf('Expected AUC: 20-30 mg·h/L\\n');
fprintf('Simulated AUC: %.1f mg·h/L\\n', results.metrics.AUC_central_mg_h_L);

if results.metrics.AUC_central_mg_h_L >= 20 && results.metrics.AUC_central_mg_h_L <= 30
    fprintf('✓ PASS: AUC within acceptable range\\n');
else
    fprintf('✗ FAIL: AUC outside range, check parameters\\n');
end