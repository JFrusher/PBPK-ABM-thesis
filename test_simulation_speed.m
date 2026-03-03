function test_simulation_speed()
% TEST_SIMULATION_SPEED  Benchmark the optimized simulation speed
%
% Usage: test_simulation_speed()
%
% Runs a single simulation and reports timing information.
%
% Expected results after optimization:
%   - Should complete in 1-2 minutes (vs ~15 min before)
%   - AUC values should be consistent
%   - Cmax values should be consistent

fprintf('\n╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║          SIMULATION SPEED BENCHMARK TEST                     ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');

% Create a simple dosing regimen
dosing_table = table([0], [30], {'bolus'}, [400], ...
    'VariableNames', {'start_time_min','end_time_min','dosing_type','dose_mg_per_m2'});

% Write to temporary file
temp_dosing = fullfile(tempdir, 'test_dosing.csv');
writetable(dosing_table, temp_dosing);

% Set output directory
output_dir = fullfile(tempdir, 'test_sim_output');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

fprintf('Test Configuration:\n');
fprintf('  Dosing: 400 mg/m² bolus at t=0-30 min\n');
fprintf('  Simulation duration: ~2 hours\n');
fprintf('  Output directory: %s\n\n', output_dir);

% Time the simulation
tic;
fprintf('Running simulation...\n');

try
    results = run5FU_PBPK_Simulation(temp_dosing, fullfile(output_dir, 'test_sim'), struct());
    
    elapsed = toc;
    
    fprintf('\n╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('║                    RESULTS                                   ║\n');
    fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');
    
    fprintf('Elapsed time: %.1f seconds (%.1f minutes)\n\n', elapsed, elapsed/60);
    
    if elapsed < 120
        fprintf('✓ EXCELLENT: Simulation completed in <2 min\n');
        fprintf('  This is the optimized speed!\n');
    elseif elapsed < 300
        fprintf('✓ GOOD: Simulation completed in <5 min\n');
        fprintf('  Optimization is working, but timesteps could be coarser.\n');
    else
        fprintf('⚠ WARNING: Simulation took >5 min\n');
        fprintf('  Optimization may not be active or timesteps are too fine.\n');
    end
    
    fprintf('\nKey pharmacokinetic results:\n');
    fprintf('  AUC (mg·h/L):           %.2f\n', results.metrics.AUC_central_mg_h_L);
    fprintf('  Cmax (µM):              %.2f\n', results.metrics.Cmax_central);
    fprintf('  Time to Cmax (min):     %.1f\n', results.metrics.Tmax_central);
    fprintf('  Terminal half-life (h): %.2f\n', results.metrics.t_half_terminal);
    
    fprintf('\nFiles generated:\n');
    csv_files = dir(fullfile(output_dir, '*.csv'));
    for i = 1:length(csv_files)
        fprintf('  - %s\n', csv_files(i).name);
    end
    
    fprintf('\n✓ Test completed successfully!\n\n');
    
catch ME
    elapsed = toc;
    fprintf('\n✗ ERROR: Simulation failed after %.1f seconds\n', elapsed);
    fprintf('  Message: %s\n\n', ME.message);
    rethrow(ME);
end

% Cleanup
% delete(temp_dosing);

end
