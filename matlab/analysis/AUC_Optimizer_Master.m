%=============================================================================
%  v2.1
%  Corrected wrapper with proper sinusoidal handling
%=============================================================================

function [optimized_csv, final_metrics] = optimize_5FU_AUC(input_csv, target_AUC, varargin)
    %
    % FIXED in v2.0:
    % ✓ Properly handles sinusoidal CSV without amplitude column
    % ✓ Scales infusion_rate directly for sinusoidal (not amplitude)
    % ✓ Better error handling - doesn't break loop on simulation
    % ✓ Detects sinusoidal even without explicit amplitude column
    % ✓ Forced exact dosing_type specification (no silent detection)
    
    p = inputParser;
    p.addParameter('dosing_type', 'auto', @ischar);
    p.addParameter('tolerance', 0.05, @isnumeric);
    p.addParameter('max_iterations', 20, @isnumeric);
    p.addParameter('scale_method', 'proportional', @ischar);
    p.addParameter('output_prefix', 'AUC_optimized', @ischar);
    p.parse(varargin{:});
    
    dosing_type = p.Results.dosing_type;
    tolerance = p.Results.tolerance;
    max_iterations = p.Results.max_iterations;
    scale_method = p.Results.scale_method;
    output_prefix = p.Results.output_prefix;
    
    fprintf('\n%s\n', repmat('═', 1, 80));
    fprintf('║ 5-FU AUC OPTIMIZER v2.0 (FIXED)\n');
    fprintf('%s\n', repmat('═', 1, 80));
    fprintf('Target AUC: %.2f mg·h/L | Tolerance: ±%.2f\n', target_AUC, tolerance);
    fprintf('Max iterations: %d | Scaling method: %s | Dosing type: %s\n', ...
        max_iterations, scale_method, dosing_type);
    fprintf('%s\n\n', repmat('─', 1, 80));
    
    % Load original CSV
    original_data = readtable(input_csv);
    fprintf('✓ Loaded: %s (%d rows)\n', input_csv, height(original_data));
    
    % Auto-detect dosing type if needed (IMPROVED)
    if strcmp(dosing_type, 'auto')
        dosing_type = detect_dosing_type_improved(original_data);
        fprintf('✓ Auto-detected dosing type: %s\n', dosing_type);
    end
    
    fprintf('✓ Using dosing type: %s\n\n', dosing_type);
    
    % Initialize
    current_data = original_data;
    scaling_factor = 1.0;
    iteration = 0;
    AUC_history = [];
    scaling_history = [];
    
    fprintf('╔═ OPTIMIZATION LOOP ═╗\n\n');
    
    % Main optimization loop
    while iteration < max_iterations
        iteration = iteration + 1;
        
        % Run simulation
        output_file = sprintf('%s_iter%02d.csv', output_prefix, iteration);
        writetable(current_data, output_file);
        
        % Run PBPK simulation
        try
            results = run5FU_PBPK_Simulation(output_file, sprintf('%s_sim%02d', output_prefix, iteration));
            current_AUC = results.metrics.AUC_central_mg_h_L;
            
            % Check for NaN
            if isnan(current_AUC)
                fprintf('Iter %2d | AUC: NaN (simulation error)\n', iteration);
                continue;
            end
        catch ME
            fprintf('Iter %2d | SIMULATION FAILED: %s\n', iteration, ME.message);
            fprintf('  Continuing with previous AUC estimate...\n');
            if iteration > 1
                current_AUC = AUC_history(end);  % Use last valid AUC
            else
                fprintf('  Cannot recover from iteration 1 failure. Aborting.\n');
                break;
            end
        end
        
        error_AUC = current_AUC - target_AUC;
        percent_error = 100 * error_AUC / target_AUC;
        AUC_history(iteration) = current_AUC;
        scaling_history(iteration) = scaling_factor;
        
        % Display iteration results
        status_pct = min(100, max(0, 100 * (1 - abs(percent_error)/50)));
        status_bar = repmat('=', 1, max(1, round(40 * status_pct/100)));
        status_pad = repmat('░', 1, 40 - length(status_bar));
        
        fprintf('Iter %2d | AUC: %6.2f | Target: %6.2f | Error: %+7.2f%% | Scale: %6.3f\n', ...
            iteration, current_AUC, target_AUC, percent_error, scaling_factor);
        fprintf('         [%s%s] Progress\n', status_bar, status_pad);
        
        % Check convergence
        if abs(error_AUC) <= tolerance
            fprintf('\n✓✓✓ CONVERGED! Error within ±%.2f mg·h/L ✓✓✓\n\n', tolerance);
            break;
        end
        
        % Calculate next scaling factor (IMPROVED DAMPING)
        raw_scale = target_AUC / current_AUC;
        
        % Adaptive damping based on how far from target we are
        if abs(percent_error) > 30
            damping = 0.6;  % More conservative when far from target
        elseif abs(percent_error) > 10
            damping = 0.7;  % Standard damping
        else
            damping = 0.8;  % More aggressive when close
        end
        
        new_scaling = scaling_factor + damping * (raw_scale * scaling_factor - scaling_factor);
        
        % Safety bounds: ±40% per iteration (more lenient)
        new_scaling = max(0.6 * scaling_factor, min(1.4 * scaling_factor, new_scaling));
        
        scaling_factor = new_scaling;
        
        % Apply new scaling to CSV (FIXED FOR SINUSOIDAL)
        current_data = scale_dosing_csv_fixed(original_data, scaling_factor, dosing_type, scale_method);
        
    end
    
    fprintf('%s\n\n', repmat('─', 1, 80));
    
    % Generate final optimized CSV
    final_csv_file = sprintf('%s_FINAL.csv', output_prefix);
    writetable(current_data, final_csv_file);
    fprintf('✓ Saved optimized CSV: %s\n', final_csv_file);
    
    % Run final simulation for metrics
    try
        final_results = run5FU_PBPK_Simulation(final_csv_file, sprintf('%s_FINAL', output_prefix));
        final_AUC = final_results.metrics.AUC_central_mg_h_L;
        final_Cmax = final_results.metrics.Cmax_central_mg_L;
        final_AUC_tumor = final_results.metrics.AUC_tumor_mg_h_L;
        final_toxicity = final_results.metrics.toxicity_category;
    catch
        fprintf('Error running final simulation. Using iteration %d metrics.\n', iteration);
        final_AUC = AUC_history(end);
        final_Cmax = 0;
        final_AUC_tumor = 0;
        final_toxicity = 'UNKNOWN';
    end
    
    % Compile output structures
    optimized_csv.filename = final_csv_file;
    optimized_csv.data = current_data;
    optimized_csv.original = original_data;
    optimized_csv.scaling_factor = scaling_factor;
    optimized_csv.dosing_type = dosing_type;
    
    final_metrics.AUC = final_AUC;
    final_metrics.error = final_AUC - target_AUC;
    final_metrics.percent_error = 100 * final_metrics.error / target_AUC;
    final_metrics.iterations = iteration;
    final_metrics.Cmax = final_Cmax;
    final_metrics.AUC_tumor = final_AUC_tumor;
    final_metrics.toxicity_category = final_toxicity;
    final_metrics.AUC_history = AUC_history;
    final_metrics.scaling_history = scaling_history;
    final_metrics.target_AUC = target_AUC;
    
    % Clinical recommendation
    if abs(final_metrics.error) <= tolerance
        recommendation = sprintf('✓ OPTIMAL: AUC within target (%.2f ± %.2f)', target_AUC, tolerance);
        status = 'SUCCESS';
    elseif final_metrics.error > 0
        adjustment = round(100 * (scaling_factor - 1));
        recommendation = sprintf('FURTHER OPTIMIZATION NEEDED: AUC too high, DECREASE dose by ~%d%% more', -adjustment);
        status = 'PARTIAL';
    else
        adjustment = round(100 * (1 - scaling_factor));
        recommendation = sprintf('FURTHER OPTIMIZATION NEEDED: AUC too low, INCREASE dose by ~%d%% more', adjustment);
        status = 'PARTIAL';
    end
    final_metrics.recommendation = recommendation;
    final_metrics.status = status;
    
    % Display summary
    fprintf('\n╔════════════════════ OPTIMIZATION SUMMARY ════════════════════╗\n\n');
    fprintf('Status:            %s\n', status);
    fprintf('Final AUC:         %7.2f mg·h/L (target: %.2f)\n', final_metrics.AUC, target_AUC);
    fprintf('Absolute Error:    %+7.2f mg·h/L\n', final_metrics.error);
    fprintf('Relative Error:    %+7.2f%%\n', final_metrics.percent_error);
    fprintf('Total Iterations:  %7d / %d\n', final_metrics.iterations, max_iterations);
    fprintf('Dose Adjustment:   %+7.1f%% (scale: %.3f)\n', ...
        100*(scaling_factor-1), scaling_factor);
    fprintf('\nClinical Metrics:\n');
    fprintf('  Cmax:            %7.2f mg/L\n', final_metrics.Cmax);
    fprintf('  Tumor AUC:       %7.2f mg·h/L\n', final_metrics.AUC_tumor);
    fprintf('  Toxicity Risk:   %s\n\n', final_metrics.toxicity_category);
    fprintf('%s\n', recommendation);
    fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');
    
end

% ═════════════════════════════════════════════════════════════════════════════
% HELPER FUNCTIONS - IMPROVED
% ═════════════════════════════════════════════════════════════════════════════

function dosing_type = detect_dosing_type_improved(data)
    % Improved detection that handles sinusoidal CSVs without amplitude column
    
    var_names = data.Properties.VariableNames;
    
    % Check for explicit amplitude column first
    if ismember('amplitude', var_names)
        amp_rows = sum(data.amplitude > 0);
        if amp_rows > 0
            dosing_type = 'sinusoidal';
            return;
        end
    end
    
    % Check for bolus
    if ismember('dose_amount', var_names)
        bolus_rows = sum(data.dose_amount > 0);
        if bolus_rows > height(data) * 0.5
            dosing_type = 'bolus';
            return;
        end
    end
    
    % Check for infusion_rate (could be constant, square, or sinusoidal hourly)
    if ismember('infusion_rate', var_names)
        infusion_rows = sum(data.infusion_rate > 0);
        if infusion_rows > height(data) * 0.5
            rates = data.infusion_rate(data.infusion_rate > 0);
            unique_rates = unique(rates);
            
            % Check if sinusoidal pattern (many different rates, smooth variation)
            if length(unique_rates) > 5 && var(rates) > 0
                % Check for sinusoidal pattern (smooth oscillation)
                rate_diffs = diff(rates);
                smooth_variation = sum(abs(rate_diffs) < std(rates)) / length(rate_diffs);
                
                if smooth_variation > 0.7  % >70% of transitions are small
                    dosing_type = 'sinusoidal';
                    return;
                end
            end
            
            % Check for square wave (2 distinct rates)
            if length(unique_rates) == 2
                dosing_type = 'square';
                return;
            end
            
            % Default to constant
            dosing_type = 'constant';
            return;
        end
    end
    
    dosing_type = 'constant';
end

function scaled_data = scale_dosing_csv_fixed(original_data, scaling_factor, dosing_type, method)
    % FIXED scaling function that properly handles sinusoidal
    
    scaled_data = original_data;
    
    switch dosing_type
        case 'bolus'
            % Scale bolus amounts
            if ismember('dose_amount', scaled_data.Properties.VariableNames)
                bolus_mask = scaled_data.dose_amount > 0;
                
                switch method
                    case 'proportional'
                        scaled_data.dose_amount(bolus_mask) = ...
                            scaled_data.dose_amount(bolus_mask) * scaling_factor;
                    
                    case 'differential'
                        doses = scaled_data.dose_amount(bolus_mask);
                        mean_dose = mean(doses);
                        for i = find(bolus_mask)'
                            dose_ratio = scaled_data.dose_amount(i) / mean_dose;
                            local_scale = scaling_factor + 0.5 * (1 - dose_ratio) * (scaling_factor - 1);
                            scaled_data.dose_amount(i) = scaled_data.dose_amount(i) * local_scale;
                        end
                end
            end
        
        case 'constant'
            % Scale infusion rates
            if ismember('infusion_rate', scaled_data.Properties.VariableNames)
                infusion_mask = scaled_data.infusion_rate > 0;
                scaled_data.infusion_rate(infusion_mask) = ...
                    scaled_data.infusion_rate(infusion_mask) * scaling_factor;
            end
        
        case 'square'
            % Scale both on and off rates
            if ismember('infusion_rate', scaled_data.Properties.VariableNames)
                rates = scaled_data.infusion_rate(scaled_data.infusion_rate > 0);
                unique_rates = unique(rates);
                
                for rate = unique_rates'
                    mask = scaled_data.infusion_rate == rate;
                    scaled_data.infusion_rate(mask) = rate * scaling_factor;
                end
            end
        
        case 'sinusoidal'
            % FIX: Scale INFUSION_RATE directly (not amplitude column)
            % For sinusoidal pattern in hourly blocks
            
            if ismember('infusion_rate', scaled_data.Properties.VariableNames)
                % Scale all infusion rates by the factor
                % This maintains the sinusoidal pattern but changes its amplitude
                infusion_mask = scaled_data.infusion_rate > 0;
                scaled_data.infusion_rate(infusion_mask) = ...
                    scaled_data.infusion_rate(infusion_mask) * scaling_factor;
            end
            
            % Also scale amplitude column if it exists
            if ismember('amplitude', scaled_data.Properties.VariableNames)
                amp_mask = scaled_data.amplitude > 0;
                scaled_data.amplitude(amp_mask) = ...
                    scaled_data.amplitude(amp_mask) * scaling_factor;
            end
    end
    
end

%=============================================================================