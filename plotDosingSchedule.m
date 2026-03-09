function plotDosingSchedule(csvFilename)
% ================================================================================
% STANDALONE DOSING SCHEDULE PLOTTING SCRIPT
% ================================================================================
%
% FUNCTION: plotDosingSchedule(csvFilename)
%
% PURPOSE:  Standalone script to load a 5-FU dosing regimen CSV file and
%           immediately plot the infusion rate profile. Outputs a PNG figure
%           with filename derived from input CSV name. Displays figure in
%           MATLAB GUI.
%
% INPUT:    csvFilename (string) - Path or filename of dosing regimen CSV
%           Examples: 'dosing_schedule.csv'
%                     'C:\Data\patient_123_protocol.csv'
%                     '../schedules/bolus_regimen.csv'
%
% OUTPUT:   1. PNG file: [csvFilename_without_extension]_dosing_plot.png
%           2. Figure displayed in MATLAB figure window
%           3. Console output with file save confirmation
%
% USAGE:    plotDosingSchedule('my_dosing_schedule.csv');
%           plotDosingSchedule('C:\Clinical_Data\Protocol_A_24hr.csv');
%
% CSV FORMAT REQUIREMENTS:
%   Required columns:
%   - start_time_min (numeric, minutes from t=0)
%   - end_time_min (numeric, minutes from t=0)
%   - dosing_type (string: 'bolus', 'continuous', or 'sinusoidal')
%
%   Type-specific columns:
%   - For 'bolus': dose_mg or dose_amount (numeric)
%   - For 'continuous': infusion_rate (numeric, mg/min)
%   - For 'sinusoidal': mean_rate, amplitude, frequency_per_min (numeric)
%
% ================================================================================
    % ===== INPUT VALIDATION =====
    if ~isstring(csvFilename) && ~ischar(csvFilename)
        error('Input must be a string filename or character array');
    end
    
    % Check if file exists
    if ~isfile(csvFilename)
        error('File not found: %s', csvFilename);
    end
    
    fprintf('Loading dosing regimen from: %s\n', csvFilename);
    
    % ===== LOAD CSV =====
    try
        dosingRegimen = readtable(csvFilename);
        fprintf('Successfully loaded %d dosing regimen(s)\n', height(dosingRegimen));
    catch ME
        error('Failed to read CSV file: %s\n%s', csvFilename, ME.message);
    end
    
    % ===== EXTRACT OUTPUT PREFIX FROM FILENAME =====
    % Remove path and extension
    [~, filename, ~] = fileparts(csvFilename);
    outputPrefix = filename;  % Use CSV filename (without extension) as prefix
    
    % Determine simulation end time
    simulationEndTime = max(dosingRegimen.end_time_min);
    
    fprintf('Simulation duration: %.0f minutes (%.1f hours)\n', ...
            simulationEndTime, simulationEndTime/60);
    
    % ===== CREATE FIGURE =====
    fprintf('Creating dosing strategy plot...\n');
    
    % Time discretization: 1 point every minute for smooth visualization
    timeVector = (0:1:simulationEndTime)';
    infusionRate = zeros(size(timeVector));
    
    % Calculate infusion rate at each time point
    for t = 1:length(timeVector)
        currentTime = timeVector(t);
        infusionRate(t) = calculateDosingRate(currentTime, dosingRegimen);
    end
    
    % Convert time to hours for more readable x-axis
    timeHours = timeVector / 60;
    
    % ===== CREATE FIGURE WINDOW =====
    fig = figure('Position', [100 100 1200 600], 'Color', 'white', ...
                 'Name', sprintf('Dosing Schedule: %s', outputPrefix), ...
                 'NumberTitle', 'off');
    
    % ===== MAIN PLOT =====
    ax1 = axes(fig, 'Position', [0.1 0.15 0.75 0.75]);
    
    % Plot infusion rate
    fill(timeHours, infusionRate, [0.2 0.6 1.0], ...
     'FaceAlpha', 0.6, 'EdgeColor', [0.1 0.3 0.8], 'LineWidth', 2);
    hold on;
    plot(timeHours, infusionRate, 'Color', [0.1 0.3 0.8], 'LineWidth', 2.5);
    
    % Formatting
    xlabel('Time (hours)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Infusion Rate (mg/min)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('%s - 5-FU Dosing Schedule', outputPrefix), ...
          'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    grid(ax1, 'minor');
    set(ax1, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
    
    % ===== STATISTICS BOX =====
    totalDose = trapz(timeHours, infusionRate);  % Integrate rate to get dose (mg)
    maxRate = max(infusionRate);
    durationActive = sum(infusionRate > 0.01);  % Minutes with non-negligible infusion
    
    statsText = sprintf('DOSING STATISTICS\nTotal Dose: %.0f mg\nPeak Rate: %.2f mg/min\nDuration: %.1f h\nNum. Regimens: %d', ...
        totalDose, maxRate, durationActive/60, height(dosingRegimen));

    
    annotation('textbox', [0.78 0.45 0.20 0.30], ...
        'String', statsText, ...
        'FontSize', 10, 'FontName', 'monospace', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', ...
        'LineWidth', 1.5, 'FitBoxToText', 'on');
    
    % ===== ANNOTATE DOSING EVENTS =====
    annotateDosingEvents(ax1, dosingRegimen);
    
    % ===== FORMATTING =====
    xlim([0 simulationEndTime/60]);
    ylim([0 max(infusionRate) * 1.1]);
    
    % Better x-axis ticks (every 6 hours if long, every 1 hour if short)
    if simulationEndTime > 360
        xticks(0:6:ceil(simulationEndTime/60));
    else
        xticks(0:1:ceil(simulationEndTime/60));
    end
    
    % Add light background
    set(ax1, 'Color', [0.98 0.98 0.99]);
    
    % ===== SAVE FIGURE =====
    outputFilename = sprintf('%s_dosing_plot.png', outputPrefix);
    saveas(fig, outputFilename);
    fprintf('\n✓ Dosing strategy plot saved: %s\n', outputFilename);
    
    % Show figure
    fprintf('✓ Figure displayed in MATLAB window\n\n');
    
end

% ================================================================================
% HELPER FUNCTION: Annotate Dosing Events
% ================================================================================

function annotateDosingEvents(ax, dosingRegimen)
    % Add color-coded legend for dosing types
    
    colors_by_type = containers.Map(...
        {'bolus', 'continuous', 'sinusoidal'}, ...
        {[1 0.3 0.3], [0.3 1 0.3], [1 0.7 0.3]});  % Red, Green, Orange
    
    % Create legend entries
    h1 = line(nan, nan, 'Color', [1 0.3 0.3], 'LineStyle', '-', 'LineWidth', 3);
    h2 = line(nan, nan, 'Color', [0.3 1 0.3], 'LineStyle', '-', 'LineWidth', 3);
    h3 = line(nan, nan, 'Color', [1 0.7 0.3], 'LineStyle', '-', 'LineWidth', 3);
    
    legend([h1 h2 h3], {'Bolus Events', 'Continuous Infusion', 'Sinusoidal/Chronotherapy'}, ...
           'Location', 'northeast', 'FontSize', 10, 'EdgeColor', 'black');
    
end

% ================================================================================
% HELPER FUNCTION: Calculate Dosing Rate
% ================================================================================

function rate = calculateDosingRate(currentTime, dosingRegimen)
    % Calculate the current 5-FU dosing rate (mg/min) at a given time point
    
    rate = 0;
    MW_5FU = 130.08;
    
    for i = 1:height(dosingRegimen)
        startTime = dosingRegimen.start_time_min(i);
        endTime_raw = dosingRegimen.end_time_min(i);
        
        % Robust string handling for dosing_type
        if iscell(dosingRegimen.dosing_type)
            dosingType = lower(strtrim(dosingRegimen.dosing_type{i}));
        else
            dosingType = lower(strtrim(string(dosingRegimen.dosing_type(i))));
        end
        
        switch dosingType
            case 'bolus'
                % Bolus: rapid injection
                duration = max(endTime_raw - startTime, 0.1);  % Minimum 0.1 min
                effectiveEndTime = startTime + duration;
                
                if currentTime >= startTime && currentTime < effectiveEndTime
                    % Extract dose amount
                    dose_mg = NaN;
                    if ismember('dose_mg', dosingRegimen.Properties.VariableNames)
                        dose_mg = dosingRegimen.dose_mg(i);
                    elseif ismember('dose_amount', dosingRegimen.Properties.VariableNames)
                        dose_mg = dosingRegimen.dose_amount(i);
                    end
                    
                    % Calculate and add rate
                    if ~isnan(dose_mg) && dose_mg > 0
                        dose_rate_mg_per_min = dose_mg / duration;
                        rate = rate + dose_rate_mg_per_min;
                    end
                end
                
            case {'constant', 'continuous'}
                % Continuous infusion: steady rate
                if currentTime >= startTime && currentTime <= endTime_raw
                    rate_mg_per_min = dosingRegimen.infusion_rate(i);
                    rate = rate + rate_mg_per_min;
                end
                
            case 'sinusoidal'
                % Sinusoidal/chronotherapy: time-varying rate
                if currentTime >= startTime && currentTime <= endTime_raw
                    mean_rate = dosingRegimen.mean_rate(i);
                    amplitude = dosingRegimen.amplitude(i);
                    frequency = dosingRegimen.frequency_per_min(i);
                    
                    % Sinusoidal function: rate(t) = mean + amplitude*sin(2πft)
                    rate_mg_per_min = mean_rate + amplitude * sin(2 * pi * frequency * currentTime);
                    rate_mg_per_min = max(rate_mg_per_min, 0);  % No negative rates
                    
                    rate = rate + rate_mg_per_min;
                end
        end
    end
end

