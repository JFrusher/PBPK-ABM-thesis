function circadianFactor = calculateCircadianDPD(hourOfDay, params)
% Calculate circadian modulation factor for DPD enzyme activity
%
% Based on Harris et al. (1990) Cancer Research:
%   - Peak DPD activity at 1 AM (acrophase)
%   - Trough at 1 PM (13:00)
%   - 1.74-fold variation (peak/trough ratio)
%
% Uses cosine function: Activity = mean + amplitude * cos(2π(t - acrophase)/24)
%
% Inputs:
%   hourOfDay - Hour of day (0-24, where 0 = midnight)
%   params    - Parameter structure with DPD circadian parameters
%
% Output:
%   circadianFactor - Multiplier for DPD Vmax (range ~0.63 to 1.37)

    % Cosine function with 24-hour period
    % Phase shift so peak occurs at acrophase (1 AM = hour 1)
    circadianFactor = params.DPD_mean + ...
                      params.DPD_amplitude * cos(2 * pi * (hourOfDay - params.DPD_acrophase) / 24);

    % Ensure factor is always positive (safety check)
    circadianFactor = max(circadianFactor, 0.1);

end
