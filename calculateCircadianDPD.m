function factor = calculateCircadianDPD(hourOfDay, params)
% calculateCircadianDPD  Circadian modulation factor for DPD activity
% hourOfDay: hour in 0-24 (can be fractional)
if nargin < 2 || isempty(params)
    params = initialise5FUParameters();
end
if ~isfield(params,'DPD_mean'), params.DPD_mean = 1.0; end
if ~isfield(params,'DPD_amplitude'), params.DPD_amplitude = 0.37; end
if ~isfield(params,'DPD_acrophase'), params.DPD_acrophase = 1; end
factor = params.DPD_mean + params.DPD_amplitude * cos(2*pi*(hourOfDay - params.DPD_acrophase)/24);
end