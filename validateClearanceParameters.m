function validateClearanceParameters(params, logger)
% validateClearanceParameters  Ensure hepatic clearance respects flow limitation
if nargin < 2, logger = []; end
% compute CL_int = Vmax_DPD / Km_DPD (units L/min when Vmax in µmol/min and Km in mg/L)
% Note: we assume Vmax_DPD is in µmol/min and Km_DPD in mg/L; convert using MW if needed is performed upstream
if ~isfield(params,'Vmax_DPD') || ~isfield(params,'Km_DPD') || ~isfield(params,'Q_liver') || ~isfield(params,'CO')
    if ~isempty(logger), logger.warn('validate_params_incomplete', struct()); end
    return;
end
CL_int = params.Vmax_DPD / params.Km_DPD; % rough proxy (units L/min)
Q_liver = params.Q_liver * params.CO;
CL_hep = (Q_liver * CL_int) / (Q_liver + CL_int);
% If hepatic clearance is greater than blood flow (physically impossible) raise error
if CL_hep > (Q_liver + 1e-9)
    if ~isempty(logger), logger.error('hepatic_cl_violation', struct('CL_hep', CL_hep, 'Q_liver', Q_liver)); end
    error('MATLAB:err','Hepatic clearance violates blood flow constraint');
end
end