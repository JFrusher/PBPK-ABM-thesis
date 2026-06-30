function C = updateConcentrations(C, dCdt, t, dt)
% Update all concentrations using Euler method
%
% Inputs:
%   C    - Current concentration structure
%   dCdt - Rate of change structure
%   t    - Current time index
%   dt   - Time step (minutes)
%
% Output:
%   C - Updated concentration structure

    % Update all compartments and metabolites
    fields = fieldnames(dCdt);

    for i = 1:length(fields)
        fieldName = fields{i};
        C.(fieldName)(t) = C.(fieldName)(t-1) + dCdt.(fieldName) * dt;

        % Ensure non-negative concentrations (numerical stability)
        C.(fieldName)(t) = max(C.(fieldName)(t), 0);
    end

    % Cap dNTP_pool_fraction at 1.0 (pool cannot exceed baseline)
    C.dNTP_pool_fraction(t) = min(1.0, C.dNTP_pool_fraction(t));
end
