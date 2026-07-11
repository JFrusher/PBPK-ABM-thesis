function [time_min, concentrations] = integrate5FU(params, dosingRegimen, logger)
%integrate5FU Integrate the 5-FU PBPK/PD system on a uniform output grid.
%   [TIME_MIN, C] = integrate5FU(PARAMS, DOSINGREGIMEN, LOGGER) drives the
%   shared RHS (odes5FU_rhs) with either a stiff solver or a fixed-step Euler
%   loop, selected by PARAMS.solver_method:
%       'ode15s' - MATLAB stiff solver, integrated segment-by-segment across
%                  dosing breakpoints so infusion on/off edges are never
%                  stepped over; NonNegative on all states.
%       'fixed'  - forward Euler at PARAMS.fixed_timestep_min (reference).
%   Output is always returned on the uniform grid 0:dt_out:maxTime so the ABM
%   coupling and AUC integration see deterministic, evenly spaced rows.
%
%   Both paths call the identical RHS, so any difference in results reflects
%   the integrator, not the model.

    if nargin < 3, logger = []; end

    allowEval = isfield(params, 'allow_custom_function_eval') && params.allow_custom_function_eval;
    doseFn  = @(tt) calculateDosingRate(tt, dosingRegimen, allowEval);
    maxTime = max(dosingRegimen.end_time_min) + params.post_dose_observation_min;
    dt_out  = params.fixed_timestep_min;
    time_min = (0:dt_out:maxTime)';

    nStates = 18;
    y0 = zeros(nStates, 1);
    y0(17) = 1.0;   % dNTP pool starts fully replenished

    method = lower(string(params.solver_method));

    switch method
        case "ode15s"
            opts = odeset('NonNegative', 1:nStates, 'RelTol', 1e-5, 'AbsTol', 1e-7);
            % Breakpoints at every dosing edge so the solver restarts on
            % discontinuities instead of stepping across a bolus.
            edges = [0; dosingRegimen.start_time_min; dosingRegimen.end_time_min];
            if ismember('effective_end_time_min', dosingRegimen.Properties.VariableNames)
                edges = [edges; dosingRegimen.effective_end_time_min];
            end
            edges = unique(edges(edges >= 0 & edges <= maxTime));
            if edges(1) > 0,        edges = [0; edges];        end
            if edges(end) < maxTime, edges = [edges; maxTime]; end

            tsAll = 0; ysAll = y0';
            ycur = y0;
            for s = 1:numel(edges)-1
                a = edges(s); b = edges(s+1);
                if b <= a, continue; end
                [ts, ys] = ode15s(@(tt, yy) odes5FU_rhs(tt, yy, params, doseFn), [a b], ycur, opts);
                tsAll = [tsAll; ts(2:end)];      %#ok<AGROW>
                ysAll = [ysAll; ys(2:end, :)];   %#ok<AGROW>
                ycur  = ys(end, :)';
            end
            % Resample solver output onto the uniform reporting grid.
            [tsAll, iu] = unique(tsAll, 'stable');
            ysAll = ysAll(iu, :);
            Y = interp1(tsAll, ysAll, time_min, 'pchip');

        case "fixed"
            n = numel(time_min);
            Y = zeros(n, nStates);
            Y(1, :) = y0';
            for k = 2:n
                dydt = odes5FU_rhs(time_min(k-1), Y(k-1, :)', params, doseFn);
                yk = Y(k-1, :)' + dydt * dt_out;
                yk = max(yk, 0);            % Euler-only nonnegativity guard
                yk(17) = min(yk(17), 1.0);  % dNTP fraction cannot exceed baseline
                Y(k, :) = yk';
            end

        otherwise
            error('integrate5FU:UnknownSolver', ...
                'params.solver_method must be ''ode15s'' or ''fixed'' (got "%s").', method);
    end

    concentrations = unpackState5FU(Y, params);

    if ~isempty(logger)
        logger.info('integration_complete', struct( ...
            'solver_method', char(method), 'nTimePoints', numel(time_min), 'maxTime', maxTime));
    end
end

function C = unpackState5FU(Y, params)
% Map the integrated state matrix back to the named concentrations struct
% (the concentrations struct layout used downstream), filling the inert kidney/muscle/fat
% compartments with zeros and deriving the algebraic TS inhibition fraction.
    n = size(Y, 1);
    C.C_central     = Y(:, 1);
    C.C_peripheral  = Y(:, 2);
    C.C_liver       = Y(:, 3);
    C.C_kidney      = zeros(n, 1);
    C.C_tumor       = Y(:, 4);
    C.C_muscle      = zeros(n, 1);
    C.C_fat         = zeros(n, 1);
    C.C_FdUMP       = Y(:, 5);
    C.C_FdUTP       = Y(:, 6);
    C.C_FUTP        = Y(:, 7);
    C.C_tumor_FdUMP = Y(:, 8);
    C.C_tumor_FdUTP = Y(:, 9);
    C.C_tumor_FUTP  = Y(:, 10);
    C.C_DHFU        = Y(:, 11);
    C.C_FBAL        = Y(:, 12);
    C.excreted_5FU          = Y(:, 13);
    C.excreted_FBAL         = Y(:, 14);
    C.excreted_metabolites  = Y(:, 15);
    C.cumulative_input_5FU  = Y(:, 16);
    C.dNTP_pool_fraction    = Y(:, 17);
    C.DNA_damage_index      = Y(:, 18);
    C.TS_inhibition_fraction = C.C_tumor_FdUMP ./ (params.IC50_TS_tumor + C.C_tumor_FdUMP);
end
