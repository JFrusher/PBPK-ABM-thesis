function report = audit_pbpk_numerics(inputCsv, outputPrefix, paramOverrides)
%AUDIT_PBPK_NUMERICS Intensive numerical continuity and sanity audit for PBPK run.
%   report = audit_pbpk_numerics(inputCsv)
%   report = audit_pbpk_numerics(inputCsv, outputPrefix)
%   report = audit_pbpk_numerics(inputCsv, outputPrefix, paramOverrides)

    if nargin < 2 || isempty(outputPrefix)
        [~, base, ~] = fileparts(char(inputCsv));
        outputPrefix = [base '_audit'];
    end
    if nargin < 3
        paramOverrides = struct();
    end

    results = run5FU_PBPK_Simulation(inputCsv, outputPrefix, paramOverrides);

    t = results.time_min(:);
    C = results.concentrations;

    dt = diff(t);
    report.time.min_dt = min(dt);
    report.time.max_dt = max(dt);
    report.time.mean_dt = mean(dt);
    report.time.n_nonpositive_dt = sum(dt <= 0);

    keySignals = {
        'C_central', 'C_tumor', 'C_liver', ...
        'C_FdUMP', 'C_tumor_FdUMP', 'C_tumor_FUTP'};

    jumpThresholdAbs = 25;   % µM/min equivalent step slope threshold
    jumpThresholdRel = 1.5;  % 150% relative step change

    jumpSummary = struct();
    for i = 1:numel(keySignals)
        name = keySignals{i};
        y = C.(name)(:);

        dy = diff(y);
        slope = zeros(size(dy));
        valid = dt > 0;
        slope(valid) = dy(valid) ./ dt(valid);

        yPrev = y(1:end-1);
        rel = zeros(size(dy));
        denom = max(abs(yPrev), 1e-9);
        rel = abs(dy) ./ denom;

        jumpMask = abs(slope) > jumpThresholdAbs | rel > jumpThresholdRel;
        zeroClampMask = (y(2:end) == 0) & (y(1:end-1) > 0);

        jumpSummary.(name).max_abs_slope = max(abs(slope));
        jumpSummary.(name).max_rel_step = max(rel);
        jumpSummary.(name).n_jump_flags = sum(jumpMask);
        jumpSummary.(name).n_zero_clamp_like = sum(zeroClampMask);
    end
    report.jumps = jumpSummary;

    % Monotonic checks for cumulative compartments
    report.cumulative.nonmonotonic_input = any(diff(C.cumulative_input_5FU(:)) < -1e-9);
    report.cumulative.nonmonotonic_excr_5fu = any(diff(C.excreted_5FU(:)) < -1e-9);
    report.cumulative.nonmonotonic_excr_fbal = any(diff(C.excreted_FBAL(:)) < -1e-9);
    report.cumulative.nonmonotonic_excr_met = any(diff(C.excreted_metabolites(:)) < -1e-9);

    report.metrics.AUC_central_mg_h_L = results.metrics.AUC_central_mg_h_L;
    report.metrics.AUC_tumor_mg_h_L = results.metrics.AUC_tumor_mg_h_L;
    report.metrics.Cmax_central_mg_L = results.metrics.Cmax_central_mg_L;

    report.mass.accounted_fraction = results.mass_balance.accounted_fraction;
    report.mass.unaccounted_umol = results.mass_balance.unaccounted_umol;

    % Aggregate flags
    flags = {};
    if report.time.n_nonpositive_dt > 0
        flags{end+1} = 'nonpositive_dt_detected'; %#ok<AGROW>
    end
    if report.time.max_dt > 0.5
        flags{end+1} = 'large_timestep_detected'; %#ok<AGROW>
    end
    if report.metrics.AUC_central_mg_h_L < 10
        flags{end+1} = 'auc_unusually_low'; %#ok<AGROW>
    elseif report.metrics.AUC_central_mg_h_L > 80
        flags{end+1} = 'auc_unusually_high'; %#ok<AGROW>
    end
    if report.mass.accounted_fraction < 0.90 || report.mass.accounted_fraction > 1.10
        flags{end+1} = 'mass_balance_outside_10pct'; %#ok<AGROW>
    end

    totalJumpFlags = 0;
    totalClampFlags = 0;
    for i = 1:numel(keySignals)
        totalJumpFlags = totalJumpFlags + report.jumps.(keySignals{i}).n_jump_flags;
        totalClampFlags = totalClampFlags + report.jumps.(keySignals{i}).n_zero_clamp_like;
    end
    report.aggregate.total_jump_flags = totalJumpFlags;
    report.aggregate.total_zero_clamp_like = totalClampFlags;
    if totalClampFlags > 0
        flags{end+1} = 'possible_numerical_clipping'; %#ok<AGROW>
    end

    report.flags = flags;

    fprintf('\n=== PBPK NUMERICAL AUDIT ===\n');
    fprintf('AUC central: %.4f mg·h/L | AUC tumor: %.4f mg·h/L\n', ...
        report.metrics.AUC_central_mg_h_L, report.metrics.AUC_tumor_mg_h_L);
    fprintf('dt min/mean/max: %.4f / %.4f / %.4f min\n', ...
        report.time.min_dt, report.time.mean_dt, report.time.max_dt);
    fprintf('Mass accounted fraction: %.4f\n', report.mass.accounted_fraction);
    fprintf('Jump flags: %d | Zero-clamp-like events: %d\n', ...
        report.aggregate.total_jump_flags, report.aggregate.total_zero_clamp_like);

    if isempty(report.flags)
        fprintf('Flags: none\n');
    else
        fprintf('Flags: %s\n', strjoin(report.flags, ', '));
    end
end
