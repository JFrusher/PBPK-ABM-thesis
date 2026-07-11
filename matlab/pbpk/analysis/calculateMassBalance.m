function mass_balance = calculateMassBalance(time_min, C, params)
% Summarize mass balance at the final time point (µmol)

V_central = params.V_central;
V_peripheral = params.V_peripheral;
V_liver = params.V_liver;   % absolute L (was x BW -> 127 L bug)
V_tumor = params.V_tumor;   % absolute L

parent_in_compartments = C.C_central(end) * V_central + ...
    C.C_peripheral(end) * V_peripheral + ...
    C.C_liver(end) * V_liver + ...
    C.C_tumor(end) * V_tumor;

metabolites_systemic = (C.C_FdUMP(end) + C.C_FdUTP(end) + C.C_FUTP(end)) * V_central;
metabolites_tumor = (C.C_tumor_FdUMP(end) + C.C_tumor_FdUTP(end) + C.C_tumor_FUTP(end)) * V_tumor;
catabolites = (C.C_DHFU(end) + C.C_FBAL(end)) * V_central;

excreted_total = C.excreted_5FU(end) + C.excreted_FBAL(end) + C.excreted_metabolites(end);
total_input = C.cumulative_input_5FU(end);
total_accounted = parent_in_compartments + metabolites_systemic + metabolites_tumor + catabolites + excreted_total;

mass_balance.time_end_min = time_min(end);
mass_balance.total_input_umol = total_input;
mass_balance.parent_in_compartments_umol = parent_in_compartments;
mass_balance.metabolites_systemic_umol = metabolites_systemic;
mass_balance.metabolites_tumor_umol = metabolites_tumor;
mass_balance.catabolites_umol = catabolites;
mass_balance.excreted_umol = excreted_total;
mass_balance.total_accounted_umol = total_accounted;
mass_balance.unaccounted_umol = total_input - total_accounted;
mass_balance.accounted_fraction = total_accounted / max(total_input, 1e-9);

% Strict closed-system check: with the clean RHS (no clamps) input should equal
% compartments + metabolites + catabolites + excreted to within integration error.
if total_input > 1e-6 && abs(mass_balance.accounted_fraction - 1) > 0.01
    warning('calculateMassBalance:NotClosed', ...
        'Mass balance off by %.2f%% (accounted %.4f of %.3f umol input) - check kinetics/params.', ...
        100 * (1 - mass_balance.accounted_fraction), mass_balance.accounted_fraction, total_input);
end
end
