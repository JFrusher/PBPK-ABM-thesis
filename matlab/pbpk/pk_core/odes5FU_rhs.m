function dydt = odes5FU_rhs(t, y, params, doseFn)
%odes5FU_rhs Solver-agnostic right-hand side for the 5-FU PBPK/PD system.
%   dydt = odes5FU_rhs(T, Y, PARAMS, DOSEFN) returns the time derivatives of
%   the 18-state 5-FU model at time T (minutes) given state vector Y and a
%   dosing-rate function DOSEFN(T) -> umol/min. Pure function: no persistent
%   state, no I/O, no time-step assumptions. Usable by ode15s/ode45 and by a
%   fixed-step Euler loop alike (see integrate5FU).
%
%   State vector order (see stateIndex5FU / integrate5FU):
%     1 C_central      2 C_peripheral  3 C_liver      4 C_tumor
%     5 C_FdUMP        6 C_FdUTP       7 C_FUTP
%     8 C_tumor_FdUMP  9 C_tumor_FdUTP 10 C_tumor_FUTP
%    11 C_DHFU        12 C_FBAL
%    13 excreted_5FU  14 excreted_FBAL 15 excreted_metabolites
%    16 cumulative_input_5FU          17 dNTP_pool_fraction  18 DNA_damage_index
%
%   Concentrations in uM; excretion/input accumulators in umol; PD states
%   dimensionless. Kidney/muscle/fat carry no dynamics (dC/dt = 0) and are
%   held out of the integrated vector.

    MW_5FU = 130.08;

    % --- unpack ---
    C_central     = y(1);  C_peripheral = y(2);  C_liver = y(3);  C_tumor = y(4);
    C_FdUMP       = y(5);  C_FdUTP      = y(6);  C_FUTP  = y(7);
    C_tumor_FdUMP = y(8);  C_tumor_FdUTP= y(9);  C_tumor_FUTP = y(10);
    C_DHFU        = y(11); C_FBAL       = y(12);
    dNTP          = y(17); DNA          = y(18);

    dosingRate = doseFn(t);   % umol/min

    % --- circadian DPD modulation ---
    hourOfDay    = mod(t, 1440) / 60;
    DPD_factor   = calculateCircadianDPD(hourOfDay, params);
    Vmax_current = params.Vmax_DPD * DPD_factor;         % umol/min (canonical)

    % --- hepatic + extra-hepatic DPD (Michaelis-Menten, mg/L basis) ---
    C_central_mg_L      = C_central * MW_5FU / 1000;
    C_liver_mg_L        = C_liver   * MW_5FU / 1000;
    Vmax_current_mg_min = Vmax_current * MW_5FU / 1000;

    rate_DPD_liver_mg_min    = (Vmax_current_mg_min * C_liver_mg_L) / (params.Km_DPD + C_liver_mg_L);
    rate_DPD_liver           = rate_DPD_liver_mg_min / (MW_5FU / 1000);
    rate_DPD_systemic_mg_min = (Vmax_current_mg_min * params.k_extra_hepatic_fraction * C_central_mg_L) / ...
                               (params.Km_DPD + C_central_mg_L);
    rate_DPD_systemic        = rate_DPD_systemic_mg_min / (MW_5FU / 1000);

    % --- flow-limited distribution (umol/min) ---
    Q_liver_L = params.Q_liver * params.CO;   V_liver_L = params.V_liver;
    Q_tumor_L = params.Q_tumor * params.CO;   V_tumor_L = params.V_tumor;

    rate_central_to_periph = (params.k_cp * params.V_central)    * C_central;
    rate_periph_to_central = (params.k_pc * params.V_peripheral) * C_peripheral;
    rate_central_to_liver  = Q_liver_L * C_central;
    rate_liver_to_central  = (Q_liver_L * C_liver) / params.Kp_liver;   % efflux respects Kp
    rate_central_to_tumor  = Q_tumor_L * C_central;
    rate_tumor_to_central  = (Q_tumor_L * C_tumor) / params.Kp_tumor;

    % --- systemic anabolite formation (saturable), cell-cycle + TS circadian ---
    rate_form_FdUMP_sys = (params.Vmax_UMPS * C_central) / (params.Km_UMPS + C_central);
    rate_form_FUTP_sys  = (params.Vmax_RR   * C_central) / (params.Km_RR   + C_central);
    rate_form_FdUTP_sys = (params.Vmax_CDA  * C_central) / (params.Km_CDA  + C_central);
    rate_form_FdUMP_sys = rate_form_FdUMP_sys * params.cycle_modulation_factor;
    rate_form_FdUTP_sys = rate_form_FdUTP_sys * params.cycle_modulation_factor;
    rate_form_FUTP_sys  = rate_form_FUTP_sys  * params.cycle_modulation_factor;
    TS_circadian        = 1.0 + params.TS_acrophase * cos(2*pi*(hourOfDay - params.TS_peak_hour)/24);
    rate_form_FdUMP_sys = rate_form_FdUMP_sys * TS_circadian;

    % --- tumor anabolite formation (mass rate, umol/min) ---
    rate_form_FdUMP_tumor = params.k_form_FdUMP * params.tumor_metabolite_factor * C_tumor * V_tumor_L;
    rate_form_FdUTP_tumor = params.k_form_FdUTP * params.tumor_metabolite_factor * C_tumor * V_tumor_L;
    rate_form_FUTP_tumor  = params.k_form_FUTP  * params.tumor_metabolite_factor * C_tumor * V_tumor_L;

    % --- anabolite elimination (mass rate, umol/min) ---
    rate_elim_FdUMP       = params.k_elim_FdUMP * C_FdUMP       * params.V_central;
    rate_elim_FdUTP       = params.k_elim_FdUTP * C_FdUTP       * params.V_central;
    rate_elim_FUTP        = params.k_elim_FUTP  * C_FUTP        * params.V_central;
    rate_elim_FdUMP_tumor = params.k_elim_FdUMP * C_tumor_FdUMP * V_tumor_L;
    rate_elim_FdUTP_tumor = params.k_elim_FdUTP * C_tumor_FdUTP * V_tumor_L;
    rate_elim_FUTP_tumor  = params.k_elim_FUTP  * C_tumor_FUTP  * V_tumor_L;

    % --- renal excretion (CL_renal in mL/min; the /V*V cancels to CL*C) ---
    sum_metabolites_sys      = C_FdUMP + C_FdUTP + C_FUTP;
    rate_excrete_metabolites = (params.CL_renal_metabolites / 1000) * sum_metabolites_sys;
    rate_excrete_5FU         = (params.CL_renal_5FU  / 1000) * C_central;
    rate_excrete_FBAL        = (params.CL_renal_FBAL / 1000) * C_FBAL;

    % --- catabolite chain 5-FU -> DHFU -> FBAL ---
    rate_form_DHFU = rate_DPD_liver + rate_DPD_systemic;
    rate_form_FBAL = params.k_DHFU_to_FBAL * C_DHFU * params.V_central;

    if sum_metabolites_sys > 0
        rate_excrete_FdUMP = rate_excrete_metabolites * (C_FdUMP / sum_metabolites_sys);
        rate_excrete_FdUTP = rate_excrete_metabolites * (C_FdUTP / sum_metabolites_sys);
        rate_excrete_FUTP  = rate_excrete_metabolites * (C_FUTP  / sum_metabolites_sys);
    else
        rate_excrete_FdUMP = 0; rate_excrete_FdUTP = 0; rate_excrete_FUTP = 0;
    end
    rate_clear_metabolites = rate_excrete_metabolites + rate_elim_FdUMP + rate_elim_FdUTP + rate_elim_FUTP + ...
                             rate_elim_FdUMP_tumor + rate_elim_FdUTP_tumor + rate_elim_FUTP_tumor;

    % --- derivatives ---
    Vc = params.V_central;
    dC_central = dosingRate/Vc ...
        + rate_periph_to_central/Vc - rate_central_to_periph/Vc ...
        + rate_liver_to_central/Vc  - rate_central_to_liver/Vc ...
        + rate_tumor_to_central/Vc  - rate_central_to_tumor/Vc ...
        - rate_DPD_systemic/Vc ...
        - rate_form_FdUMP_sys/Vc - rate_form_FdUTP_sys/Vc - rate_form_FUTP_sys/Vc ...
        - rate_excrete_5FU/Vc;
    dC_peripheral = rate_central_to_periph/params.V_peripheral - rate_periph_to_central/params.V_peripheral;
    dC_liver      = (rate_central_to_liver - rate_liver_to_central - rate_DPD_liver) / V_liver_L;
    dC_tumor      = (rate_central_to_tumor - rate_tumor_to_central ...
                     - rate_form_FdUMP_tumor - rate_form_FdUTP_tumor - rate_form_FUTP_tumor) / V_tumor_L;

    dC_FdUMP = (rate_form_FdUMP_sys - rate_elim_FdUMP - rate_excrete_FdUMP) / Vc;
    dC_FdUTP = (rate_form_FdUTP_sys - rate_elim_FdUTP - rate_excrete_FdUTP) / Vc;
    dC_FUTP  = (rate_form_FUTP_sys  - rate_elim_FUTP  - rate_excrete_FUTP)  / Vc;

    dC_tumor_FdUMP = (rate_form_FdUMP_tumor - rate_elim_FdUMP_tumor) / V_tumor_L;
    dC_tumor_FdUTP = (rate_form_FdUTP_tumor - rate_elim_FdUTP_tumor) / V_tumor_L;
    dC_tumor_FUTP  = (rate_form_FUTP_tumor  - rate_elim_FUTP_tumor)  / V_tumor_L;

    dC_DHFU = (rate_form_DHFU - rate_form_FBAL) / Vc;
    dC_FBAL = (rate_form_FBAL - rate_excrete_FBAL) / Vc;

    d_excreted_5FU         = rate_excrete_5FU;
    d_excreted_FBAL        = rate_excrete_FBAL;
    d_excreted_metabolites = rate_clear_metabolites;
    d_cumulative_input     = dosingRate;

    % --- pharmacodynamics ---
    TS_inhib = C_tumor_FdUMP / (params.IC50_TS_tumor + C_tumor_FdUMP);
    d_dNTP   = updateDNTPPool(dNTP, TS_inhib, params);
    d_DNA    = updateDNADamage(DNA, dNTP, params);

    dydt = [dC_central; dC_peripheral; dC_liver; dC_tumor; ...
            dC_FdUMP; dC_FdUTP; dC_FUTP; ...
            dC_tumor_FdUMP; dC_tumor_FdUTP; dC_tumor_FUTP; ...
            dC_DHFU; dC_FBAL; ...
            d_excreted_5FU; d_excreted_FBAL; d_excreted_metabolites; d_cumulative_input; ...
            d_dNTP; d_DNA];
end
