function damage_rate = updateDNADamage(DNA_damage_index, dNTP_pool_fraction, params)
% Rate of change for DNA damage index: S-phase accumulation minus saturable repair.
    repair = (params.Vmax_repair * DNA_damage_index) / (params.Km_repair + DNA_damage_index);
    damage_rate = params.k_damage * (1 - dNTP_pool_fraction) * params.S_phase_fraction - repair;
end
