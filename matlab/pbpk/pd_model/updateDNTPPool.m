function dNTP_rate = updateDNTPPool(dNTP_pool_fraction, TS_inhibition_fraction, params)
% Rate of change for dNTP pool fraction: recovery minus TS-inhibition-driven depletion.
    dNTP_rate = params.k_repletion * (1 - dNTP_pool_fraction) ...
              - params.k_depletion * TS_inhibition_fraction * dNTP_pool_fraction;
end
