cd('c:/Users/hp/OneDrive - University of Southampton/Year 3/4 IP/PBPK-ABM-thesis');
try
    s = evalc("r = run5FU_PBPK_Simulation('Mayo.csv','diag_mayo_silent');"); %#ok<NASGU>
    fprintf('AUC=%.6f\n', r.metrics.AUC_central_mg_h_L);
    fprintf('TumorAUC=%.6f\n', r.metrics.AUC_tumor_mg_h_L);
catch ME
    disp(getReport(ME,'extended'));
    exit(1);
end
exit(0);
