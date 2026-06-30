function ts_inhib = calculateTSInhibition(C_tumor_FdUMP, params)
% Emax (Hill n=1) fraction of TS bound by FdUMP in tumor tissue.
    ts_inhib = C_tumor_FdUMP / (params.IC50_TS_tumor + C_tumor_FdUMP);
end
