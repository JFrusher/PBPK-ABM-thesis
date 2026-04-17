function params = initialise5FUParameters()
% Initialise a default parameter set for 5-FU PBPK simulation
% This is a standalone copy extracted to be discoverable by unit tests.

fprintf('Initialising physiological parameters from literature...\n');

% Basic anthropometrics
params.BW = 70;              % Body weight (kg)
params.height_cm = 175;      % Default height (cm) - used to compute BSA
params.BSA = calculateBSA(params.BW, params.height_cm); % Body surface area (m^2)
params.CO = 6.5;             % Cardiac output (L/min) - typical resting value

% Volumes (fractions or L where appropriate)
params.V_blood = 0.0743;
params.V_liver = 0.026;
params.V_kidney = 0.004;
params.V_heart = 0.005;
params.V_brain = 0.02;
params.V_muscle = 0.4;
params.V_fat = 0.2;
params.V_skin = 0.05;
params.V_tumor = 0.0005;

% Blood flow fractions
params.Q_liver = 0.25;
params.Q_kidney = 0.19;
params.Q_heart = 0.04;
params.Q_brain = 0.12;
params.Q_muscle = 0.17;
params.Q_fat = 0.05;
params.Q_skin = 0.05;
params.Q_tumor = 0.02;

% Distribution volumes
params.V_central = params.BW * 0.2;
params.V_peripheral = params.BW * 0.3;
params.Kp_liver = 1.2;
params.Kp_kidney = 1.1;
params.Kp_muscle = 0.8;
params.Kp_fat = 0.3;
params.Kp_tumor = 0.61;
params.Kp_brain = 0.4;
params.k_cp = 0.15;
params.k_pc = 0.08;

% DPD kinetics (validated values)
MW_5FU = 130.08; % g/mol
params.Vmax_DPD_mg_per_h = 1220;
params.Vmax_DPD = (params.Vmax_DPD_mg_per_h / MW_5FU) / 60 * 1000; % µmol/min
params.Km_DPD = 5.0; % mg/L

% Diagnostics
fprintf('\nDPD KINETIC PARAMETERS (Perfusion-Limited Model)\n');
fprintf('Vmax_DPD: %.3f µmol/min = %.1f mg/h\n', params.Vmax_DPD, params.Vmax_DPD_mg_per_h);
fprintf('Km_DPD: %.1f mg/L\n', params.Km_DPD);

% Circadian
params.DPD_mean = 1.0;
params.DPD_amplitude = 0.37;
params.DPD_acrophase = 1; % 1 AM

% Some additional defaults used elsewhere
params.MW_5FU = MW_5FU;
params.k_form_FdUMP = 0.01;
params.k_form_FdUTP = 0.01;
params.k_form_FUTP = 0.01;
params.k_elim_FdUMP = 0.05;
params.k_elim_FdUTP = 0.05;
params.k_elim_FUTP = 0.05;
params.k_DHFU_to_FBAL = 0.1;
params.CL_renal_5FU = 0.05;
params.CL_renal_FBAL = 0.1;
params.CL_renal_metabolites = 0.05;

fprintf('Parameter initialisation complete.\n\n');
end