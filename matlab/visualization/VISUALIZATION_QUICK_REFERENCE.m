% =========================================================================
% QUICK REFERENCE: Advanced Pharmacometric Visualizations
% =========================================================================
%
% This file provides copy-paste examples for common visualization tasks
% using the enhanced MC_5FU_PK_sensitivity.m functions.
%
% =========================================================================

%% BASIC USAGE - Automatic Visualizations
% =========================================================================
% Simply run the main function - all plots are automatically generated!

MC_5FU_PK_sensitivity(100, './my_results/')

% This creates:
%   - CDF plots with PTA markers (AUC and Cmax)
%   - Visual Predictive Checks (VPC)
%   - Box-and-whisker plots
%   - Probability of Target Attainment (PTA) plots

%% LOAD EXISTING RESULTS
% =========================================================================

% Load previously saved results
load('MC_results/2026-01-07_14-30-00/MC_5FU_PK_sensitivity_20260107.mat');

% Extract data
AUC = results.AUC;
Cmax = results.Cmax;
timeseries = results.raw_timeseries;

%% CDF PLOTS WITH PTA MARKERS
% =========================================================================

% Example 1: AUC with Gamelin 2008 targets (standard)
targets_auc = get5FU_AUC_targets('gamelin2008');
options.metric_name = 'AUC';
options.units = 'mg·h/L';
options.title_text = '5-FU AUC Distribution - Standard Targets';
plotCDF_with_PTA(AUC, targets_auc, options);

% Example 2: Cmax with Thyss 1986 targets (classic toxicity thresholds)
targets_cmax = get5FU_Cmax_targets('thyss1986');
options_cmax.metric_name = 'Cmax';
options_cmax.units = 'µM';
options_cmax.title_text = 'Peak Concentration with Toxicity Thresholds';
plotCDF_with_PTA(Cmax, targets_cmax, options_cmax);

% Example 3: Conservative dosing (for elderly/frail patients)
targets_conservative = get5FU_AUC_targets('conservative');
options_cons.metric_name = 'AUC';
options_cons.units = 'mg·h/L';
options_cons.title_text = 'Conservative Dosing Strategy';
options_cons.show_percentiles = true;
plotCDF_with_PTA(AUC, targets_conservative, options_cons);

% Example 4: Aggressive dosing (for fit patients)
targets_aggressive = get5FU_AUC_targets('aggressive');
options_agg.metric_name = 'AUC';
options_agg.units = 'mg·h/L';
options_agg.title_text = 'Aggressive Dosing - Select Patients';
plotCDF_with_PTA(AUC, targets_aggressive, options_agg);

%% VISUAL PREDICTIVE CHECKS
% =========================================================================

% Example 1: Standard VPC with linear scale
time_vec = timeseries{1}.time_hr;
n_time = length(time_vec);
n_sims = length(timeseries);
conc_matrix = zeros(n_time, n_sims);
for i = 1:n_sims
    conc_matrix(:, i) = timeseries{i}.C_central;
end
observed = median(conc_matrix, 2);  % Or use actual patient data

vpc_options.title_text = '5-FU Central Compartment VPC';
vpc_options.xlabel_text = 'Time (hours)';
vpc_options.ylabel_text = 'Concentration (µM)';
vpc_options.log_scale = false;
vpc_options.show_individual_sims = true;
plotVPC(observed, conc_matrix, time_vec, vpc_options);

% Example 2: VPC with log scale (for wide concentration ranges)
vpc_log.title_text = '5-FU VPC (Log Scale)';
vpc_log.xlabel_text = 'Time (hours)';
vpc_log.ylabel_text = 'Concentration (µM)';
vpc_log.log_scale = true;  % Use log scale
vpc_log.show_individual_sims = false;  % Cleaner plot
vpc_log.percentiles = [5 50 95];
plotVPC(observed, conc_matrix, time_vec, vpc_log);

%% BOX-AND-WHISKER PLOTS
% =========================================================================

% Example 1: Compare exposure quartiles
[~, sort_idx] = sort(AUC);
n_quartiles = 4;
quartile_size = floor(length(AUC) / n_quartiles);

box_data.AUC = cell(n_quartiles, 1);
box_data.Cmax = cell(n_quartiles, 1);
box_data.labels = cell(n_quartiles, 1);

for q = 1:n_quartiles
    idx_start = (q-1) * quartile_size + 1;
    idx_end = min(q * quartile_size, length(AUC));
    indices = sort_idx(idx_start:idx_end);
    
    box_data.AUC{q} = AUC(indices);
    box_data.Cmax{q} = Cmax(indices);
    box_data.labels{q} = sprintf('Q%d', q);
end

box_options.show_points = true;
box_options.point_jitter = 0.2;
box_options.show_mean = true;
box_options.title_text = 'PK Variability Across Exposure Quartiles';
plotBoxWhisker_PK(box_data, box_options);

% Example 2: Compare two dose groups (with statistics)
% Assuming you have data from two dosing regimens:
dose_data.AUC = {auc_low_dose, auc_high_dose};
dose_data.Cmax = {cmax_low_dose, cmax_high_dose};
dose_data.labels = {'400 mg/m²', '800 mg/m²'};

dose_box_opts.show_points = true;
dose_box_opts.add_stats = true;  % Adds p-value comparison
dose_box_opts.title_text = 'Dose Comparison Study';
plotBoxWhisker_PK(dose_data, dose_box_opts);

% Example 3: Compare three patient populations
pop_data.AUC = {auc_young, auc_middle, auc_elderly};
pop_data.Cmax = {cmax_young, cmax_middle, cmax_elderly};
pop_data.labels = {'Age <50', 'Age 50-65', 'Age >65'};

pop_box_opts.show_points = true;
pop_box_opts.show_mean = true;
pop_box_opts.point_alpha = 0.4;
pop_box_opts.title_text = 'Age-Stratified PK Analysis';
pop_box_opts.box_colors = {[0.4 0.76 0.65], [0.99 0.55 0.38], [0.55 0.63 0.80]};
plotBoxWhisker_PK(pop_data, pop_box_opts);

%% PROBABILITY OF TARGET ATTAINMENT (PTA)
% =========================================================================

% Example 1: Standard PTA plot for AUC
pta_efficacy.values = [20, 28];
pta_efficacy.labels = {'Minimum efficacy', 'Optimal target'};

pta_toxicity.values = [50, 70];
pta_toxicity.labels = {'Grade 3+ risk (20%)', 'Grade 3+ risk (60%)'};

pta_options.metric_name = 'AUC';
pta_options.units = 'mg·h/L';
pta_options.title_text = 'Dose Optimization: Efficacy vs Toxicity';
pta_options.show_optimal = true;
pta_options.optimal_pta = 90;      % Require ≥90% PTA
pta_options.max_tox_prob = 10;     % Accept ≤10% toxicity
plotPTA(AUC, pta_efficacy, pta_toxicity, pta_options);

% Example 2: PTA for Cmax (acute toxicity assessment)
cmax_efficacy.values = [300, 500];
cmax_efficacy.labels = {'Minimum peak', 'Target peak'};

cmax_toxicity.values = [800, 1200];
cmax_toxicity.labels = {'Moderate toxicity', 'Severe toxicity'};

cmax_pta_opts.metric_name = 'Cmax';
cmax_pta_opts.units = 'µM';
cmax_pta_opts.title_text = 'Peak Concentration: Acute Toxicity Risk';
cmax_pta_opts.show_optimal = true;
cmax_pta_opts.optimal_pta = 85;
cmax_pta_opts.max_tox_prob = 15;
plotPTA(Cmax, cmax_efficacy, cmax_toxicity, cmax_pta_opts);

% Example 3: Strict PTA requirements (high-risk patients)
strict_pta.metric_name = 'AUC';
strict_pta.units = 'mg·h/L';
strict_pta.title_text = 'High-Risk Population: Strict Safety Requirements';
strict_pta.show_optimal = true;
strict_pta.optimal_pta = 95;       % Very high PTA required
strict_pta.max_tox_prob = 5;       % Very low toxicity tolerance
plotPTA(AUC, pta_efficacy, pta_toxicity, strict_pta);

%% USING PREDEFINED TARGET SETS
% =========================================================================

% Available AUC target sets:
targets_gamelin = get5FU_AUC_targets('gamelin2008');    % Standard
targets_kaldate = get5FU_AUC_targets('kaldate2012');    % PopPK model
targets_beumer = get5FU_AUC_targets('beumer2019');      % TDM guidelines
targets_cons = get5FU_AUC_targets('conservative');      % Safety-focused
targets_agg = get5FU_AUC_targets('aggressive');         % Efficacy-focused

% Available Cmax target sets:
targets_thyss = get5FU_Cmax_targets('thyss1986');       % Classic
targets_saif = get5FU_Cmax_targets('saif2013');         % Modern
targets_cons_cm = get5FU_Cmax_targets('conservative');  % Safety
targets_agg_cm = get5FU_Cmax_targets('aggressive');     % Efficacy

% Each target struct contains:
%   .values      - Numeric thresholds
%   .labels      - Description for each threshold
%   .colors      - RGB colors for plotting
%   .reference   - Literature citation

%% CUSTOMIZATION: EXPORT/IMPORT CONFIGURATIONS
% =========================================================================

% Export default configuration template
exportPlotOptions('my_custom_config.mat');

% Edit the file manually, then load
load('my_custom_config.mat');

% Use custom configuration
plotCDF_with_PTA(AUC, targets_auc, plot_config.cdf);
plotVPC(observed, conc_matrix, time_vec, plot_config.vpc);
plotBoxWhisker_PK(box_data, plot_config.box);
plotPTA(AUC, pta_efficacy, pta_toxicity, plot_config.pta);

%% CUSTOM TARGETS (Not using predefined sets)
% =========================================================================

% Define your own targets
custom_targets.values = [15, 25, 35, 45];
custom_targets.labels = {'Low', 'Target', 'Upper', 'Toxic'};
custom_targets.colors = {'yellow', 'green', 'orange', 'red'};

custom_opts.metric_name = 'Custom Exposure';
custom_opts.units = 'arbitrary units';
custom_opts.title_text = 'Custom Analysis';
plotCDF_with_PTA(my_data, custom_targets, custom_opts);

%% SAVING FIGURES
% =========================================================================

% PDF (vector graphics, publication-ready)
print(gcf, 'my_figure.pdf', '-dpdf', '-bestfit');

% High-resolution PNG
print(gcf, 'my_figure.png', '-dpng', '-r300');

% TIFF for journals
print(gcf, 'my_figure.tiff', '-dtiff', '-r600');

% EPS for LaTeX
print(gcf, 'my_figure.eps', '-depsc');

% MATLAB .fig for future editing
savefig('my_figure.fig');

%% BATCH PROCESSING
% =========================================================================

% Create all plots for multiple analyses
analysis_dirs = {'analysis1', 'analysis2', 'analysis3'};

for i = 1:length(analysis_dirs)
    % Load results
    result_file = fullfile(analysis_dirs{i}, 'MC_5FU_PK_sensitivity.mat');
    load(result_file);
    
    % Create output directory
    output_dir = fullfile(analysis_dirs{i}, 'plots');
    mkdir(output_dir);
    
    % Generate all plots
    targets = get5FU_AUC_targets('gamelin2008');
    plotCDF_with_PTA(results.AUC, targets, struct('metric_name','AUC','units','mg·h/L'));
    print(gcf, fullfile(output_dir, 'CDF_AUC.pdf'), '-dpdf');
    
    % ... repeat for other plot types
end

%% TROUBLESHOOTING
% =========================================================================

% Check for NaN/Inf values
sum(isnan(AUC))
sum(isinf(AUC))

% Remove invalid data
AUC_clean = AUC(isfinite(AUC));
Cmax_clean = Cmax(isfinite(Cmax));

% Check data range
fprintf('AUC range: %.2f - %.2f mg·h/L\n', min(AUC), max(AUC));
fprintf('Cmax range: %.2f - %.2f µM\n', min(Cmax), max(Cmax));

% Verify timeseries structure
if ~isempty(timeseries{1})
    disp('Timeseries fields:');
    disp(fieldnames(timeseries{1}));
end

%% GETTING HELP
% =========================================================================

% View function documentation
help plotCDF_with_PTA
help plotVPC
help plotBoxWhisker_PK
help plotPTA
help get5FU_AUC_targets
help get5FU_Cmax_targets
help exportPlotOptions

% View main function documentation
help MC_5FU_PK_sensitivity

% Check available target sets
get5FU_AUC_targets('help')     % Lists available sets
get5FU_Cmax_targets('help')    % Lists available sets

% =========================================================================
% END OF QUICK REFERENCE
% =========================================================================
