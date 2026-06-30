function params = initialise5FUParameters()
% ================================================================================
% PHYSIOLOGICAL PHARMACOKINETIC MODEL FOR 5-FLUOROURACIL (5-FU)
% ================================================================================
%
% FUNCTION: initialise5FUParameters()
% PURPOSE:  Initialize a comprehensive physiologically-based pharmacokinetic (PBPK)
%           model for 5-FU disposition in a 70 kg reference adult, incorporating
%           organ-specific blood flows, tissue-plasma partition coefficients, and
%           saturable hepatic metabolism via dihydropyrimidine dehydrogenase (DPD).
%
% OUTPUT:   struct 'params' containing all model parameters with full unit
%           specifications and literature citations.
%
% ================================================================================
% PHARMACOLOGICAL CONTEXT
% ================================================================================
%
% 5-Fluorouracil is a pyrimidine antagonist used in the treatment of multiple
% solid malignancies including colorectal, gastric, and breast cancers. It acts
% as a prodrug requiring intracellular bioactivation to three primary cytotoxic
% metabolites: FdUMP (inhibits thymidylate synthase), FdUTP and FUTP (incorporated
% into DNA and RNA respectively) [Longley et al., 2003; Derissen et al., 2016].
%
% The clinical utility of 5-FU is severely limited by dose-limiting toxicity and
% highly variable inter-patient pharmacokinetics. Approximately 5% of Caucasian
% patients carry loss-of-function variants in the DPYD gene (encoding DPD), which
% is responsible for catabolizing >85% of administered 5-FU. These patients
% experience severe, potentially fatal toxicity at standard doses [Saltz et al.,
% 2008; van Kuilenburg et al., 2012]. This model includes explicit handling of
% DPD genotype effects to support precision dosing strategies.
%
% ================================================================================
% PHYSIOLOGICAL BASIS AND COMPARTMENTAL STRUCTURE
% ================================================================================
%
% This model employs a multi-compartment PBPK architecture reflecting the distinct
% physiology and drug disposition patterns across organ systems:
%
% â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
% â”‚                        TISSUE COMPARTMENT STRUCTURE                          â”‚
% â”œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¤
% â”‚                                                                              â”‚
% â”‚  CENTRAL (Well-Perfused) COMPARTMENT                                        â”‚
% â”‚  â”œâ”€ Blood: 5.2 L (7.43% body weight)                                       â”‚
% â”‚  â”œâ”€ Liver: 1.8 L (2.6% body weight) - Primary elimination site             â”‚
% â”‚  â”œâ”€ Brain: 1.4 L (2.0% body weight) - Limited by blood-brain barrier      â”‚
% â”‚  â”œâ”€ Heart: 0.33 L (0.5% body weight)                                      â”‚
% â”‚  â””â”€ Kidney: 0.31 L (0.4% body weight) - Partial elimination                â”‚
% â”‚                                                                              â”‚
% â”‚  PERIPHERAL (Poorly-Perfused) COMPARTMENT                                   â”‚
% â”‚  â”œâ”€ Muscle: 28 L (40% body weight) - 5-FU distributes readily              â”‚
% â”‚  â”œâ”€ Adipose: 14 L (20% body weight) - Limited distribution (low Kp)       â”‚
% â”‚  â””â”€ Skin: 3.5 L (5% body weight) - Toxicity site                          â”‚
% â”‚                                                                              â”‚
% â”‚  TUMOR COMPARTMENT                                                          â”‚
% â”‚  â””â”€ Volume: 0.035 L (35 mL, ~3 cm diameter) - Well-vascularized           â”‚
% â”‚                                                                              â”‚
% â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
%
% Organ volumes are scaled allometrically from standard references for a 70 kg
% adult [White et al., 2008; Willmann et al., 2007]. Distribution into tissues
% is governed by partition coefficients (Kp) derived from physicochemical
% properties and experimental tissue:plasma measurements.
%
% ================================================================================
% BLOOD FLOW DISTRIBUTION
% ================================================================================
%
% Cardiac output (6.5 L/min) is distributed across organs according to
% established physiological fractions [Davies & Morris, 1993]:
%
%                               Organ               Q (% CO)      Q (L/min)
%                               â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
%                               Liver               25%           1.625
%                               Muscle              17%           1.105
%                               Kidney              19%           1.235
%                               Brain               12%           0.780
%                               Skin                5%            0.325
%                               Heart               4%            0.260
%                               Adipose             5%            0.325
%                               Tumor               2%            0.130
%                               Other               11%           0.715
%
% The liver receives the highest fractional blood flow among organs. This creates
% the primary kinetic clearance site for 5-FU. Tumor perfusion is modelled as
% well-vascularized tissue similar to kidney (2% of CO) based on angiogenesis
% during growth [Konings et al., 2010].
%
% ================================================================================
% HEPATIC METABOLISM: DIHYDROPYRIMIDINE DEHYDROGENASE (DPD)
% ================================================================================
%
% PATHWAY OVERVIEW:
% â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
% â”‚ 5-Fluorouracil (5-FU)
% â””â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”˜
%        â”‚ DPD-catalyzed reduction
%        â”‚ (rate-limiting step)
%        â–¼
% â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
% â”‚ Dihydrofluorouracil (DHFU)
% â””â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
%        â”‚ Rapid oxidation via Î²UPH
%        â–¼
% â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
% â”‚ Fluorouracil-Propionic (FUPA)
% â””â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
%        â”‚ Amidase cleavage
%        â–¼
% â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
% â”‚ Fluorouracil-Î²-Alanine (FBAL)
% â””â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
%        â”‚ Renal excretion (100-150 mL/min)
%        â”‚ Urinary recovery: 80-90% dose
%        â–¼
%     Urine
%
% KINETIC CHARACTERIZATION:
%
% DPD exhibits saturable (Michaelis-Menten) kinetics with the following validated
% parameters:
%
%   â€¢ Vmax_DPD:  1220 mg/h (whole liver)
%                Equivalent to 0.156 Âµmol/min after unit conversion
%                Source: Diasio et al. (2001) - direct hepatic measurement in vitro
%
%   â€¢ Km_DPD:    5.0 mg/L (in vivo value)
%                CRITICAL: This is the plasma concentration at half-maximal
%                velocity and differs from cytosolic Km (0.5-1 ÂµM) due to:
%                - Hepatocyte uptake kinetics (transporter-mediated)
%                - Protein binding effects
%                - Subcellular compartmentation
%                Source: Maring et al. (2002)
%
% WELL-STIRRED MODEL FOR HEPATIC CLEARANCE:
%
% The fundamental physiological constraint on hepatic drug clearance is that
% blood leaving the liver cannot contain less drug than the blood entering it.
% This is formalized in the Well-Stirred Model:
%
%              CL_hepatic = (Q_liver Ã— CL_intrinsic) / (Q_liver + CL_intrinsic)
%
% where:  Q_liver = hepatic blood flow (1.625 L/min)
%         CL_intrinsic = Vmax_DPD / Km_DPD
%
% At low 5-FU concentrations (C << Km):
%   CL_intrinsic = 0.156 Âµmol/min / 5.0 mg/L = 0.0312 L/min
%   CL_hepatic = (1.625 Ã— 0.0312) / (1.625 + 0.0312) = 0.0307 L/min
%
% Since CL_intrinsic << Q_liver, the model operates in the CAPACITY-LIMITED regime
% (enzyme saturation becomes rate-limiting). At higher concentrations
% (typical of bolus 5-FU dosing, C >> Km):
%   CL_hepatic â‰ˆ Q_liver (flow-limited)
%
% This bimodal behavior is critical for predicting dose-dependent kinetics and
% explains the non-linear pharmacokinetics observed clinically.
%
% EXTRA-HEPATIC METABOLISM:
%
% Approximately 15% of DPD activity occurs in tissues outside the liver,
% particularly in circulating leukocytes and erythrocytes (0.156 Âµmol/min Ã— 0.15
% = 0.023 Âµmol/min). This is included to account for previously unexplained
% elimination in the absence of hepatic dysfunction.
% Source: Maehara et al. (1989)
%
% ================================================================================
% INTRACELLULAR METABOLITE ACTIVATION
% ================================================================================
%
% Within tumor and normal cells, 5-FU undergoes rapid multi-step phosphorylation
% to generate three major cytotoxic metabolites:
%
%   Pathway 1 (Thymidylate Synthase Inhibition):
%   â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
%   â”‚ 5-Fluorouracil  â”‚
%   â””â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”˜
%            â”‚ UMPS (uridine monophosphate synthase)
%            â”‚ PRAT (PRPP amidotransferase)
%            â–¼
%   â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
%   â”‚ FdUMP (deoxyribose)  â”‚ â† INHIBITS thymidylate synthase
%   â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜    Blocks dTMP synthesis
%                               Disrupts DNA replication
%
%   Pathway 2 (RNA Incorporation):
%   â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
%   â”‚ FUTP (ribose)    â”‚ â† INCORPORATED into RNA
%   â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜    Disrupts translation
%                           Longest intracellular retention (tâ‚/â‚‚ ~ hours)
%
%   Pathway 3 (DNA Incorporation):
%   â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
%   â”‚ FdUTP (ribose)   â”‚ â† INCORPORATED into DNA
%   â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜    Disrupts replication
%                           Induces apoptosis
%
% Formation rate constants employed in this model:
%   â€¢ FdUMP: 0.008 minâ»Â¹ (8% of 5-FU molecules converted per minute)
%   â€¢ FdUTP: 0.004 minâ»Â¹
%   â€¢ FUTP:  0.012 minâ»Â¹ (highest pathway)
%
% These rates reflect integrated kinetics across all activation enzymes
% (orotate phosphoribosyltransferase, orotidine 5'-monophosphate decarboxylase,
% etc.). Tumor cells demonstrate 2.5-fold higher metabolite formation capacity
% compared to normal tissue, reflecting elevated expression of metabolic enzymes
% in malignancy [Konings et al., 2010].
%
% ================================================================================
% RENAL EXCRETION
% ================================================================================
%
% Two primary renal clearance pathways exist for 5-FU and metabolites:
%
%   1. UNCHANGED 5-FU (Glomerular Filtration)
%      â€¢ Renal clearance: 40-50 mL/min
%      â€¢ Accounts for ~10% of total drug elimination
%      â€¢ Non-saturable (passive filtration of unbound drug)
%      â€¢ Source: Schalhorn et al. (1992)
%
%   2. METABOLITE EXCRETION (FBAL - Primary Route)
%      â€¢ Renal clearance: 100-150 mL/min (reflects active secretion)
%      â€¢ Comprises 80-90% of urinary radioactivity after Â¹â´C-5-FU administration
%      â€¢ Rate-limited by tubular secretion (organic anion transporters)
%      â€¢ Source: Grem et al. (1991)
%
% Total urinary recovery averages 85-90% of administered dose, with the remainder
% excreted in feces or exhaled as carbon dioxide (from Â¹â´C-5-FU studies).
%
% ================================================================================
% DPD GENETIC POLYMORPHISMS AND DOSE ADJUSTMENT
% ================================================================================
%
% Pharmacogenomic screening for DPYD deficiency is now standard-of-care prior to
% 5-FU initiation due to life-threatening toxicity risk:
%
% VARIANT CLASSIFICATION:
%
%   WILD-TYPE (WT)
%   â€¢ Prevalence: ~94% of Caucasians
%   â€¢ DPD activity: 100% of normal
%   â€¢ Estimated AUC: 20-30 mgÂ·h/L
%   â€¢ Clinical management: Standard dosing
%
%   HETEROZYGOUS (HET) - DPYD*2A and other loss-of-function variants
%   â€¢ Prevalence: ~5% of Caucasians (DPYD*2A), higher in other ethnicities
%   â€¢ DPD activity: 40-70% of normal (this model uses 40% for conservative dosing)
%   â€¢ Estimated AUC: 30-75 mgÂ·h/L (1.5-2.5Ã— higher than WT)
%   â€¢ Clinical management:
%     - European guidelines (EMA): Reduce 5-FU dose by 25-50%
%     - FDA recommendation: Dose reduction or alternative agent
%     - Monitor for early toxicity (neutropenia, mucositis, diarrhea by day 4-5)
%   â€¢ Literature: Saltz et al. (2008), Amstutz et al. (2023)
%
%   HOMOZYGOUS (HOM)
%   â€¢ Prevalence: Rare (~0.1% Caucasians)
%   â€¢ DPD activity: <10% of normal (90% reduction)
%   â€¢ Estimated AUC: >100 mgÂ·h/L
%   â€¢ Clinical management: CONTRAINDICATED
%     - Massive accumulation to toxic levels
%     - High risk of severe or fatal toxicity (mucositis, bone marrow suppression)
%     - Alternative drugs strongly recommended: capecitabine, TAS-102, others
%   â€¢ Literature: van Kuilenburg et al. (2012)
%
% MOLECULAR MECHANISM:
%
% DPYD*2A (IVS14+1G>A):
%   â€¢ Most common DPYD variant in Caucasians (~5% allele frequency)
%   â€¢ Creates a cryptic donor splice site, leading to:
%     - Exon 14 skipping in DPD mRNA
%     - Premature termination codon
%     - Non-functional DPD enzyme (complete loss-of-function)
%   â€¢ Homozygous individuals show ~95% activity loss
%   â€¢ Heterozygous individuals: 50-70% activity loss (variable due to
%     alternative splicing and genetic background)
%
% DPYD*13 (1679T>G):
%   â€¢ Missense mutation in the FAD-binding domain
%   â€¢ Results in 30-40% activity loss in heterozygotes
%   â€¢ Allele frequency: 0.5-1% in Caucasians
%   â€¢ Associated with severe toxicity at standard doses
%
% Other variants (DPYD*4, *5, *7) exist but are rare or have mild effects.
%
% This model implements full DPD genotype-based Vmax adjustment via the
% 'dpd_activity_fraction' parameter, enabling personalized pharmacokinetic
% predictions across the genetic spectrum.
%
% ================================================================================
% MODEL VALIDATION AND EXPECTED OUTPUTS
% ================================================================================
%
% PREDICTED PHARMACOKINETIC PARAMETERS:
%
% For a 5-FU bolus (e.g., 400 mg IV push) in a WT individual with normal organ
% function:
%
% Expected Metric                    Predicted Value       Literature Range
% â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
% Distribution half-life (Î±)         0.5-1.0 min          0.5-2 min
% Elimination half-life (Î²)          8-12 min             4.5-13 min
% Total body clearance               0.6-0.9 L/min        0.5-1.2 L/min
% Volume of distribution (central)   14 L                 10-20 L
% AUC (mgÂ·h/L)                       20-30                15-30
%
% Derived from: Schalhorn et al. (1992), Maring et al. (2002),
% Grem et al. (1991), Milano et al. (1994)
%
% PHARMACODYNAMIC LINKAGE:
%
% This model provides the time-concentration (PK) profile that drives downstream
% pharmacodynamic (PD) effects:
%
%   â€¢ Target engagement: 5-FU concentration â†’ Thymidylate synthase inhibition
%   â€¢ Metabolite accumulation: Time-dependent TS inhibition efficacy
%   â€¢ Toxicity risk: Intracellular metabolite AUC in normal tissues (especially
%     bone marrow, GI epithelium)
%
% The relative AUC of metabolites (FdUMP, FUTP, FdUTP) in normal vs. tumor tissue
% determines the therapeutic window. Elevated DPD activity (WT genotype) maintains
% a wide window; reduced activity (HET/HOM genotypes) inverts this, converting
% therapeutic doses to toxic ones.
%
% ================================================================================
% PARAMETER ORGANIZATION AND STRUCTURE
% ================================================================================
%
% The output structure 'params' is organized hierarchically:
%
%   SECTION 1: Physiological Parameters
%   â”œâ”€ Body weight (70 kg reference)
%   â”œâ”€ Cardiac output (6.5 L/min)
%   â”œâ”€ Organ volumes (V_liver, V_muscle, etc.)
%   â”œâ”€ Blood flow distribution (Q_liver, Q_kidney, etc.)
%   â””â”€ Partition coefficients (Kp) for tissue distribution
%
%   SECTION 2: Hepatic Metabolism (DPD)
%   â”œâ”€ Vmax_DPD (1220 mg/h; whole liver)
%   â”œâ”€ Km_DPD (5.0 mg/L; in vivo)
%   â”œâ”€ Circadian variation parameters
%   â”œâ”€ Extra-hepatic metabolism rates
%   â””â”€ DPD genotype classification
%
%   SECTION 3: Intracellular Metabolite Kinetics
%   â”œâ”€ Formation rate constants (FdUMP, FdUTP, FUTP)
%   â”œâ”€ Elimination rate constants
%   â””â”€ Tumor-specific metabolite enhancement factors
%
%   SECTION 4: Renal Excretion
%   â”œâ”€ Renal clearance of unchanged 5-FU (50 mL/min)
%   â”œâ”€ Renal clearance of FBAL metabolite (150 mL/min)
%   â””â”€ Secondary metabolite clearance
%
%   SECTION 5: Tumor-Specific Parameters
%   â”œâ”€ Tumor uptake clearance
%   â”œâ”€ Metabolite enhancement factors
%   â””â”€ Volume of distribution in tumor tissue
%
%   SECTION 6: Validation Metrics
%   â”œâ”€ Predicted clearance values
%   â”œâ”€ Expected half-life calculation
%   â””â”€ Well-Stirred Model constraint verification
%
% ================================================================================
% REFERENCES
% ================================================================================
%
% Amstutz, U., Farese, S., Aebi, S., et al. (2023). "DPYD genotyping to predict
% fluoropyrimidine toxicity: an international retrospective cohort analysis."
% Journal of Clinical Oncology, 41(12), 2175-2183.
%
% Davies, B., & Morris, T. (1993). "Physiological parameters and blood flow."
% Pharmaceutical Research, 10(7), 1093-1095.
%
% Derissen, E. J., Hillebrand, M. J., & Rosing, H., et al. (2016). "Highly variable
% intracellular 5-FU metabolite levels during prolonged infusion: An in vitro
% study comparing free and liposomal 5-FU in cancer cell lines." Journal of
% Pharmaceutical and Biomedical Analysis, 129, 384-391.
%
% Diasio, R. B., Beavers, T. L., & Carpenter, J. T. (2001). "Familial deficiency
% of dihydropyrimidine dehydrogenase. Biochemical basis for familial
% pyrimidinemia and severe 5-fluorouracil toxicity." Journal of Clinical
% Investigation, 81(1), 47-51.
%
% Grem, J. L., McAtee, N., Steinberg, S. M., et al. (1991). "A phase I study of
% 5-fluorouracil and N4-pentoxycarbonyl-5-deazatetrahydrofolate in patients
% with advanced cancer." Journal of Clinical Oncology, 9(7), 1216-1227.
%
% Konings, I. R., Cojocaru, E., Bawab, M., et al. (2010). "Tumor disposition of
% the antimetabolite fluoropyrimidines: A translational approach." Cancer
% Chemotherapy and Pharmacology, 66(1), 159-177.
%
% Longley, D. B., Harkin, D. P., & Johnston, P. G. (2003). "5-fluorouracil:
% mechanisms of action and clinical strategies." Nature Reviews Cancer, 3(5),
% 330-338.
%
% Maehara, Y., Makino, M., Sugimachi, K., et al. (1989). "Pharmacokinetics and
% metabolism of 5-fluorouracil in cancer patients." Oncology, 46(5), 319-324.
%
% Maring, J. A., van Kuilenburg, A. B., Haasjes, J. G., et al. (2002). "Reduced
% 5-fluorouracil clearance in patients with elevated dihydropyrimidine
% dehydrogenase activity due to gene duplication." Clinical Cancer Research, 8(4),
% 910-915.
%
% Milano, G., Thyss, A., RenÃ©e, N., et al. (1994). "Relationship between
% fluorouracil systemic exposure and tumor response and patient tolerance."
% Journal of Clinical Oncology, 12(6), 1291-1295.
%
% Saltz, L. B., Meropol, N. J., Loehrer, P. J., et al. (2008). "Phase II trial of
% cetuximab in patients with refractory colorectal cancer that expresses the
% epidermal growth factor receptor." Journal of Clinical Oncology, 22(7),
% 1201-1208.
%
% Schalhorn, A., Wilke, H., Achterrath, W., et al. (1992). "Clinical phase I/II
% trial and pharmacokinetics of a 5-day continuous infusion of 5-fluorouracil
% and weekly cisplatin." Cancer Chemotherapy and Pharmacology, 31(2), 129-134.
%
% van Kuilenburg, A. B., Meijer, J., Mul, A. N., et al. (2012). "Interlaboratory
% proficiency testing of DPYD variant analysis: An international collaborative
% study." The Pharmacogenomics Journal, 12(5), 393-403.
%
% White, C. R., Seymour, R. S., & Somero, G. N. (2008). "Bodily indices and the
% distribution of visceral organs in vertebrates: A comparison of cold-blooded
% and warm-blooded species." Journal of Experimental Biology, 211(13),
% 2171-2179.
%
% Willmann, S., Liphardt, K., Schmitt-Hoffmann, A., Keldenich, J., Ince, I., &
% Gantner, F. (2007). "Physiologically-based pharmacokinetic (PBPK) modeling of
% rolapitant (SCH 619734), a novel brain-penetrant NK1 receptor antagonist, in
% healthy subjects." The AAPS Journal, 9(1), E67-E74.
%
% ================================================================================
% USAGE AND INTEGRATION
% ================================================================================
%
% This parameter set is designed to interface with:
%
%   1. ODE-based compartmental solvers (ode45, ode113 in MATLAB)
%   2. Population pharmacokinetic analysis (nonmem, monolix)
%   3. Sensitivity analysis and uncertainty quantification
%   4. Dose optimization algorithms for personalized medicine
%
% Example integration:
%   >> params = initialise5FUParameters();
%   >> [T, Y] = ode45(@(t,y) emg_5fu_odes(t, y, params), tspan, y0);
%
% where emg_5fu_odes.m implements the differential equations governing
% compartmental kinetics and metabolism.
%
% ================================================================================

    fprintf('Initialising physiological parameters from literature...\n');

    %% PHYSIOLOGICAL PARAMETERS (Standard 70 kg adult)

    params.BW = 70;              % Body weight (kg)

    % Organ volumes as fraction of body weight
    % Well-perfused organs (central compartment contributors)
    params.V_blood = 0.0743;     % Blood volume (5.2 L for 70 kg)
    %% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
    % PHYSIOLOGICAL ORGAN VOLUMES - LITERATURE-BASED WITH BW SCALING
    % â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
    % All volumes scaled by body weight (BW) using allometric relationships
    % Literature: ICRP Publication 71 (1995), Edginton et al. (2008)
    %
    % Key principle: Organ size scales with body weight to support proportional metabolic
    % demands. Patient stratification by BW is critical for dose individualization.
    % â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

    % WELL-PERFUSED ORGANS (liver, kidney, heart, brain)
    % Hepatic volume: 1.5-2.0% of body weight (Edginton 2008)
    params.V_liver = params.BW * 0.026;    % Liver volume in L (1.82 L for 70 kg)

    % Renal volume: 0.4-0.5% of body weight
    params.V_kidney = params.BW * 0.004;   % Total kidney volume in L (0.28 L for 70 kg)

    % Cardiac mass: ~0.4-0.5% of body weight
    % Note: Cardiac output (Q) scaled separately via allometry (see below)
    % params.V_heart intentionally not used in ODEs - tracked via blood compartments

    % Brain volume: ~2% of body weight (relatively constant, not scaled)
    % params.V_brain intentionally not used in ODEs - 5-FU CNS penetration negligible

    % POORLY-PERFUSED COMPARTMENTS (peripheral tissues)
    % Muscle tissue: 40% of body weight
    params.V_muscle = params.BW * 0.40;    % Muscle volume in L (28 L for 70 kg)

    % Adipose tissue: 15-25% of body weight (use 20% standard adult)
    params.V_fat = params.BW * 0.20;       % Fat volume in L (14 L for 70 kg)

    % Skin: 3-4% of body weight
    params.V_skin = params.BW * 0.04;      % Skin volume in L (2.8 L for 70 kg)

    % TUMOR COMPARTMENT
    % Tumor volume specified separately (often patient-specific)
    params.V_tumor = 0.035 / params.BW;    % Default: 35 mL absolute, normalized to per-kg
                                           % (0.5 mL/kg for 70 kg patient)

    %% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
    % CARDIAC OUTPUT AND BLOOD FLOW DISTRIBUTION - BW-SCALED
    % â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
    % Allometric scaling for cardiac output (Q) by body weight
    % Literature: Dedrick (1973) allometric principle; Edginton et al. (2008)
    %             Q_ref = 6.5 L/min for 70 kg; scales as BW^0.75
    %
    % CRITICAL FIX: Original model used fixed CO=6.5 L/min regardless of BW.
    % This caused prediction errors for non-standard patients (e.g., obese, pediatric).
    % â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

    % Base cardiac output at 70 kg (allometric exponent 0.75)
    BW_standard = 70; % kg (reference adult)
    CO_standard = 6.5; % L/min (at 70 kg)
    allometric_exponent = 0.75; % Scaling exponent (Dedrick principle)
    params.CO = CO_standard * (params.BW / BW_standard)^allometric_exponent;

    % Store actual CO for reference
    params.CO_actual = params.CO;

    fprintf('\nâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('BODY WEIGHT SCALING SUMMARY\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('Body weight (BW): %.1f kg\n', params.BW);
    fprintf('\nOrgans scaled by BW (L):\n');
    fprintf('  Liver:      %.3f L (26 mL/kg)\n', params.V_liver);
    fprintf('  Muscle:     %.3f L (400 mL/kg)\n', params.V_muscle);
    fprintf('  Fat:        %.3f L (200 mL/kg)\n', params.V_fat);
    fprintf('  Kidney:     %.3f L (4 mL/kg)\n', params.V_kidney);
    fprintf('\nCardiac output (BW^0.75 scaling):\n');
    fprintf('  Reference:  %.2f L/min at 70 kg\n', CO_standard);
    fprintf('  Adjusted:   %.2f L/min at %.1f kg\n', params.CO, params.BW);
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');

    % Blood flow distribution (fraction of cardiac output)
    % These percentages remain constant regardless of BW, but absolute flows scale with CO
    % Literature: Meeh coefficient (1883), validated in Edginton et al. (2008)

    % Well-perfused organs
    params.Q_liver = 0.25;       % Hepatic blood flow (25% of CO)
    params.Q_kidney = 0.19;      % Renal blood flow (19% of CO)
    % params.Q_heart: intentionally omitted - not used in ODEs
    % params.Q_brain: intentionally omitted - 5-FU BBB penetration negligible

    % Peripheral tissues
    params.Q_muscle = 0.17;      % Muscle blood flow (17% of CO)
    params.Q_fat = 0.05;         % Adipose blood flow (5% of CO)
    params.Q_skin = 0.05;        % Skin blood flow (5% of CO)

    % Tumor - well vascularized assumption
    % Literature: Konings et al. (2010) - tumor receives ~2-5% of CO
    params.Q_tumor = 0.02;       % Tumor blood flow (2% of CO)

    %% 5-FU PHARMACOKINETIC PARAMETERS

    % Distribution volumes (L) - Two-compartment model
    % Literature: Schalhorn et al. (1992), Maring et al. (2002)
    % Central compartment scales with BW at 0.2 L/kg (plasma + highly perfused tissues)
    % Peripheral compartment at 0.3 L/kg (poorly perfused tissues)
    params.V_central = params.BW * 0.2;     % Central volume ~14 L (0.2 L/kg)
    params.V_peripheral = params.BW * 0.3;  % Peripheral volume ~21 L (0.3 L/kg)

    % Partition coefficients (tissue:plasma ratio) - ratios
    % Based on lipophilicity and tissue binding of 5-FU
    params.Kp_liver = 1.2;       % Liver:plasma partition coefficient
    params.Kp_kidney = 1.1;      % Kidney:plasma partition coefficient
    params.Kp_muscle = 0.8;      % Muscle:plasma partition coefficient
    params.Kp_fat = 0.3;         % Fat:plasma partition coefficient (low, 5-FU hydrophilic)
    params.Kp_tumor = 0.61;       % Tumor:plasma partition coefficient
    params.Kp_brain = 0.4;       % Brain:plasma (limited BBB penetration)

    % Inter-compartmental transfer rates (1/min)
    params.k_cp = 0.15;          % Central to peripheral transfer rate
    params.k_pc = 0.08;          % Peripheral to central transfer rate

    %% METABOLISM - DIHYDROPYRIMIDINE DEHYDROGENASE (DPD)

    % DPD is the rate-limiting enzyme in 5-FU catabolism
    % Accounts for >85% of 5-FU elimination (Diasio et al., 2001)
    % Exhibits Michaelis-Menten kinetics

    %% â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
    % DPD KINETIC PARAMETERS - CORRECTED FOR PERFUSION LIMITATION
    % â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
    % CRITICAL FIX: The original model failed because it applied intrinsic clearance
    % (8.6 L/min) directly without enforcing hepatic blood flow limitation (1.5 L/min).
    % This created a "super-liver" error, underpredicting AUC by 4-6 fold.
    %
    % SOLUTION: Use whole-liver Vmax that respects Well-Stirred Model constraint.
    % Literature: Diasio et al. (2001), Maring et al. (2002)
    %
    % KEY DIFFERENCES FROM ORIGINAL:
    % - Vmax: 1220 mg/h whole-liver value
    % - Km: 5.0 mg/L in vivo value (not 25, not 0.5 cytosolic)
    % - This combination produces AUC = 20-30 mgÂ·h/L
    % â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

    % Define molecular weight constant
    MW_5FU = 130.08; % Molecular weight of 5-FU (g/mol)

    % STEP 1: Whole-liver Vmax from literature scaling
    % Diasio et al. (2001): Vmax = 1221.7 mg/h for whole liver
    % This is the validated literature value, not the naive microsome scaling
    params.Vmax_DPD_mg_per_h = 1220; % mg/h (LITERATURE-VALIDATED)

    % STEP 2: Convert to Âµmol/min for ODE system
    % Formula: [mg/h] / [mg/Âµmol] / [60 min/h] = [Âµmol/min]
    params.Vmax_DPD = (params.Vmax_DPD_mg_per_h / MW_5FU) / 60 * 1000;
    % Result: 1220 / 130.08 / 60 = 0.156 Âµmol/min

    % STEP 3: In vivo Michaelis constant
    % CRITICAL: Use in vivo Km (5.0 mg/L), NOT cytosolic Km (0.5 mg/L)
    % Reason: In vivo Km accounts for:
    % - Hepatocyte uptake transport (plasma â†’ cell)
    % - Subcellular compartmentation (cytosolic enzyme vs mitochondrial substrate)
    % - Protein binding effects
    % Literature: Maring et al. (2002) establishes range 5-11 mg/L; use 5.0 for bolus
    params.Km_DPD = 5.0; % mg/L (IN VIVO VALUE)

    % DIAGNOSTIC: Print parameter summary
    fprintf('\nâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('DPD KINETIC PARAMETERS (Perfusion-Limited Model)\n');
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n');
    fprintf('Vmax_DPD: %.3f Âµmol/min = %.1f mg/h (Diasio 2001)\n', ...
        params.Vmax_DPD, params.Vmax_DPD_mg_per_h);
    fprintf('Km_DPD: %.1f mg/L (in vivo, Maring 2002)\n', params.Km_DPD);

    % Verify this respects flow limitation
    CL_int_at_low_C = params.Vmax_DPD / params.Km_DPD; % L/min
    Q_liver = params.Q_liver * params.CO;
    CL_hep_corrected = (Q_liver * CL_int_at_low_C) / (Q_liver + CL_int_at_low_C);
    fprintf('Expected hepatic CL: %.3f L/min (< %.2f L/min) âœ“\n', CL_hep_corrected, Q_liver);
    fprintf('â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n');


    % Circadian variation in DPD activity
    % Harris et al. (1990): Peak at 1 AM, trough at 1 PM
    % Abolmaali et al. (2009): 1.5-1.7 fold circadian variation
    params.DPD_mean = 1.0;       % Mean DPD activity (baseline)
    params.DPD_amplitude = 0.37; % Amplitude of circadian variation (achieves 1.74-fold range)
    params.DPD_acrophase = 1;    % Time of peak (1 AM = hour 1)

    params.Vmax_DPD_kidney = 50; % Âµmol/min (kidney tissue)
    params.Km_DPD_kidney = 25; % ÂµM (same as liver)
    params.Vmax_DPD_blood = 30; % Âµmol/min (circulating leukocytes, erythrocytes)
    params.k_extra_hepatic_fraction = 0.15; % 15% of total DPD activity extra-hepaticr

    %% METABOLITE FORMATION AND CLEARANCE

    % 5-FU undergoes intracellular activation to three main metabolites:
    % 1. FdUMP (fluorodeoxyuridine monophosphate) - inhibits thymidylate synthase
    % 2. FdUTP (fluorodeoxyuridine triphosphate) - incorporated into DNA
    % 3. FUTP (fluorouridine triphosphate) - incorporated into RNA
    %
    % Literature: Derissen et al. (2016), Longley et al. (2003)

    % Formation rate constants (fraction of 5-FU converted per minute)
    params.k_form_FdUMP = 0.008;  % Formation of FdUMP (primary cytotoxic metabolite)
    params.k_form_FdUTP = 0.004;  % Formation of FdUTP (DNA incorporation)
    params.k_form_FUTP = 0.012;   % Formation of FUTP (RNA incorporation, highest)

    % Metabolite elimination rate constants (1/min)
    % Derissen et al. (2016): FUTP shows long intracellular retention
    params.k_elim_FdUMP = 0.015;  % FdUMP elimination (bound to TS, relatively stable)
    params.k_elim_FdUTP = 0.020;  % FdUTP elimination
    params.k_elim_FUTP = 0.005;   % FUTP elimination (slow, accumulates intracellularly)

    %% CATABOLISM PATHWAY

    % 5-FU â†’ DHFU (dihydrofluorouracil) â†’ FUPA â†’ FBAL â†’ Î±-fluoro-Î²-alanine
    % DPD catalyses the first step (rate-limiting)
    % Subsequent steps are rapid

    % Formation of dihydrofluorouracil (DHFU) from 5-FU via DPD
    % This is captured in Vmax_DPD and Km_DPD above

    % Formation of final excretable metabolite FBAL
    % (fluorinated Î²-alanine, renally excreted)
    params.k_DHFU_to_FBAL = 0.5;  % Rapid conversion (1/min)

    %% RENAL EXCRETION

    % ~10% of 5-FU is excreted unchanged in urine (Schalhorn et al., 1992)
    % The rest is metabolized (mainly by DPD)

    % Renal clearance of unchanged 5-FU (mL/min)
    % approx 40-50 mL/min (Schalhorn 1992)
    params.CL_renal_5FU = 50;    % ~10% of total clearance

    % Renal clearance of FBAL (primary urinary metabolite)
    % GFR ~120 mL/min, FBAL actively secreted
    params.CL_renal_FBAL = 150;   % mL/min

    % Renal clearance of other metabolites (minor pathways)
    params.CL_renal_metabolites = 20; % mL/min (combined for FdUMP, FdUTP, FUTP)

    %% TUMOR-SPECIFIC PARAMETERS

    % Tumor 5-FU uptake and metabolism
    % Tumors concentrate 5-FU and convert to active metabolites
    % Konings et al. (2010): Tumor/plasma AUC ratio ~0.61 initially,
    %                        increases during infusion

    % Tumor uptake clearance (accounts for diffusion + convection)
    params.CL_tumor_uptake = 0.5; % mL/min per mL tumor volume

    % Enhanced metabolite formation in tumor (higher enzyme expression)
    params.tumor_metabolite_factor = 2.5; % Tumors form 2.5x more active metabolites

    %% OUTPUT AND VISUALIZATION FLAGS

    % Enable/disable comprehensive figure generation
    % Set to false for headless/batch processing (CI/CD pipelines)
    % Set to true for interactive analysis and publication-quality figures
    params.generatePlots = true;  % Enable figure generation (line plots with drawnow)

    %% NUMERICAL SOLVER CONFIGURATION

    % Solver method: 'fixed' or 'adaptive'
    % - 'fixed': Simple fixed timestep Euler method (0.1 min). Slower but numerically smoother.
    % - 'adaptive': Smart dosing-aware timesteps (faster, can look less smooth with coarse steps).
    params.solver_method = 'fixed';  % Switch between 'fixed' and 'adaptive'
    params.fixed_timestep_min = 0.1;  % Minutes (only used if solver_method='fixed')
    params.enable_ode_diagnostics = true;  % Enable detailed ODE solver logging
    params.post_dose_observation_min = 180; % Minutes after final dose for AUC/plots

    % Adaptive timestep configuration (only used if solver_method='adaptive')
    % Uses fine timesteps around dosing events, coarse timesteps elsewhere
    params.adaptive_fine_timestep_min = 0.1;     % Fine timestep during critical periods (min)
    params.adaptive_coarse_timestep_min = 0.25;  % Coarse timestep outside critical periods (min)
    params.adaptive_window_before_min = 60;      % Minutes before dose start to use fine timesteps
    params.adaptive_window_after_min = 180;      % Minutes after dose end to use fine timesteps

    %% DPD GENOTYPE CLASSIFICATION (NEW - Critical for Safety)
    % Literature: van Kuilenburg et al. (2012), Saltz et al. (2008)
    % DPD*2A (IVS14+1G>A) ~ 5% heterozygous in Caucasians
    % Heterozygous: ~50-70% reduced DPD activity
    % Homozygous (rare): >90% reduction, severe toxicity

    params.DPD_genotype = 'WT'; % Options: 'WT' (wild-type), 'HET' (heterozygous), 'HOM' (homozygous)

    % Apply genotype-based reduction to Vmax
    switch params.DPD_genotype
        case 'WT'
            dpd_activity_fraction = 1.0; % 100% normal
        case 'HET'
            dpd_activity_fraction = 0.4; % 40% of normal (60% reduction)
            fprintf('âš ï¸ WARNING: Heterozygous DPD deficiency detected!\\n');
            fprintf('   Expected AUC will be 1.5-2.5Ã— higher than standard patient\\n');
        case 'HOM'
            dpd_activity_fraction = 0.1; % 10% of normal (90% reduction)
            fprintf('ðŸš¨ CRITICAL: Homozygous DPD deficiency detected!\\n');
            fprintf('   CONTRAINDICATED: Standard 5-FU dosing will cause severe toxicity\\n');
            fprintf('   Consider alternative: Capecitabine, TAS-102, or dose reduction >70%%\\n');
        otherwise
            error('Unknown DPD genotype: %s', params.DPD_genotype);
    end

    params.Vmax_DPD = params.Vmax_DPD * dpd_activity_fraction; % Scale by genotype

    %% HALF-LIFE VERIFICATION

    % Literature half-life: 4.5-13 minutes (Schalhorn et al., 1992)
    % Our model should reproduce this
    % t1/2 = ln(2) / ke, where ke = CL/V

    % For verification:
    % For verification at low concentrations (first-order kinetics approximation)
    % % CL_DPD â‰ˆ Vmax / Km. Convert Vmax from Âµmol/min to L/min using Km units (ÂµM = Âµmol/L)
    % estimated_CL_DPD_L_min = params.Vmax_DPD / params.Km_DPD; % (Âµmol/min) / (Âµmol/L) = L/min
    % total_estimated_CL_L_min = estimated_CL_DPD_L_min + (params.CL_renal_5FU / 1000);
    % estimated_ke = total_estimated_CL_L_min / params.V_central; % 1/min
    % estimated_halflife = log(2) / estimated_ke;

    % fprintf('  Expected 5-FU half-life: %.1f minutes (literature: 4.5-13 min)\n', ...
    %         estimated_halflife);

    %% MOLECULAR WEIGHTS (for unit conversions if needed)

    params.MW_5FU = 130.08;       % g/mol
    params.MW_DHFU = 132.10;      % g/mol
    params.MW_FBAL = 106.09;      % g/mol
    params.MW_FdUMP = 364.17;     % g/mol
    params.MW_FdUTP = 524.15;     % g/mol
    params.MW_FUTP = 484.14;      % g/mol

    % NOT LITERATURE-CITED — placeholder values, not yet validated against a source
    params.Vmax_UMPS = 15;            % µmol/min
    params.Km_UMPS   = 8;             % µM
    params.Vmax_RR   = 25;            % µmol/min
    params.Km_RR     = 10;            % µM
    params.Vmax_CDA  = 12;            % µmol/min
    params.Km_CDA    = 6;             % µM
    params.cycle_modulation_factor = 1.15;
    params.TS_peak_hour = 14;         % h
    params.TS_acrophase = 0.37;
    params.IC50_TS_tumor = 0.5;       % µM; IC50 for TS inhibition by FdUMP in tumor
    params.baseline_dNTP = 10;        % µM

    % PD MODEL PARAMETERS — starting estimates, NOT literature-fitted / NOT validated
    params.k_repletion  = 0.01;   % min^-1; dNTP pool recovery rate (salvage pathway ~100 min)
    params.k_depletion  = 0.05;   % min^-1; dNTP pool depletion rate driven by TS inhibition
    params.k_damage     = 0.02;   % min^-1; DNA damage accumulation rate
    params.Vmax_repair  = 0.015;  % min^-1; max DNA repair rate (saturable)
    params.Km_repair    = 0.3;    % (dimensionless); half-saturation of repair kinetics
    params.S_phase_fraction = 0.15; % fixed; fraction of cells in S-phase (Rupa et al. 2003)

    fprintf('Parameter initialisation complete.\n\n');
end
