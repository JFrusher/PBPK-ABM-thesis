# PBPK-ABM Dissertation Compendium

Research Repository

![MATLAB](https://img.shields.io/badge/MATLAB-R2024a%2B-0076A8?logo=mathworks)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-[INSERT_DESCRIPTION]-276DC3?logo=r)
![PhysiCell](https://img.shields.io/badge/Engine-PhysiCell%2FBioFVM-005f73)
![License](https://img.shields.io/badge/License-BSD--3--Clause-94d2bd)
![Reproducibility](https://img.shields.io/badge/Reproducibility-Documented-success)

---

## Abstract

This repository is a research compendium for a dissertation on integrated physiologically based pharmacokinetic and agent-based modeling of 5-fluorouracil dosing. The workflow links simulation, parameter exploration, and post-processing pipelines across C++, MATLAB, and Python to evaluate exposure, target attainment, and regimen behavior. The structure and scripts are organized to support transparent, auditable, and reproducible computational methods.

---

## Research Overview

The central research objective is to identify and evaluate physiological and regimen parameters that govern 5-FU exposure profiles and downstream treatment-relevant metrics across central and tumor compartments. In practical terms, this compendium supports:

1. Mechanistic PBPK simulation with circadian modulation and metabolite tracking.
2. Multicellular tumor-environment simulation through PhysiCell/BioFVM.
3. Monte Carlo and regimen-comparison analyses for target attainment and sensitivity interpretation.

---

## Directory Map (Navigation First)

| Path | Role in Workflow |
| --- | --- |
| main.cpp | Primary PhysiCell simulation entry point. |
| core/ | PhysiCell core dynamics and model engine source files. |
| modules/ | PhysiCell standard output, geometry, and settings modules. |
| custom_modules/ | Project-specific custom model behavior (for example custom.cpp/custom.h). |
| BioFVM/ | Diffusion and transport framework dependencies used by PhysiCell. |
| config/ | Simulation configuration and cell-rule inputs (XML/CSV). |
| run5FU_PBPK_Simulation.m | Main PBPK solver and PK output generation routine. |
| MC_5FU_PK_sensitivity.m | Monte Carlo sensitivity simulation logic. |
| analyze_MC10k_comprehensive.m | Comprehensive analysis of merged 10k Monte Carlo outputs. |
| compare_MC10k_regimens.m | Comparative analysis between regimen cohorts. |
| merge_MC_chunk_results.py | Merges chunked Monte Carlo outputs into unified datasets. |
| export_outputs_to_excel.py | Converts PhysiCell XML outputs to analysis-ready CSV tables. |
| tests/ and unit_tests/ | Test assets and validation scaffolding. |
| run_MC_10k_array.sh and submit_MC_10k.sh | SLURM/HPC batch orchestration scripts. |
| MC_10K_HPC_Operations_Guide.md | Operational runbook for 10k Monte Carlo execution on HPC. |

---

## Repository Structure

```text
PBPK-ABM-thesis/
├─ main.cpp
├─ Makefile
├─ core/
├─ modules/
├─ custom_modules/
├─ BioFVM/
├─ config/
│  ├─ PhysiCell_settings.xml
│  ├─ cell_rules.csv
│  └─ colorectal_tumor_cells.csv
├─ matlab/
├─ tests/
│  ├─ cases/
│  ├─ system/
│  ├─ timing/
│  └─ unit/
├─ unit_tests/
├─ protocols/
├─ WP2/
├─ run5FU_PBPK_Simulation.m
├─ MC_5FU_PK_sensitivity.m
├─ analyze_MC10k_comprehensive.m
├─ compare_MC10k_regimens.m
├─ merge_MC_chunk_results.py
├─ export_outputs_to_excel.py
├─ RUN5FU_PBPK_DISSERTATION_GUIDE.md
├─ MC_10K_HPC_Operations_Guide.md
├─ CITATION.cff
└─ LICENSE
```

---

## Tech Stack

Primary languages and computational tooling used in-repo:

- MATLAB: PBPK simulation, parameter analysis, statistical summaries, and publication figures.
- Python: output extraction, dataset merging, tabular transformations, and utility scripts.
- C++ (PhysiCell/BioFVM): multicellular simulation and microenvironment dynamics.
- R: [INSERT DESCRIPTION: if/where R analysis is executed outside the tracked repository files].

Primary Python libraries visible in scripts include numpy and pandas, with optional pcdl integration for pyMCDS parsing.

---

## Methodology and Key Scripts

| Research Question | Script(s) Used | Output Artifacts |
| --- | --- | --- |
| How does a given dosing regimen affect 5-FU and metabolite concentration-time profiles? | run5FU_PBPK_Simulation.m, readDosingRegimen.m, calculateDosingRate.m | Compartment and metabolite CSV outputs, PK summary metrics, figures. |
| Which parameter uncertainties most influence AUC and Cmax at scale? | MC_5FU_PK_sensitivity.m, run_MC_10k_array.sh, submit_MC_10k.sh, merge_MC_chunk_results.py | Chunk outputs, merged Monte Carlo tables, sensitivity-ready datasets. |
| What global patterns emerge from 10k Monte Carlo runs? | analyze_MC10k_comprehensive.m, check_mc_quality.m | QC reports, distribution plots, target-attainment tables, tornado and correlation figures. |
| How do alternative regimens compare for exposure and target attainment? | compare_MC10k_regimens.m | Comparative summaries, percentile shifts, regimen-delta sensitivity plots. |
| How are ABM raw outputs transformed into analyzable time-series tables? | export_outputs_to_excel.py, add_regimen_time_column.py | Time-aligned CSV outputs, curve-fitter tables, metadata logs. |

---

## Reproducibility Guide

### 1. Core Environment Requirements

- OS: Linux or Windows for analysis; Linux cluster recommended for large batch campaigns.
- Compiler: GNU g++ with OpenMP support.
- MATLAB: version compatible with scripts in this repository (cluster guide references matlab/2025a).
- Python: 3.10+ with pandas and numpy.

### 2. Build the Simulation Engine

```bash
make
```

This builds the executable named project from main.cpp and linked PhysiCell/BioFVM sources.

### 3. Python Environment Setup (Conda Example)

```bash
conda create -n pbpk-abm python=3.10 -y
conda activate pbpk-abm
pip install numpy pandas
pip install pcdl
```

If pcdl is unavailable in your environment, export_outputs_to_excel.py can fall back to beta/pyMCDS.py.

### 4. Run Typical Workflow

1. Execute simulation batches with configuration from config/ (for example PhysiCell_settings.xml).
1. Export raw output XML files to tabular CSV:

```bash
python export_outputs_to_excel.py --folder simulation_results_2 --outdir Physicell_results_2
```

1. Run PBPK and Monte Carlo analyses in MATLAB:

```matlab
results = run5FU_PBPK_Simulation('dosing_example.csv', 'runA');
analyze_MC10k_comprehensive('MC_10k_merged_outputs.csv');
```

1. For HPC operation and restart logic, follow MC_10K_HPC_Operations_Guide.md.

### 5. Reproducibility Checklist for Reporting

1. Record commit hash and branch for every analysis pack.
2. Archive exact command invocations and parameter overrides.
3. Store generated CSV, figures, and logs together per run.
4. Document MATLAB version, Python version, and compiler details.
5. Confirm quality checks (mass-balance checks, NaN checks, target attainment summaries).

---

## Key Results (Figure Placeholders)

Use this section to surface dissertation headline outputs.

### Figure 1. Target Attainment

[INSERT FIGURE: Target Attainment Summary]

### Figure 2. Tornado Plot Sensitivity

[INSERT FIGURE: Tornado Plot of Parameter Influence]

### Figure 3. Exposure Space

[INSERT FIGURE: Exposure Space (AUC vs Cmax or equivalent)]

---

## How to Cite

Repository-level citation metadata is available in CITATION.cff and CITATION.txt.

Dissertation BibTeX placeholder:

```bibtex
@phdthesis{[INSERT_CITATION_KEY],
  author  = {[INSERT AUTHOR NAME]},
  title   = {[INSERT DISSERTATION TITLE]},
  school  = {[INSERT UNIVERSITY]},
  year    = {[INSERT YEAR]},
  address = {[INSERT LOCATION]},
  note    = {Computational workflow repository: PBPK-ABM-thesis}
}
```

---

## License

This repository is distributed under the BSD 3-Clause License. See LICENSE for full terms.
