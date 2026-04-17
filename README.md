# PBPK-ABM Dissertation Compendium

Integrated PBPK-ABM workflow for 5-FU simulation, sensitivity analysis, and dissertation reporting.

![MATLAB](https://img.shields.io/badge/MATLAB-R2024a%2B-0076A8?logo=mathworks)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![PhysiCell](https://img.shields.io/badge/Engine-PhysiCell%2FBioFVM-005f73)
![License](https://img.shields.io/badge/License-BSD--3--Clause-94d2bd)

## Scope

This repository supports three linked workflows:

1. PhysiCell/BioFVM multicellular simulations.
2. MATLAB PBPK simulation and Monte Carlo sensitivity studies.
3. Python-based data extraction and table generation.

## Repository Layout

```text
PBPK-ABM-thesis/
|- main.cpp
|- Makefile
|- config/
|- core/
|- modules/
|- custom_modules/
|- BioFVM/
|- matlab/
|  |- pbpk/
|  |- analysis/
|  |- visualization/
|  |- diagnostics/
|  |- tests/
|- scripts/
|  |- hpc/
|  |- python/
|- docs/
|  |- guides/
|  |- project/
|  |- references/
|- archive/
|- CITATION.cff
|- LICENSE
`- README.md
```

## Key Entrypoints

| Path | Purpose |
| --- | --- |
| main.cpp | Primary PhysiCell model entry point. |
| matlab/pbpk/run5FU_PBPK_Simulation.m | Main PBPK simulation function. |
| matlab/analysis/MC_5FU_PK_sensitivity.m | Monte Carlo sensitivity workflow. |
| matlab/analysis/analyze_MC10k_comprehensive.m | Post-merge 10k analysis. |
| scripts/python/export_outputs_to_excel.py | Convert PhysiCell output to analysis CSVs. |
| scripts/python/merge_MC_chunk_results.py | Merge chunked MC outputs. |
| scripts/hpc/run_MC_10k_array.sh | Chunk worker for SLURM arrays. |
| scripts/hpc/submit_MC_10k.sh | Array submit helper and manifest writer. |

## Quick Start

### 1. Build engine

```bash
make
```

### 2. Python setup

```bash
conda create -n pbpk-abm python=3.10 -y
conda activate pbpk-abm
pip install numpy pandas pcdl
```

### 3. Export PhysiCell results

```bash
python scripts/python/export_outputs_to_excel.py --folder simulation_results_2 --outdir Physicell_results_2
```

### 4. Run MATLAB PBPK analysis

```matlab
addpath(genpath('matlab'));
results = run5FU_PBPK_Simulation('dosing_example.csv', 'runA');
summary = analyze_MC10k_comprehensive('MC_10k_merged_outputs.csv');
```

### 5. Submit chunked MC run on HPC

```bash
bash scripts/hpc/submit_MC_10k.sh
```

## Reproducibility Checklist

1. Record commit hash and branch for every analysis pack.
2. Save exact run commands and environment versions.
3. Keep outputs grouped per run (CSV, figures, logs).
4. Document MATLAB, Python, and compiler versions.
5. Run quality checks before reporting (mass balance, NaN checks, target attainment summaries).

## Key Results (Figure Placeholders)

### Figure 1. Target Attainment

[INSERT FIGURE: Target Attainment Summary]

### Figure 2. Tornado Plot Sensitivity

[INSERT FIGURE: Tornado Plot of Parameter Influence]

### Figure 3. Exposure Space

[INSERT FIGURE: Exposure Space (AUC vs Cmax or equivalent)]

## Citation

Repository citation metadata:

1. CITATION.cff
2. docs/references/CITATION.txt

## License

BSD 3-Clause. See LICENSE.
