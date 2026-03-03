# PBPK-ABM Thesis Model (5-FU / Tumor Microenvironment)

A research-focused PhysiCell-based modeling repository for integrated PBPK + ABM simulation, parameter exploration, and post-processing workflows (MATLAB + Python).

This repository is tailored to your thesis-specific model structure, configuration files, scripts, and analysis pipeline.

---

## Overview

This project extends a PhysiCell codebase into a custom, end-to-end modeling workflow that includes:

- custom C++ simulation logic in `main.cpp` and `custom_modules/`
- scenario configuration under `config/`
- simulation orchestration scripts (`run_sim.sh`, `run_array.sh`, `run_MC.sh`)
- MATLAB analysis and plotting scripts (`run5FU_PBPK_Simulation.m`, `analyze_batch_simulations.m`, etc.)
- Python export utilities for converting PhysiCell outputs to analysis-ready tables (`export_outputs_to_excel.py`)

The goal is fast, reproducible simulation studies for 5-FU dynamics and tumor behavior analysis.

---

## Key Capabilities

- **Integrated PBPK-ABM workflow** for treatment and tumor response studies.
- **Batch execution support** for parameter sweeps and repeated runs.
- **Monte Carlo / sensitivity scripts** for uncertainty exploration.
- **Automated CSV export pipeline** from `output*.xml` (via `pyMCDS` parsing).
- **Publication-oriented plotting utilities** and figure styling helpers.

---

## Repository Structure (high level)

- `main.cpp`, `main-backup.cpp`  
  Main simulation entry points.
- `custom_modules/`  
  Project-specific cell rules and simulation logic.
- `config/`  
  PhysiCell settings, cell definitions, and rules tables.
- `core/`, `BioFVM/`, `modules/`, `addons/`  
  Engine / framework internals and optional integrations.
- `output/`  
  Default simulation outputs.
- `export_outputs_to_excel.py`  
  Converts PhysiCell outputs to time-aligned CSV summaries.
- MATLAB scripts in repo root (`*.m`)  
  PBPK simulation, diagnostics, post-processing, and plotting.

---

## Requirements

### Core simulation

- C++ compiler with OpenMP support
- `make`
- PhysiCell-compatible build toolchain

### Python export pipeline

- Python 3.10+ (tested in your environment with newer Python as well)
- `numpy`, `pandas`
- `pcdl` (optional; script auto-falls back to local `beta/pyMCDS.py` if unavailable)

### MATLAB workflow

- MATLAB for PBPK scripts, analyses, and plotting

---

## Quick Start

### 1) Build simulation

```bash
make
```

### 2) Run simulation

```bash
./project
```

(Executable name depends on your active Makefile target/config.)

### 3) Export PhysiCell outputs to CSV

From extracted output folder:

```bash
python export_outputs_to_excel.py --folder simulation_results_2 --outdir Physicell_results_2
```

From zip and persist extraction folder:

```bash
python export_outputs_to_excel.py --zip simulation_results_2.zip --extract-zip-to-folder --outdir Physicell_results_2
```

---

## Fast Export / Performance Options

`export_outputs_to_excel.py` supports speed-oriented flags for large studies:

- `--every-nth N` : process every Nth timestep
- `--max-files N` : cap number of processed XML files
- `--workers N` : parallel processing workers
- `--fast-mode` : skip expensive microenvironment/spatial/attribute stats

Example (fast smoke test):

```bash
python export_outputs_to_excel.py --folder simulation_results_2 --every-nth 5 --max-files 20 --workers 4 --fast-mode --outdir Physicell_results_fast
```

---

## Typical Analysis Workflow

1. Run simulation(s) to generate PhysiCell outputs.
2. Export time-aligned CSVs using `export_outputs_to_excel.py`.
3. Run MATLAB analysis scripts for summary statistics and visualizations.
4. Generate publication figures with helper scripts:
   - `set_publication_figure_style.m`
   - `save_publication_figure.m`
   - `MC_5FU_publication_helpers.m`

---

## Reproducibility Notes

- Keep scenario-specific inputs under `config/` version-controlled.
- Use explicit output directories per run batch.
- Archive run metadata (date, parameters, script command, commit hash) with results.
- Prefer scripted runs (`run_array.sh`, `run_MC.sh`) over manual execution for sweeps.

---

## Troubleshooting

### Python parser import issues (`pcdl`)

If `pcdl` fails to import dependencies, the exporter falls back automatically to local `beta/pyMCDS.py`.

### Large zip extraction interruptions

If zip extraction is slow/interrupted, extract once to folder and run with `--folder` mode to avoid repeated extraction.

### Performance constraints

Start conservative (`--workers 2`) and increase gradually based on RAM/CPU and dataset size.

---

## Attribution and Upstream Contributions

This repository is built on top of the **PhysiCell** ecosystem and includes adapted framework components.

Please acknowledge PhysiCell in derivative scientific work:

- A. Ghaffarizadeh, R. Heiland, S. H. Friedman, S. M. Mumenthaler, P. Macklin,  
  *PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems*,  
  PLoS Computational Biology 14(2): e1005991 (2018).  
  DOI: https://doi.org/10.1371/journal.pcbi.1005991

Upstream project:

- PhysiCell: https://github.com/MathCancer/PhysiCell

This README intentionally focuses on your thesis-specific model and workflow while preserving explicit credit to upstream contributors.

---

## Citation

If this repository is used in reports, thesis chapters, or publications, include:

- the thesis/project citation details
- the PhysiCell core citation above
- any domain-specific model references in `CITATION.txt`

---

## Contact / Project Notes

For internal project usage, consider maintaining:

- a changelog of experimental scenarios
- a run registry (parameter sets, seed strategy, commit hash)
- a short methods note linking scripts to figures/tables

This keeps thesis results auditable and easier to defend/reproduce.
