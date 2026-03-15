# PBPK-ABM Thesis Model

[![PhysiCell](https://img.shields.io/badge/engine-PhysiCell%201.14.2-005f73)](https://github.com/MathCancer/PhysiCell)
[![Language](https://img.shields.io/badge/core-C%2B%2B%20%7C%20MATLAB%20%7C%20Python-0a9396)](#technology-stack)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-94d2bd)](LICENSE)

Integrated pharmacokinetic and agent-based simulation for **5-FU treatment dynamics** in the tumor microenvironment, combining:

- PhysiCell/BioFVM simulation core (C++)
- PBPK and analysis workflows (MATLAB)
- high-throughput post-processing and export tooling (Python)

This repository is structured to support both scientific review (reproducibility, citation, model traceability) and engineering review (software quality, automation, maintainability).

## Why This Repository Matters

This project demonstrates an end-to-end computational oncology workflow:

- translating biological hypotheses into executable simulation models
- integrating mechanistic PBPK and multicellular ABM pipelines
- scaling experiments through batch and cluster workflows
- building analysis outputs suitable for publication and dissertation defense

For academic reviewers, this repo emphasizes methodological transparency.
For industry/hiring reviewers, it highlights practical scientific software engineering across C++, MATLAB, Python, reproducibility, and data workflow design.

## Highlights

- **Hybrid PBPK-ABM modeling** for 5-FU transport, metabolism, and tumor exposure.
- **Circadian-aware PK support** via dosing and metabolism timing workflows.
- **Batch and Monte Carlo execution** using cluster-friendly scripts.
- **Multi-format analysis outputs** from XML snapshots to aligned CSV tables.
- **Publication-oriented figure tooling** and reproducible analysis helpers.

## Repository Map

| Path | Purpose |
|---|---|
| `main.cpp` | PhysiCell simulation entrypoint and runtime loop |
| `custom_modules/` | Model-specific ABM logic and substrate effects |
| `config/` | XML/csv settings, rules, and simulation configuration |
| `run5FU_PBPK_Simulation.m` | Core PBPK simulation and PK metric generation |
| `export_outputs_to_excel.py` | Extracts `output*.xml` into analysis-ready CSVs |
| `tests/`, `unit_tests/` | Validation assets and reference test material |
| `run_sim.sh`, `run_array.sh`, `run_MC.sh` | HPC/SLURM execution scripts |
| `RUN5FU_PBPK_DISSERTATION_GUIDE.md` | Methods-focused guide for dissertation reporting |

## Technology Stack

- C++11 with OpenMP (`make`-based build)
- MATLAB (PBPK simulation, diagnostics, plotting)
- Python 3.10+ (`numpy`, `pandas`, optional `pcdl`)
- PhysiCell/BioFVM framework foundation

## Quick Start

### 1. Build the simulation binary

```bash
make
```

### 2. Run with default configuration

```bash
./project ./config/PhysiCell_settings.xml
```

### 3. Export simulation outputs to CSV

From extracted output folders:

```bash
python export_outputs_to_excel.py --folder simulation_results_2 --outdir Physicell_results_2
```

From zip archives (persistent extraction enabled):

```bash
python export_outputs_to_excel.py --zip simulation_results_2.zip --extract-zip-to-folder --outdir Physicell_results_2
```

Performance-oriented export example:

```bash
python export_outputs_to_excel.py --folder simulation_results_2 --every-nth 5 --max-files 20 --workers 4 --fast-mode --outdir Physicell_results_fast
```

## Reproducible Workflow

1. Build and run PhysiCell experiments from versioned config inputs.
2. Export XML outputs to normalized CSV artifacts.
3. Run PBPK and aggregate analyses in MATLAB.
4. Generate publication/dissertation figures and metrics.
5. Archive run metadata: command, parameters, environment, commit hash.

## Running on HPC

SLURM scripts are included for parallel and high-throughput execution:

- `run_sim.sh`: single optimized simulation run
- `run_array.sh`: array jobs for repeated scenario execution
- `run_MC.sh`: Monte Carlo sensitivity workflow

Before cluster submission, verify account/partition/module settings in each script.

## Quality and Testing

- GitHub Actions workflows are available in `.github/workflows/` for build/test automation.
- Additional MATLAB checks and benchmark-style tests are available in root scripts and test folders.
- Use smoke tests for large datasets with exporter options `--max-files`, `--every-nth`, and `--workers`.

## Citation

If you use this repository in academic work:

- cite this project using `CITATION.cff`
- cite PhysiCell and BioFVM as listed in `CITATION.txt`
- include dissertation-specific methodological references where relevant

## License and Third-Party Notice

This repository is distributed under the **BSD 3-Clause License**. See `LICENSE`.

It includes and extends components from PhysiCell/BioFVM ecosystems; corresponding notices are retained in:

- `licenses/PhysiCell.txt`
- `licenses/BioFVM.txt`
- `licenses/MaBoSS.txt`
- `licenses/pugixml.txt`

## Professional Review Checklist

For collaborators, examiners, or hiring reviewers, start with:

1. `README.md` for architecture and workflow orientation
2. `RUN5FU_PBPK_DISSERTATION_GUIDE.md` for methods detail
3. `export_outputs_to_excel.py` and `run5FU_PBPK_Simulation.m` for analysis pipeline depth
4. `.github/workflows/tests.yml` for automated verification scope

## Contact

For reproducibility questions, feature suggestions, or collaboration, open a GitHub issue with:

- scenario/config details
- command used
- expected vs observed behavior
- relevant logs or output snippets
