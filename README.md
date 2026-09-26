# Module-Selection Balance in the Evolution of Modular Organisms

**Authors:** Minkyu Kim, Sarah M. Ardell, Sergey Kryazhimskiy  
**Affiliations:** Cornell University; University of California San Diego

---

## Overview

This repository contains MATLAB code and a Python notebook for the simulations and analyses in:

> Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025). Module-Selection Balance in the Evolution of Modular Organisms.

Genotype-phenotype-fitness map (GPFM) models are implemented in the directories below — pleiotropic, modular, discordant-module and nested FGM in the main tree, with a constant-supply model, a restricted mutational cone and a ten-module version added as extensions. Each is simulated under the successive mutations (SSWM) and concurrent mutations (CM) regimes. A Python script reproduces the LTEE metagenomic analysis (Figure 7).

---

## Data Availability

Pre-computed simulation output files (`.mat`) will be deposited on Zenodo upon acceptance. To regenerate all figures directly from these files without rerunning simulations, download the archive, extract the contents into the project root so that `results/` and `results/Supplementary/` are present, then run:

```matlab
cd /path/to/project
addpath(genpath(pwd));          % the Run_ scripts live in run_scripts/
Run_mainFigures('reproduce')
Run_supplementaryFigures('reproduce')
```

To reproduce the simulations from scratch (~10 hours on 10 cores), see [Reproducing Paper Results](#reproducing-paper-results) below.

---

## Repository Structure

```
project_root/
├── README.md
├── LICENSE
├── reproduce_all.m        # The one entry point: all simulations and figures
│
├── run_scripts/           # One Run_ script per model or figure set
├── simulation_scripts/    # Wright-Fisher and Gillespie simulators
├── analysis_scripts/      # Trajectory averaging, analytical predictions, diagnostics
├── figure_scripts/        # Figure producers
├── utils/                 # Parameters, genome initialization, helpers
├── docs/                  # One note per model, plus the figure manifest
├── LTEE_analysis/         # Figures 7 and S6 (Python, with its own data and results)
│
├── results/               # All MATLAB output
│   ├── SSWM/  CM_Asexual/  CM_Sexual/     # pleiotropic and modular
│   ├── ConstantSupply/                     # constant-supply control
│   ├── RestrictedTheta_sampledinit/        # restricted cone, walk-and-swap init
│   ├── RestrictedTheta_maxentinit/         # restricted cone, maximum-entropy init
│   ├── Generalization/NestedFGM/
│   ├── nModule/                            # ten-module simulations
│   ├── Supplementary/                      # parameter-grid data and its figures
│   ├── Figures/
│   └── main_figures/                       # the manuscript set, collected
│
└── _relegated/            # Retained but inactive; see _relegated/MANIFEST.md
```

`reproduce_all.m` sets up its own path. Every other entry point is a function
inside one of these directories, so a bare `cd` to the project root is not
enough - run `addpath(genpath(pwd));` once per MATLAB session first.

Every model is laid out the same way: its runner in `run_scripts/`, its
simulators in `simulation_scripts/`, its predictions and diagnostics in
`analysis_scripts/`, its figures in `figure_scripts/`, and its output in a
subdirectory of `results/`. `docs/` holds one note per model explaining what it
is and why it exists; `docs/figures.md` maps manuscript figure numbers onto the
scripts that produce them.

---

## Requirements

### MATLAB (Figures 2-6, S1-S5, S7-S8)
- MATLAB R2025a or later

| Toolbox | Purpose |
|---------|---------|
| Statistics and Machine Learning | `mnrnd`, `poissrnd`, `datasample` |
| Symbolic Math | `syms`, `solve` in `findInitialPhenotypes` |
| Parallel Computing | `parfor` acceleration |

### Python (Figures 7 and S6)
- Python 3.8 or later
- Dependencies: `numpy`, `pandas`, `matplotlib`
- LTEE metagenomic data (see below)

---

## LTEE Data Setup (Figure 7)

The LTEE analysis requires data from Good et al. (2017), which is not included in this repository.

1. Download the repository ZIP from https://github.com/benjaminhgood/LTEE-metagenomic
2. Unzip and place the resulting `LTEE-metagenomic-master` folder inside `LTEE_analysis/`.

The expected directory layout is:

```
LTEE_analysis/
├── LTEE-metagenomic-master/   ← place downloaded data here
│   └── data_files/
├── ltee_analysis.py
└── results/
```

Then:

```
cd LTEE_analysis
python ltee_analysis.py --download    # first run, fetches the data files
python ltee_analysis.py               # subsequent runs
```

Pass `--data-dir` to point at a copy held elsewhere.

---

## Reproducing Paper Results

To verify the pipeline before a full run, at test scale, every simulation and
every figure script:

```matlab
cd /path/to/project
test_all
```

`test_all` mirrors `reproduce_all` stage for stage at reduced parameter sets and
writes everything to `results_test/`, so nothing in `results/` is touched. It
reports each stage's status, wall-clock time, and whether that stage actually
wrote a fresh output; one failure does not stop the rest. Minutes, not hours.

To reproduce all simulations and figures (~10 hours on 10 cores):

```matlab
cd /path/to/project
reproduce_all
```

Output files are written to `results/`, `results/Supplementary/`, and their respective `Figures/` subdirectories.

---

## Running Individual Models

Each driver script accepts `'reproduce'` (default) or `'test'` as the mode argument.

```matlab
cd /path/to/project
addpath(genpath(pwd));
Run_pleiotropicFGM('reproduce', 'sampled')   % pleiotropic arm of Figures 2-4
Run_modularFGM('reproduce')                  % modular arm of Figures 2-4
Run_nestedFGM('reproduce')                   % Figure 5, S4  - symmetric [n1=10, n2=10]
Run_nestedFGM('reproduce', {}, [10, 20])     % Figure S5     - asymmetric [n1=10, n2=20]
Run_modularConstantSupply                    % Figure 6
Run_pleiotropicRestrictedTheta('init','maxent')  % restricted cone
Run_supplementary('reproduce')               % Figures S1, S2
Run_ThresholdDAnalysis('reproduce')          % Figure S3
Run_nModule                                  % Figures S7, S8
```

After simulations complete, generate figures with:

```matlab
addpath(genpath(pwd));                 % if not already done
Run_allMainFigures                     % Figures 2-6 and the restricted cone
Run_supplementaryFigures('reproduce')  % Figures S1-S5
```

Individual figures can be regenerated without rerunning the others:

```matlab
Run_mainFigures('reproduce', 3)           % Figure 3 only
Run_supplementaryFigures('reproduce', 3)  % Figure S3 only
```

---

## Figure Index

`docs/figures.md` holds the full mapping, including output paths.

### Main figures

| Figure | Description | Produced by |
|--------|-------------|-------------|
| 1 | Schematic of the two genotype-phenotype-fitness maps | drawn by hand |
| 2 | Successive mutations, pleiotropic and modular | `Run_mainFigures` |
| 3 | Concurrent mutations, linked chromosomes | `Run_mainFigures` |
| 4 | Concurrent mutations, unlinked chromosomes | `Run_mainFigures` |
| 5 | Nested FGM | `Run_mainFigures` |
| 6 | Constant supply of module-improving mutations | `Run_modularConstantSupply` |
| 7 | LTEE metagenomic analysis | `ltee_analysis.py` |

`Run_allMainFigures` redraws figures 2 to 6 in one call and collects them.
The restricted-cone figure (`makeFigure_RestrictedTheta`) is produced but its
manuscript number is not yet settled.

### Supplementary figures

| Figure | Description | Produced by |
|--------|-------------|-------------|
| S1 | Numerical validation of the two-module adaptation rates | `Run_supplementaryFigures` |
| S2 | Accuracy of the heuristic approximation under different weights | `Run_supplementaryFigures` |
| S3 | Threshold D sensitivity | `Run_supplementaryFigures` |
| S4 | Nested FGM mutation-effect distributions | `Run_supplementaryFigures` |
| S5 | Asymmetric nested FGM dynamics | `Run_supplementaryFigures` |
| S6 | LTEE multi-hit gene sensitivity | `ltee_analysis.py` |
| S7 | Ten functional modules | `Run_nModule` |
| S8 | Evolutionary stalling with ten modules | `Run_nModule` |

`Run_supplementary` and `makeFigureS_SteadyStateCM` produce a steady-state
concurrent-mutations validation that the current supplement does not cite. The
code is retained; the figures are not part of the manuscript.

---

## Key Parameters

| Parameter | Symbol | Default | Description |
|-----------|--------|---------|-------------|
| `popSize` | $N$ | $10^4$ | Population size |
| `mutationRate` | $U$ | model-dependent | Genome-wide mutation rate (see below) |
| `deltaTrait` | $\delta$ | 0.1 | Mutational step size (or mutation-vector magnitude $m = \sqrt{2\delta}$ for the nested FGM) |
| `landscapeStdDev` | $\sigma$ | 2 | Fitness landscape width |
| `ellipseRatio` | $a_1/a_2$ | $\sqrt{2}$ | Selection anisotropy |
| `geneticTargetSize` | $[L_1, L_2]$ | $[200, 200]$ | Number of loci per module (pleiotropic, modular, and discordant models only) |
| `recombinationRate` | $\rho$ | 0 or 1 | Recombination rate |

Parameters are set in each `Run_*.m` script and passed to simulations via `initializeSimParams`.

### Mutation rate parameterization

Mutations are parameterized by the per-locus rate $\mu$, set to $5\times10^{-9}$ (SSWM) and $10^{-5}$ (CM) as described in the paper. Because the four GPFM models differ in whether they have explicit loci, the genome-wide rate $U$ passed to the simulation code differs across models:

- **Pleiotropic GPFM**: has $2L = 400$ explicit loci. The code passes $U = 2\mu L$, giving $U = 2\times10^{-6}$ (SSWM) and $4\times10^{-3}$ (CM). Implemented in `Run_pleiotropicFGM` as `mutationRateSlow = 1e-7 * (L/K)` with $K = 10$, $L = 200$.
- **Modular and discordant GPFMs**: these models track traits rather than individual loci, and the per-locus rate $\mu$ is recovered internally from $U$ and the `geneticTargetSize` parameter $L_i$. In the SSWM regime, these models use the `initializeSimParams` default $U = 10^{-7}$. In the CM regime, they use $U = 2\times10^{-4} \times (L/K) = 4\times10^{-3}$.
- **Nested FGM**: has no explicit loci and uses a continuous Gaussian mutation framework. Mutations arise at rate $U/2$ per module. The code passes $U$ directly: $U = 10^{-7}$ (SSWM) and $U = 2\times10^{-4}$ (CM).

Note that in the SSWM regime, the per-locus rate $\mu$ only affects the timescale of adaptation, not the shape of the evolutionary trajectories in trait space (see equations in the paper). In the CM regime, $\mu$ enters the Desai-Fisher function and does affect trajectory shape.

### Other implementation notes

- **Beneficial mutation threshold**: all SSWM simulations use $s > 0$ as the threshold for a mutation to be considered beneficial, consistently across all four GPFM models.
- **Lattice initialization**: in the modular and discordant models, initial phenotypes (or latent states $y$) are snapped to the nearest multiple of $\delta$ and clamped to $\leq 0$ at initialization, reflecting the discrete binary locus structure. All subsequent mutations shift trait values by exactly $\pm\delta$, preserving the lattice structure throughout.
- **Output filenames** encode key simulation parameters (e.g., `ModularFGM_SSWM_N1e+04_M1e-07_d0.10_eR1.41_s2.00_L200-200.mat`).

---

## Citation

If you use this code, please cite:

> Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025). Module-Selection Balance in the Evolution of Modular Organisms. *Submitted.*

A BibTeX entry will be provided here once the paper is published.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

Minkyu Kim  
Department of Computational Biology, Cornell University  
mk2687@cornell.edu
