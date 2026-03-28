# Module-Selection Balance in the Evolution of Modular Organisms

**Authors:** Minkyu Kim, Sarah M. Ardell, Sergey Kryazhimskiy  
**Affiliations:** Cornell University; University of California San Diego

---

## Overview

This repository contains MATLAB code and a Python notebook for the simulations and analyses in:

> Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025). Module-Selection Balance in the Evolution of Modular Organisms.

Four genotype-phenotype-fitness map (GPFM) models are implemented — Pleiotropic, Modular, Discordant-module, and Nested FGM — each simulated under the Strong Selection Weak Mutation (SSWM) and Concurrent Mutations (CM) regimes. A separate Python notebook reproduces the LTEE metagenomic analysis (Figure 6).

---

## Repository Structure

```
project_root/
├── run_scripts/
│   ├── Run_pleiotropicFGM.m          # Pleiotropic FGM simulations (Figure 2)
│   ├── Run_modularFGM.m              # Modular FGM simulations (Figure 3)
│   ├── Run_discordantFGM.m           # Discordant-module FGM simulations (Figure 4)
│   ├── Run_nestedFGM.m               # Nested FGM simulations (Figure 5, S5, S6)
│   ├── Run_supplementary.m           # Steady-state CM validation (Figures S1–S3)
│   ├── Run_ThresholdDAnalysis.m      # Threshold D sensitivity analysis (Figure S4)
│   ├── Run_mainFigures.m             # Generate main figures 2–5
│   └── Run_supplementaryFigures.m    # Generate supplementary figures S1–S6
├── simulation_scripts/
│   ├── simulatePleiotropicSSWM.m
│   ├── simulatePleiotropicCM.m
│   ├── simulateModularSSWM.m
│   ├── simulateModularCM.m
│   ├── simulateDiscordantSSWM.m
│   ├── simulateDiscordantCM.m
│   ├── simulateNestedSSWM.m
│   ├── simulateNestedCM.m
│   └── simulateSteadyStateCM.m
├── analysis_scripts/
│   ├── computeAverageTrajectory.m
│   ├── predictModularSSWM.m
│   ├── predictModularCM.m
│   ├── predictFullRecomb.m
│   └── predictPleiotropicSSWM.m
├── utils/
│   ├── initializeSimParams.m
│   ├── findInitialPhenotypes.m
│   ├── freezeParam.m
│   ├── preRunSimulation.m
│   └── initializeGenomeTheta.m
├── figure_scripts/
│   ├── makeFigure2_Generations.m
│   ├── makeFigure3_Generations.m
│   ├── makeFigure4_Generations.m
│   ├── makeFigure5_Generations.m
│   ├── makeFigureS_SteadyStateCM.m
│   ├── makeFigureS_ThresholdDTrajectories.m
│   └── makeFigureS_NestedFGMDistributions.m
├── LTEE_analysis/
│   └── notebooks/
│       └── ltee_analysis_revised.ipynb   # LTEE metagenomic analysis (Figure 6)
├── reproduce_all.m               # Reproduce all simulations and figures (~10 hours, 10 cores)
├── test_all.m                    # Pipeline verification (~15 minutes)
├── README.md
└── LICENSE
```

---

## Requirements

### MATLAB (Figures 2–5, S1–S6)
- MATLAB R2025a

| Toolbox | Purpose |
|---------|---------|
| Statistics and Machine Learning | `mnrnd`, `poissrnd`, `datasample` |
| Symbolic Math | `syms`, `solve` in `findInitialPhenotypes` |
| Parallel Computing | `parfor` acceleration |

### Python (Figure 6)
- Python 3.8 or later
- Dependencies: `numpy`, `pandas`, `matplotlib`
- LTEE metagenomic data (see below)

---

## LTEE Data Setup (Figure 6)

The LTEE analysis notebook requires data from Good et al. (2017), which is not included in this repository.

1. Download the repository ZIP from https://github.com/benjaminhgood/LTEE-metagenomic
2. Unzip and place the resulting `LTEE-metagenomic-master` folder at `LTEE_analysis/LTEE-metagenomic-master/` (one level above the notebook).

The expected directory layout is:

```
LTEE_analysis/
├── LTEE-metagenomic-master/   ← place downloaded data here
│   └── data_files/
└── notebooks/
    └── ltee_analysis_revised.ipynb
```

Alternatively, update `LTEE_REPO_PATH` in the Setup cell of the notebook to point to your local copy.

---

## Reproducing Paper Results

To verify the pipeline before a full run (~15 minutes, reduced parameter sets):

```matlab
cd /path/to/project
test_all
```

To reproduce all simulations and figures (~10 hours on 10 cores):

```matlab
cd /path/to/project
reproduce_all
```

Output files are written to `results/`, `results_supplementary/`, and their respective `Figures/` subdirectories.

---

## Running Individual Models

Each driver script accepts `'reproduce'` (default) or `'test'` as the mode argument.

```matlab
Run_pleiotropicFGM('reproduce')           % Figure 2
Run_modularFGM('reproduce')               % Figure 3
Run_discordantFGM('reproduce')            % Figure 4
Run_nestedFGM('reproduce')               % Figure 5, S5  — symmetric [n1=10, n2=10]
Run_nestedFGM('reproduce', {}, [10, 20]) % Figure S6     — asymmetric [n1=10, n2=20]
Run_supplementary('reproduce')           % Figures S1–S3
Run_ThresholdDAnalysis('reproduce')      % Figure S4
```

After simulations complete, generate figures with:

```matlab
Run_mainFigures('reproduce')           % Figures 2–5
Run_supplementaryFigures('reproduce')  % Figures S1–S6
```

Individual figures can be regenerated without rerunning the others:

```matlab
Run_mainFigures('reproduce', 3)           % Figure 3 only
Run_supplementaryFigures('reproduce', 4)  % Figure S4 only
```

---

## Figure Index

### Main Figures

| Figure | Description | Source |
|--------|-------------|--------|
| 2 | Evolutionary dynamics on the pleiotropic GPFM | `Run_pleiotropicFGM` |
| 3 | Evolutionary dynamics on the modular GPFM | `Run_modularFGM` |
| 4 | Evolutionary dynamics on the discordant-module GPFM | `Run_discordantFGM` |
| 5 | Evolutionary dynamics on the nested FGM | `Run_nestedFGM` |
| 6 | LTEE metagenomic analysis | `ltee_analysis_revised.ipynb` |

### Supplementary Figures

| Figure | Description | Simulation required |
|--------|-------------|---------------------|
| S1 | Steady-state CM rate-of-adaptation: main validation | `Run_supplementary` |
| S2 | Steady-state CM: parameter ranges | `Run_supplementary` |
| S3 | Steady-state CM: weighting scheme comparisons | `Run_supplementary` |
| S4 | Threshold D sensitivity (D = 10, 100, 1000, 10000) | `Run_ThresholdDAnalysis` |
| S5 | Nested FGM mutation-effect distributions | `Run_nestedFGM` |
| S6 | Asymmetric nested FGM dynamics (n1=10, n2=20) | `Run_nestedFGM('reproduce', {}, [10,20])` |

---

## Key Parameters

| Parameter | Symbol | Default | Description |
|-----------|--------|---------|-------------|
| `popSize` | $N$ | $10^4$ | Population size |
| `mutationRate` | $U$ | $U = U_\mathrm{ref} \times (L/K)$ | Genome-wide mutation rate, scaled to keep the per-locus rate $\mu = U_\mathrm{ref}/(2K)$ fixed as $L$ changes. $U_\mathrm{ref} = 10^{-7}$ (SSWM), $2\times10^{-4}$ (CM) |
| `deltaTrait` | $\delta$ | 0.1 | Mutational step size |
| `landscapeStdDev` | $\sigma$ | 2 | Fitness landscape width |
| `ellipseRatio` | $a_1/a_2$ | $\sqrt{2}$ | Selection anisotropy |
| `geneticTargetSize` | $[L_1, L_2]$ | $[200, 200]$ | Number of loci per module; controls the deleterious mutation supply (upper boundary = $L\delta$) |
| `K` | — | 10 | Reference locus count at which $U_\mathrm{ref}$ was originally calibrated; defines the fixed per-locus rate $\mu = U_\mathrm{ref}/(2K) = 5\times10^{-9}$ (SSWM), $10^{-5}$ (CM) |
| `recombinationRate` | $\rho$ | 0 or 1 | Recombination rate |

Parameters are set in each `Run_*.m` script and passed to simulations via `initializeSimParams`.

### Implementation Notes

- **Beneficial mutation threshold**: all SSWM simulations use $s > 0$ as the threshold for a mutation to be considered beneficial, consistently across all four GPFM models.
- **Lattice initialization**: in the modular and discordant models, initial phenotypes (or latent states $y$) are snapped to the nearest multiple of $\delta$ and clamped to $\leq 0$ at initialization, reflecting the discrete binary locus structure. All subsequent mutations shift trait values by exactly $\pm\delta$, preserving the lattice structure throughout.
- **Per-locus mutation rate**: for the pleiotropic, modular, and discordant models, $\mu = U/(2L)$ is held fixed at $5\times10^{-9}$ (SSWM) and $10^{-5}$ (CM) regardless of $L$, by scaling $U = U_\mathrm{ref} \times (L/K)$. For the pleiotropic model, $L$ refers to the per-module locus count with total loci $= 2L = 400$. The nested FGM uses a continuous Gaussian mutation framework and is parameterized directly by $U$.
- **Output filenames** encode key simulation parameters including locus count as `L{value}` (e.g., `ModularFGM_SSWM_N1e+04_M2e-06_d0.10_eR1.41_s2.00_L200-200.mat`).

---

## Citation

```bibtex
@article{kim2025module,
  title={Module-Selection Balance in the Evolution of Modular Organisms},
  author={Kim, Minkyu and Ardell, Sarah M. and Kryazhimskiy, Sergey},
  year={2025}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

Minkyu Kim  
Department of Computational Biology, Cornell University  
mk2687@cornell.edu
