# Module-Selection Balance in the Evolution of Modular Organisms

**Authors:** Minkyu Kim, Sarah M. Ardell, Sergey Kryazhimskiy  
**Affiliations:** Cornell University; University of California, San Diego

---

## Overview

This repository contains MATLAB code implementing Fisher's Geometric Model (FGM) variants for studying how organismal modularity affects evolutionary trajectories.

**Reference:**
> Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025). Module-Selection Balance in the Evolution of Modular Organisms. *Genetics* / *PNAS*.

---

## Models Implemented

| Model | Description |
|-------|-------------|
| **Modular FGM** | Mutations affect single modules with varying target sizes |
| **Pleiotropic FGM** | Mutations affect both traits with random pleiotropic angles |
| **Nested FGM** | Each module is itself a multi-dimensional FGM |
| **Standard FGM** | Isotropic mutations (baseline comparison) |

Each model is simulated under:
- **SSWM** (Strong Selection Weak Mutation) - Sequential fixation regime
- **CM** (Concurrent Mutations) - Multiple segregating mutations

---

## Repository Structure

```
project_root/
├── main/                     # Driver scripts
│   ├── Run_modularFGM.m      # Modular FGM (Figures 3-6)
│   ├── Run_pleiotropicFGM.m  # Pleiotropic FGM (Figure 2)
│   ├── Run_standardFGM.m     # Standard FGM baseline
│   ├── Run_supplementary.m   # All supplementary analyses
│   └── Run_ThresholdDAnalysis.m  # Threshold D sensitivity
├── simulation/               # Core simulation functions
│   ├── simulateModularSSWM.m
│   ├── simulateModularCM.m
│   ├── simulateNestedSSWM.m
│   ├── simulateNestedCM.m
│   ├── simulatePleiotropicSSWM.m
│   ├── simulatePleiotropicCM.m
│   ├── simulateStandardSSWM.m
│   ├── simulateStandardCM.m
│   ├── simulateSteadyStateCM.m      # CM steady-state validation
│   └── simulateSteadyStateRecomb.m  # Recombination variance
├── analysis/                 # Analysis and prediction functions
│   ├── computeAverageTrajectory.m
│   ├── predictModularSSWM.m
│   ├── predictModularCM.m
│   ├── predictFullRecomb.m
│   ├── predictPleiotropicSSWM.m
│   └── predictCM_ThresholdD.m   # Threshold D predictions
├── utils/                    # Utility functions
│   ├── initializeSimParams.m
│   ├── findInitialPhenotypes.m
│   ├── freezeParam.m
│   ├── preRunSimulation.m
│   └── initializeGenomeTheta.m
├── figures/                  # Figure generation scripts
│   ├── makeFigure2.m - makeFigure6.m
│   ├── makeFigure_StandardFGM.m
│   ├── makeFigureS_SteadyStateCM.m       # Supp: CM validation
│   ├── makeFigureS_SteadyStateRecomb.m   # Supp: Variance
│   ├── makeFigureS_NestedFGMDistributions.m  # Supp: Distributions
│   └── makeFigureS_ThresholdDTrajectories.m  # Supp: D threshold
├── setup.m                   # Path configuration
├── README.md
└── LICENSE
```

---

## Requirements

### Software
- **MATLAB R2020a** or later (developed on R2025a)

### Required Toolboxes
| Toolbox | Required | Purpose |
|---------|----------|---------|
| Statistics and Machine Learning | Yes | `mnrnd`, `poissrnd`, `datasample` |
| Symbolic Math | Yes | `solve`, `syms` in predictions |
| Parallel Computing | Optional | `parfor` acceleration |

### System Requirements
- Memory: 8 GB minimum, 16+ GB recommended
- Disk: ~1 GB for results

---

## Quick Start

```matlab
% 1. Setup paths
cd /path/to/project
setup

% 2. Run demo (5-10 minutes)
Run_modularFGM('demo')

% 3. Run full simulations (hours)
Run_modularFGM('full')
```

---

## Reproducing Paper Figures

| Figure | Command |
|--------|---------|
| Figure 2 | `Run_pleiotropicFGM('full')` |
| Figures 3-6 | `Run_modularFGM('full')` |

### Supplementary Figures

| Figure | Description | Command |
|--------|-------------|---------|
| Fig S1-S6 | CM steady-state validation | `Run_supplementary('full')` |
| Fig S_Recomb | Variance at recombination balance | `Run_supplementary('full')` |
| Fig S_Nested | Nested FGM distributions | `Run_supplementary('full')` |
| Fig S_ThresholdD | Threshold D sensitivity | `Run_supplementary('full')` |

To generate all supplementary figures:
```matlab
Run_supplementary('full')    % Full analysis
Run_supplementary('demo')    % Quick demo
Run_supplementary('figures') % Figures only (requires existing data)
```

---

## Parameters

| Parameter | Symbol | Default | Description |
|-----------|--------|---------|-------------|
| `popSize` | N | 10^4 | Population size |
| `mutationRate` | μ | 10^-7 (SSWM) / 2×10^-4 (CM) | Mutation rate |
| `deltaTrait` | δ | 0.1 | Mutational step size |
| `landscapeStdDev` | σ | 2 | Fitness landscape width |
| `ellipseRatio` | a₁/a₂ | √2 | Selection anisotropy |
| `geneticTargetSize` | [K₁,K₂] | [10, 10] | Loci per module |
| `recombinationRate` | ρ | 0-1 | Recombination rate |

---

## Citation

```bibtex
@article{kim2025module,
  title={Module-Selection Balance in the Evolution of Modular Organisms},
  author={Kim, Minkyu and Ardell, Sarah M. and Kryazhimskiy, Sergey},
  journal={Genetics},
  year={2025},
  doi={10.XXXX/genetics.XXX.XXXXXX}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) file.

---

## Contact

**Minkyu Kim**  
Department of Computational Biology, Cornell University  
Email: mk2687@cornell.edu  
GitHub: [@mkkim1894](https://github.com/mkkim1894)
