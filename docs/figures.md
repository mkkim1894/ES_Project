# Figure manifest

Which script produces which manuscript figure, and where the file lands. The
manuscript's `main_fig/` and `Supplementary_fig/` directories are maintained
separately; the paths below are the sources to copy from.

## The two kinds of script

`run_scripts/` holds both, and the distinction matters because one costs hours
and the other costs seconds.

**Simulators** write `.mat` files into `results/`. Each takes `'reproduce'`
(paper scale) or `'test'` (fast, every code path, minutes).

| script | writes | feeds |
| --- | --- | --- |
| `Run_pleiotropicFGM` | `results/{SSWM,CM_Asexual,CM_Sexual}/PleiotropicFGM_*` | F2, F3, F4 |
| `Run_modularFGM` | `results/{SSWM,CM_Asexual,CM_Sexual}/ModularFGM_*` | F2, F3, F4, S3 |
| `Run_nestedFGM` | `results/Generalization/NestedFGM/<REGIME>/` | F5, S4, S5 |
| `Run_modularConstantSupply` | `results/ConstantSupply/<REGIME>/` | F6 |
| `Run_pleiotropicRestrictedTheta` | `results/RestrictedTheta_<init>init/<REGIME>/` | restricted cone |
| `Run_supplementary` | `results/Supplementary/SteadyStateCM_N_*.mat` | S1, S2 |
| `Run_nModule` | `results/nModule/` | S7, S8 |
| `LTEE_analysis/ltee_analysis.py` | `LTEE_analysis/results/` | F7, S6 |

`Run_maxEntSimulations` is a wrapper, not a simulator: it calls
`Run_pleiotropicRestrictedTheta` and `Run_pleiotropicFGM` with
`'init', 'maxent'`, so any difference in the output is attributable to the
initial genotypes alone. It is a script driven by workspace variables
(`MODE`, `MODELS`, `STAGES`, `SEED`), not a function.

**Figure producers** read those `.mat` files and write `.pdf`. They run no
simulations, so they are cheap to repeat and are what most reruns need.

| script | writes |
| --- | --- |
| `Run_allMainFigures` | F2-F6 and the restricted cone, collected into `results/main_figures/` |
| `Run_mainFigures` | F2-F5 into `results/Figures/` |
| `Run_supplementaryFigures` | S1-S5 into `results/Supplementary/Figures/` |
| `Run_ThresholdDAnalysis` | `results/ThresholdD_Analysis.mat`, which S3 needs |
| `Run_nModule('stages', {'figures'})` | S7, S8 into `results/Figures/` |

`Run_ThresholdDAnalysis` sits between the two: it runs no simulation but it
writes a `.mat`, so S3 needs it before `Run_supplementaryFigures`.

## Main figures

| ms | data from | drawn by | output |
| --- | --- | --- | --- |
| F1 | — | drawn by hand | — |
| F2 | `Run_pleiotropicFGM`, `Run_modularFGM` | `Run_mainFigures`, figure 2 | `results/Figures/Figure2_SSWM.pdf` |
| F3 | same | `Run_mainFigures`, figure 3 | `results/Figures/Figure3_CMLinked.pdf` |
| F4 | same | `Run_mainFigures`, figure 4 | `results/Figures/Figure4_CMUnlinked.pdf` |
| F5 | `Run_nestedFGM` | `Run_mainFigures`, figure 5 | `results/Figures/Figure_NestedFGM_Generations.pdf` |
| F6 | `Run_modularConstantSupply` | same, `'stages', {'figures'}` | `results/Figures/Figure_ModularConstantSupply_Generations.pdf` |
| F7 | LTEE metagenomics | `ltee_analysis.py` | `LTEE_analysis/results/F7.pdf` |

The restricted-cone figure is drawn by `makeFigure_RestrictedTheta` into
`results/Figures/Figure_RestrictedTheta_Generations.pdf`. Its `init` argument
selects the results subtree, `sampled` or `maxent`; the parameter still defaults
to `sampled` for backward compatibility, but `Run_allMainFigures` passes
`maxent`, which is the published initializer. Its top row is the distribution of
phenotypic effects at three of the six initial conditions, computed by
`computeDPEAtPhenotype` and drawn by `drawDPEPanel` — no simulation is involved
in that row.

## Supplementary figures

| ms | data from | drawn by | output |
| --- | --- | --- | --- |
| S1 | `Run_supplementary` | `Run_supplementaryFigures`, 1 | `results/Supplementary/Figures/FigureS_CM_Main.pdf` |
| S2 | `Run_supplementary` | `Run_supplementaryFigures`, 2 | `results/Supplementary/Figures/FigureS_CM_Weights.pdf` |
| S3 | `Run_modularFGM` + `Run_ThresholdDAnalysis` | `Run_supplementaryFigures`, 3 | `results/Supplementary/Figures/FigureS_ThresholdDTrajectories.pdf` |
| S4 | `Run_nestedFGM` | `Run_supplementaryFigures`, 4 | `results/Supplementary/Figures/FigureS_NestedFGMDistributions.pdf` |
| S5 | `Run_nestedFGM(..., [10,20])` | `Run_supplementaryFigures`, 5 | `results/Supplementary/Figures/FigureS_NestedFGM_Asymmetric.pdf` |
| S6 | LTEE metagenomics | `ltee_analysis.py` | `LTEE_analysis/results/FigureS_LTEE_multihit.pdf` |
| S7 | `Run_nModule` | `Run_nModule`, figures stage | `results/Figures/Figure_nModule.pdf` |
| S8 | `Run_nModule` | same | `results/Figures/Figure_nModule_Reference.pdf` |

Panels C and D of S1 are split at `v'_2/v'_1 = 100`, the threshold `D` of the
heuristic, so that the panel titles `> 100` and `<= 100` are literally true.
`makeFigureS_SteadyStateCM` encodes that split directly; the earlier variant that
put the ratio 95.6 in panel C is in `_relegated/published_grid_duplicate/`.

`Run_supplementaryFigures` numbers S1 to S5 only. The ten-module figures are
produced by `Run_nModule` and the LTEE figures by Python, so neither passes
through it. `results/Figures/Figure_nModule_RarefactionCheck.pdf` is a diagnostic,
not a manuscript figure.

## Order to run from nothing

```matlab
cd ~/ES_Project
reproduce_all            % every simulation and every MATLAB figure, ~10 h
```

`reproduce_all` is in the project root and sets up its own path. Every other
entry point below is a function in one of the subdirectories, so a MATLAB
session needs `addpath(genpath(pwd));` once before calling them.

Or, to redraw figures from the simulation output already on disk:

```matlab
cd ~/ES_Project
addpath(genpath(pwd));
Run_allMainFigures                     % F2-F6 and the restricted cone
Run_supplementaryFigures('reproduce')  % S1-S5
Run_nModule('stages', {'figures'})     % S7, S8
```

## Theory curves

Trait-space panels carry the analytical prediction in every regime. Log-ratio
panels carry it only in the successive-mutations regime, where a closed form
exists for both the pleiotropic and the modular model. Under concurrent
mutations a time-resolved prediction exists for the modular model but not for
the pleiotropic one, so drawing it would fill one column of the row and leave
the other permanently empty. `predictLogRatioTrajectory` still returns the
concurrent-mutations curves; the figure code declines to draw them.

The constant-supply figure carries no theory on its log-ratio row for the same
reason: of its three regimes only the successive-mutations prediction is
available in closed form. Its closed form is documented in the header of
`predictModularSSWM_ConstantSupply`; the manuscript section it is cited from
derives the declining-supply model instead.
