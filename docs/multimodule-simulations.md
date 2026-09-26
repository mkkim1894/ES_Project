# Ten-module simulations

The modular and pleiotropic GPFMs with more than two functional modules,
simulated in the concurrent-mutations regime with complete linkage. Produces
Supplementary Figures S7 and S8.

| | |
| --- | --- |
| entry point | `Run_nModule` |
| simulators | `runModularND`, `runPleiotropicND` |
| figure | `makeFigure_nModule` |
| output | `results/nModule` (simulations), `results/Figures` (figures) |

## Why this model exists

The two-module models establish the mechanism but cannot be compared directly
with genomic data, where the number of functional modules contributing to
fitness is larger. These simulations ask whether the distinction between the two
architectures survives in a higher-dimensional trait space, and whether the
biphasic pattern the modular model predicts — adaptation concentrated in a few
lagging modules, then broadening — is what the ten-module version shows.

## Summary statistics

Two effective-number statistics replace the module performance ratio, which has
no direct analogue beyond two modules.

- `n_lag`, the effective number of lagging modules, summarises how the distances
  of the traits from their optima are distributed. It is 1 when a single module
  is far behind and *n* when all are equally far.
- `n_sel`, the effective number of module targets, measures how broadly
  phenotypic improvement is distributed within a sliding window. It is 1 when a
  single module improves and *n* when all improve equally.

`n_sel` is computed the same way as the effective number of selected genes in
the LTEE analysis: observations are pooled across replicates, rarefied to a
common depth, and Simpson's index is averaged across subsamples before being
inverted. The two analyses are therefore directly comparable.

## Running it

```matlab
cd ~/ES_Project
addpath(genpath('.'));
Run_nModule                      % defaults: n = 10, 30 replicates, 1e4 generations
Run_nModule('n', 5)
```

Copyright (c) 2025 Minkyu Kim, Cornell University. MIT License.
