# Maximum-entropy initialization

Initial genotypes drawn from the maximum-entropy (microcanonical) ensemble at a
given phenotype, and the figure set rebuilt on top of it.

Method: maximum-entropy (microcanonical) sampling at a fixed phenotype, by
exponential tilting.

## Why

`utils/initializeGenomeSampled` builds an initial genotype by a random-order walk
that fixes the number of loci at allele 1, followed by constrained swaps. The
swap stage is a symmetric-proposal Metropolis chain with an indicator
acceptance, so it is uniform over genotypes at fixed *n*; the choice of *n* is a
heuristic with no guarantee. The maximum-entropy construction removes the
heuristic: it samples from the uniform distribution over **all** genotypes at the
target phenotype, with *n* unconstrained.

Plain rejection sampling cannot do this. With the sign convention of this
codebase a locus at allele 1 displaces the phenotype by
`-delta*(cos theta, sin theta)`, so writing `e_l = (cos theta_l, sin theta_l)`
and `y = -x/delta`, a genotype's phenotype is `y(g) = sum_l g_l e_l`. Inside a
cone of width `pi/2` the direction of that sum concentrates near 45 degrees with
a standard deviation of about 4.7 degrees at `n = 33`, so the `R_0 = 0.16` target
at 9.5 degrees is 7.6 standard deviations out — roughly `2e-14` per draw. That is
the obstacle.

## The construction

Tilt the ensemble:

```
P_lambda(g) = exp(<lambda, y(g)>) / Z(lambda)
```

Every genotype with the same `y` gets the same probability whatever `lambda`, so
conditioned on `y` this *is* the uniform distribution over genotypes at that
phenotype. `Z` factorises over loci, so under `P_lambda` the loci are independent
Bernoulli with `p_l = sigmoid(<lambda, e_l>)`. Choosing `lambda` to solve

```
sum_l sigmoid(<lambda, e_l>) e_l = y_0
```

makes the target the *mean* of the ensemble rather than a rare event. The
Jacobian `sum_l p_l(1-p_l) e_l e_l'` is positive definite, so the solution is
unique and Newton converges in a handful of steps. Acceptance then runs at
0.04–1.4 percent instead of `1e-14`.

### One addition to the note

`y` is a sum of `L` Bernoulli terms, so the set of genotypes hitting `y_0`
*exactly* is generally empty and one accepts a window `|y - y_0| <= tol/delta`.
Inside that window `P_lambda` is not flat — it varies as
`exp(<lambda, y - y_0>)`. Accepted draws are therefore importance weighted by
`exp(-<lambda, y - y_0>)` before one is chosen, which restores the uniform
ensemble over the window. At the default tolerance of 0.02 the correction spans a
factor of 3–27; at a tolerance of 0.06 it reaches `2e4` for targets near the cone
edge, which is why the default is tight. The diagnostic `weightRange` reports it
per condition.

## Files

| file | what it is |
| --- | --- |
| `utils/initializeGenomeMaxEnt.m` | the sampler; drop-in for `initializeGenomeSampled` (same four arguments, same output fields) |
| `utils/solveMaxEntTilt.m` | the Newton solve for the tilt `lambda` |
| `analysis_scripts/computeDPEAtPhenotype.m` | the distribution of available mutational directions at a phenotype, pooled over the ensemble |
| `run_scripts/Run_maxEntSimulations.m` | drives the three regimes with `'init', 'maxent'` |

The vector mean of the distribution of available directions is the direction to
the optimum for *any* genotype at that phenotype — it follows from
`sum_l delta e_l = -x` — so no two initializers can disagree about where the
distribution is centred. They can disagree about `n`, and the spread follows from
it, since the mean resultant length is pinned at `|x|/(delta*n)`.


## Running the simulations

`Run_maxEntSimulations.m` drives the three evolutionary regimes - SSWM, CM
asexual, CM sexual - with max-entropy initial genotypes. It does not contain a
simulator: it calls the same `Run_pleiotropicRestrictedTheta` and
`Run_pleiotropicFGM` that produce the published runs, with `'init', 'maxent'`,
so any difference in the results comes from the initial genotypes alone. Those
two runners each gained a third `init` case; their defaults are unchanged.

```matlab
cd ~/ES_Project
addpath(genpath('.'));
MODE = 'test';      Run_maxEntSimulations    % 8 replicates, minutes
MODE = 'reproduce'; Run_maxEntSimulations    % 250 replicates, hours
```

Defaults are `MODELS = {'restricted'}` and all three regimes; add
`MODELS = {'restricted', 'unrestricted'}` for the universal-pleiotropy control.
Run the test configuration first - it exercises every code path in minutes, and
each regime prints the initializer's per-condition line before simulating, so an
unreachable initial condition surfaces immediately rather than after hours.

Output filenames carry the initializer (`RestrictedThetaFGM-maxentinit_*`,
`PleiotropicFGM-maxentinit_*`), so nothing can overwrite an existing run.

### A pre-existing filename collision

In the main tree, `greedy` and `sampled` pleiotropic runs both write files named
`PleiotropicFGM_*`, so a sampled run overwrites a greedy one in place - which is
why `results/*/_relegated_greedyinit/` exists. `maxent` is tagged so it cannot do
the same. Giving `sampled` its own tag would be the tidier fix, but it would
orphan the existing September files, so it has been left alone.

### Feasibility at the production initial conditions

All six production initial conditions are reachable inside the cone at
`L_total = 400`, checked over six independent draws of the angles:

| x_2/x_1 | direction | reachable | E[n] | acceptance |
| --- | --- | --- | --- | --- |
| 0.156 | 8.9 deg | 6/6 | 33.2 | 0.80% |
| 0.313 | 17.3 deg | 6/6 | 33.5 | 0.27% |
| 0.625 | 32.0 deg | 6/6 | 32.9 | 0.18% |
| 1.25 | 51.3 deg | 6/6 | 29.8 | 0.18% |
| 2.5 | 68.2 deg | 6/6 | 26.3 | 0.27% |
| 5.0 | 78.7 deg | 6/6 | 24.3 | 0.55% |

The margin is thinnest at the two extremes, where the target direction
approaches the edge of the cone. At `L_total = 200` the most extreme condition
becomes marginal, and at 100 it is unreachable; `assertReachable` catches this
before Newton runs.

## What to expect

From the Python prototype of the same construction, at the three initial
conditions of the trajectory investigation:

| condition | `n` now | `n` max-ent | spread now | spread max-ent |
| --- | --- | --- | --- | --- |
| `R_0 = 0.16`, blue | 33 | 33.0 | 9.0 deg | 8.3 deg |
| `R_0 = 0.625`, orange | 32 | 32.8 | 23.4 deg | 26.4 deg |
| `R_0 = 5`, green | 28 | 25.0 | 29.5 deg | 10.8 deg |

Blue and orange should reproduce. Green is the informative case: the
maximum-entropy ensemble puts about three fewer loci at allele 1 there, which
tightens the cone of available directions from roughly 30 degrees to 11. Since
the selective shift away from the radial direction scales with the square of that
spread, the mutational-bias result should come out *stronger* at green, not
weaker. Note also that green is the condition whose current genotype landed 0.082
from its target, outside the 0.06 tolerance of the old sampler — the
maximum-entropy draws land inside 0.02 by construction.

For the unrestricted map the ensemble puts about 200 of the 400 loci at allele 1.
That has now been checked against the published pleiotropic SSWM result, whose
saved `genomeParams.initDiagnostics` records `nOnes = [148 293 212 230 140 178]`,
mean 200, from `initializeGenomeSampled` - so on the count the two agree, and
max-entropy mainly removes the 140-293 scatter. (`nOnesGreedy` in the same file
is `[48 48 44 37 39 41]`, which is what the relegated March runs used.) The
published Figures 2-4 therefore do not need redoing on account of the
initializer.

Copyright (c) 2025 Minkyu Kim, Cornell University. MIT License.
