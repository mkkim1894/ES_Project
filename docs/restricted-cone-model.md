# Restricted-cone pleiotropic GPFM

The fourth cell of the stress test on the necessary conditions for module-selection balance.

## Why this model exists

Module-selection balance might be a property of a **changing mutational supply** rather than of
**variational modularity**. Distinguishing the two requires crossing the factors rather than arguing
about which one matters.

The supply in the modular GPFM declines because of assumption (ii) in *Necessary conditions for
module-selection balance*: the optimal trait value is realized by a single genetic sequence. Low
genotypic redundancy near the optimum is what makes `b_i = |x_i|/δ` shrink as module *i* improves.

|  | **low redundancy** (supply declines) | **high redundancy** (supply constant) |
|---|---|---|
| **modular GPM** | main model → MSB | constant-supply model → MSB only under complete linkage |
| **pleiotropic GPM** | **this model → ?** | standard pleiotropic GPFM → no MSB |

The constant-supply model relaxes assumption (ii) while keeping modularity. This model tightens
assumption (ii) while keeping universal pleiotropy. With both cells filled, the alternative
explanation is tested directly rather than argued around.

## What restricting the cone actually does

Each locus carries a fixed pleiotropic angle θ_ℓ. The genotype at phenotype **x** satisfies

```
sum_{g=1} delta * cos(theta) = |x_1|
sum_{g=1} delta * sin(theta) = |x_2|
```

Under **isotropic** θ these sums contain terms of both signs, so **cancellation is possible**: a
near-optimal phenotype can be reached by hundreds of large opposing contributions. Redundancy near the
optimum stays high, the number of available beneficial mutations does not track |x_i|, and the supply
never declines. That is why the standard pleiotropic GPFM behaves like gradient ascent.

Confine θ to a single quadrant and every term becomes non-negative. **Cancellation is impossible.** The
number of allele-1 loci is tightly bounded by |x|, few genotypes map to near-optimal phenotypes, and as
x₂ approaches its optimum the surviving allele-1 loci are forced toward θ ≈ 0 — that is, toward being
x₁-improving. The redundancy of the optimum collapses and the supply declines per axis, with no modular
encoding anywhere.

Crucially, **every mutation still moves both traits.** The cone constrains the sign pattern, not the
pleiotropy: a quadrant is symmetric about θ = π/4, so no mutation is axis-aligned and none affects only
one trait. What it removes is antagonistic mutations — those that improve one trait while degrading the
other.

`diagnoseRestrictedTheta` confirms the mechanism is present before any simulation is run. At the
quadrant cone the beneficial supply falls from 24 available mutations to 4 as the traits improve, and
the mean displacement tilts toward x₁. Widen the cone to a half-circle and cancellation returns: the
supply only falls from 114 to 67, and the fraction of beneficial mutations that improve x₁ drops from
0.14 to 0.045 because most now trade one trait against the other. The cone width is therefore an
internal control on redundancy, not just a robustness knob.

## This is a property of the map, not of the starting point

θ_ℓ is drawn once per locus and fixed. There is no separate mutation process to restrict — the direction
a mutation moves the population is whatever angle its locus was assigned. So the restriction is
instantiated when the map is built, and it **persists for the entire run**: every mutation available at
every generation comes from the same restricted set of directions.

This deserves an explicit statement in the Methods, because "we restrict the angles at initialization" can be
misread as a transient starting-point effect. It is not: it is a permanent property of the
genotype–phenotype map.

## The angle range: [π, 3π/2] and [0, π/2] are the same model

Equation (pleiotropic GPM) in the manuscript writes `x_i = +delta * sum g_l cos(theta_l)`, so allele 1
displaces the phenotype by `+delta*(cos, sin)` and the third quadrant requires **θ ∈ [π, 3π/2]**. This
codebase uses the opposite sign — allele 1 displaces by `-delta*(cos, sin)` — so the identical set of
directions is **θ ∈ [0, π/2]**. Both give uniformly distributed unit vectors in the third quadrant.
Nothing differs scientifically; only the stored sign of θ.

Note that this is a genuine inconsistency inside the manuscript: Eq. (pleiotropic GPM) uses +δ while
Eq. (modular GPM) uses −δ, and the code follows the modular one. It is unobservable while θ is
isotropic, and load-bearing here. Fixing the sign in Eq. (pleiotropic GPM) is a one-character change
with no downstream consequences — the expression for s, the θ₀ half-plane, and the uniqueness of
`00...0` as the optimum all hold either way.

`initializeGenomeThetaRestricted` raises an explicit error if the cone points the wrong way, rather than
producing all-zero genomes and six populations sitting at the optimum while the run completes normally.

## One deliberate deviation from `initializeGenomeTheta`

The main-tree initializer makes a **single forward pass** over the loci and only ever flips 0 → 1. With
angles over the full circle that is fine — some direction always helps. Inside a cone it is not, and the
failure is silent.

Measured on the six paper initial conditions at the quadrant cone, the single-pass initializer leaves
residuals up to **1.17** and compresses the realized spread of initial ratios from a requested log R₀
range of **3.47 down to 1.19**. Populations would have started two-thirds converged, and any convergence
in the results would have been manufactured before generation 1.

`initializeGenomeThetaRestricted` iterates the greedy pass to convergence and allows flips in both
directions — ordinary coordinate descent on the distance to the target. Worst residual becomes **3.1% of
the target norm** and **3.31 of the 3.47** spread is preserved. Both numbers are printed on every run,
with a warning if the spread is compressed.

## Files

| file | what it does |
|---|---|
| `diagnoseRestrictedTheta.m` | Pre-flight, seconds, no simulation. Checks the cone reaches the initial conditions, then walks a genome toward the optimum and reports how the available supply changes. |
| `initializeGenomeThetaRestricted.m` | `initializeGenomeTheta` with a `thetaRange`, the iterated initializer, the redundancy diagnostics, and the wrong-cone guard. |
| `Run_pleiotropicRestrictedTheta.m` | Driver. Stages: `diagnose`, `sswm`, `cm_asexual`, `cm_sexual`, `figures`. |
| `makeFigure_RestrictedTheta.m` | 2×3 figure, A–C trait space and D–F log R, matching Figures 2–4 and 6. Regimes not yet run are drawn as empty labelled panels. |

Nothing here modifies the main tree, but it does need `simulation_scripts/`, `analysis_scripts/` and
`utils/` from the parent directory, which the driver locates automatically.

`simulatePleiotropicSSWM` calls **`mybinornd`**, which is not in the repository. Put the folder
containing it on the MATLAB path before the SSWM stage; the driver checks and errors out early. The CM
stages do not need it.

## Running it

```matlab
% 1. Seconds, no simulation.
diagnoseRestrictedTheta

% 2. Sanity run. Prints wall-clock time per stage.
Run_pleiotropicRestrictedTheta('mode', 'test')

% 3. Paper scale, SSWM. This alone answers the question.
Run_pleiotropicRestrictedTheta('stages', {'sswm', 'figures'})

% 4. Only if a balance appears: redundancy control, then the other regimes.
Run_pleiotropicRestrictedTheta('coneHalfWidth', pi/8, 'stages', {'sswm'})
Run_pleiotropicRestrictedTheta('coneHalfWidth', pi/4, 'stages', {'sswm'})
Run_pleiotropicRestrictedTheta('stages', {'cm_asexual', 'cm_sexual', 'figures'})
```

Test mode is 8 replicates × 2 conditions = 16 runs; paper scale is 250 × 6 = 1500. **Multiply the
printed test-mode time by ~94.** SSWM is cheap for a structural reason: inside a cone every allele-1
locus is beneficial and there are only ~25–35 of them, so each replicate fixes at most that many
mutations. The CM stages are far more expensive.

**All three regimes are needed**, not just SSWM. The constant-supply model's answer turned out to be
regime-dependent — no balance under successive mutations or free reassortment, but a balance under
complete linkage, where clonal interference substitutes for the declining supply. The 2×2 above is
really a 2×2×3, and this cell needs all three entries or the comparison is incomplete in exactly the
dimension that mattered.

SSWM is run first only as sequencing: it is minutes rather than hours, and without clonal interference
it isolates the redundancy mechanism from the selection-coupling mechanism that the constant-supply
model already showed can produce a balance by itself. A balance appearing here is the
conclusion-changing outcome, so it is checked before committing to the concurrent-mutations runs.

Parameters are identical to `Run_pleiotropicFGM`: N = 10⁴, δ = 0.1, 2L = 400, ã₁ = 1, ã₂ = 1/√2, σ = 2,
250 replicates, six initial conditions on the W₀ = 0.25 contour, U = 2×10⁻⁶ (SSWM) and 4×10⁻³ (CM). The
angle distribution is the only difference.

## Reading the result

Panels D–F carry the answer. In A–C trajectories crowd together near the optimum for reasons that have
nothing to do with an attractor.

- **Trajectories hug the yellow curves** (the unrestricted isotropic FGM gradient path): reducing
  redundancy alone did not produce a balance. Modularity is doing the work, the claim stands as written,
  and this becomes a supplementary control that closes the alternative.
- **Trajectories collapse onto the orange line** (log a₂²/a₁²): a module-selection balance from low
  genotypic redundancy alone, without modular encoding.
- **Neither**: the cone changes the path without producing an attractor, with the six curves staying
  ordered by initial condition. Report it as such.

## If a balance appears

The revision is a narrowing, not a retraction. The operative condition becomes **low genotypic
redundancy of the trait optimum**, of which variational modularity is the biologically documented
realization rather than the only conceivable one. That reading already fits the rest of the section: the
nested FGM keeps the balance under a non-linear decline in supply, and the constant-supply model loses it
except under complete linkage, where clonal interference supplies the coupling instead. The cone-width
sweep then becomes the quantitative statement — how much redundancy the balance tolerates before it
disappears.
