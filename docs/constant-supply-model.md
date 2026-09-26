# Constant-supply control model

Self-contained extension to ES_Project: the modular GPFM in which the supply of
module-improving mutations does **not** decline as a module approaches its
optimum. Control model, simulated in all three evolutionary regimes, with a figure laid out exactly like the nested-FGM
figure (`makeFigure5_Generations`).

Its files are laid out like every other model: runner in `run_scripts/`,
simulators in `simulation_scripts/`, predictions in `analysis_scripts/`,
figure in `figure_scripts/`, output under `results/ConstantSupply/`.

---

## The one change from the modular GPFM of Figure 3

Per-locus mutation rate is `mu_i = U / (2*L_i)` in both models. What differs is how
many loci of module *i* are module-improving.

| | improving loci `b_i` | beneficial supply `U_i` (per genome per generation) |
|---|---|---|
| **Figure 3** (`simulateModularSSWM`) | `b_i = \|x_i\| / delta` | `U_i = U*\|x_i\| / (2*delta*L_i)` — declines to zero at the optimum |
| **This model** | `b_i = f_i * L_i` | `U_i = U * f_i / 2` — constant |

`f_1` and `f_2` are fixed beneficial fractions. They differ between the two modules
but do **not** vary across the trait space, across initial conditions, or over the
course of evolution — they are properties of the genotype-phenotype map. Note that
`L_i` cancels out of `U_i`, so the supply depends on the fraction alone.

**The supply is constant everywhere in trait space and is never switched off** —
not at the optimum, not beyond it. There is no test on `x_i` anywhere in the
mutation rates, in either regime. What stops a module at its optimum is
**selection**: a `+delta` step at `x_i = 0` has `s < 0`, so in SSWM the Kimura
factor rejects it (`Pr_fix ≈ 3e-14` at N = 10⁴) and in CM Wright-Fisher sampling
purges it. This matters — a supply that switched off at the optimum would be a
performance dependence, which is the very thing this control model exists to
remove, and it would also misrepresent the argument being tested: under constant
evolvability the *rate* of evolution decays near the optimum because the selection
gradient decays, not because the variance runs out.

This framing deliberately severs the locus bookkeeping of the declining model. In
that model `|x_i| = b_i·δ` **by definition**, so `b_i` cannot be held fixed — a
fixed `b_i` would pin `|x_i|` to a single value. Severing that link *is* the
assumption under test, so it is better stated as an assumption about genotypic
redundancy than as a literal count of loci.

In the CM simulations the complementary fraction `1 - f_i` of loci is deleterious,
so the deleterious supply is likewise constant at `U*(1-f_i)/2` — also with no
dependence on `x_i`. In the declining model it is `mu_i*(L_i*delta + x_i)/delta`,
which at the starting phenotypes is numerically very close (≈0.43 U vs 0.48 U per
module at `f = 0.05`), so the deleterious load is not what drives the difference.

Everything else — population size, landscape, anisotropy, step size, genetic target
sizes, mutation rates, recombination, initial conditions, Kimura fixation,
Wright-Fisher resampling — is unchanged. Any difference in outcome is attributable
to the mutational supply alone.

---

## Why this destroys module-selection balance

With `s_i ~= |x_i|*delta / (sigma^2*a_i^2)` and `P_fix ~= 2*s_i`, the per-module
fixation flux is `N * U_i * 2*s_i * delta`.

**Declining supply** (Figure 3): `U_i ∝ |x_i|`, so the flux is **quadratic** in
`|x_i|`:

```
dx_i/dt = alpha_i * x_i^2        x_i(t) = x_i0 / (1 - alpha_i*x_i0*t) ~ -1/(alpha_i*t)
```

Both traits decay as a *power law*, so `x_2/x_1 -> alpha_1/alpha_2 =
(L_2*a_2^2)/(L_1*a_1^2)` — a genuine attractor, independent of the initial
condition, at which `s_1 = s_2`. That is module-selection balance.

**Constant supply** (this model): `U_i` is constant, so the flux is **linear** in
`|x_i|`:

```
dx_i/dt = -beta_i * x_i          beta_i = N*U*f_i*delta^2 / (sigma^2*a_i^2)
x_i(t)  = x_i0 * exp(-beta_i*t)
```

Both traits decay *exponentially*, so

```
log( x_2(t)/x_1(t) ) = log( x_2(0)/x_1(0) ) - (beta_2 - beta_1)*t
beta_2/beta_1 = (f_2/f_1) * (a_1^2/a_2^2)
```

The log ratio drifts linearly and without bound, to `0` or to `Inf`, and never
forgets its initial condition. There is no attractor, which is the predicted
behaviour for a supply that does not decline.

**Knife-edge case.** If `f_2/f_1 = a_2^2/a_1^2` exactly, then
`beta_1 = beta_2` and the ratio is frozen at its initial value. That still is *not*
a balance: each initial condition keeps its own ratio, so the trajectories stay a
family of distinct rays instead of collapsing onto one. The contrast is therefore
about the existence of an attractor, not about whether the ratio happens to change. `Run_modularConstantSupply` prints a
note if you pick fractions that land on this knife-edge.

The default `f = [0.05, 0.10]` gives `beta_2/beta_1 = 2 * 2 = 4`, so every initial
condition drifts steadily downward and the effect is unambiguous.

---

## Contents

```
results/ConstantSupply/
├── Run_modularConstantSupply.m            driver: all three regimes + figure
├── simulateModularSSWM_ConstantSupply.m   SSWM simulation
├── simulateModularCM_ConstantSupply.m     CM simulation (rho = 0 and rho = 1)
├── predictModularSSWM_ConstantSupply.m    closed form, exponential decay
├── predictModularCM_ConstantSupply.m      Desai-Fisher; one line differs from predictModularCM
├── predictFullRecomb_ConstantSupply.m     uncoupled modules, integrated numerically
├── makeFigure_ConstantSupply.m            3x2 figure, same layout as Figure 5
├── README.md
└── results/                               created on first run
    ├── ConstantSupply/{SSWM, CM_Asexual, CM_Sexual}/*.mat
    └── Figures/Figure_ModularConstantSupply_Generations.pdf
```

`predictFullRecomb_ConstantSupply` integrates the two uncoupled module ODEs
numerically rather than using a closed form. The closed form in
`predictFullRecomb` relies on `s_i/U_i` being nearly constant along the
trajectory, which holds when `U_i ∝ |x_i|` but fails when `U_i` is constant —
there `s_i/U_i` varies by orders of magnitude. The two ODEs are uncoupled under
full reassortment, so integrating them is cheap and involves no extra
approximation.

### Dependencies on the parent project

Deliberately **not** copied here, so that this model and the modular GPFM of
Figure 3 are guaranteed to share the same code:

- `../utils/initializeSimParams.m`, `../utils/findInitialPhenotypes.m`
- `../analysis_scripts/computeAverageTrajectory.m`

`Run_modularConstantSupply` adds those two folders to the path itself and errors
early with a clear message if anything is missing. Constant-supply settings are
attached to `simParams` *after* `initializeSimParams` returns, so that shared
utility needs no modification.

---

## Running

```matlab
cd ~/ES_Project
addpath(genpath('.'));

Run_modularConstantSupply('mode', 'test')          % fast sanity run -> results_test/
Run_modularConstantSupply                          % full run, f = [0.05, 0.10]
Run_modularConstantSupply('beneficialFraction', [0.05 0.10])
Run_modularConstantSupply('stages', {'figures'})   % redraw without rerunning
```

Runtime should be the same order as the corresponding `Run_modularFGM` runs — SSWM
is quick, the two CM regimes dominate.

The SSWM stage prints the predicted `beta_1`, `beta_2` and the drift slope
`-(beta_2 - beta_1)`, so you can see the expected direction of divergence before
looking at the figure.

---

## Figure

Same 3x2 layout as `makeFigure5_Generations`:

- **A, B, C** — phenotypic trajectories in the three regimes. Blue: simulation
  averages. Yellow: analytical prediction. Dashed orange: the `s_1 = s_2` ray of
  the *declining*-supply model, drawn as the reference the trajectories are **not**
  expected to approach.
- **D, E, F** — `log(x_2/x_1)` vs generations, mean ± s.d. across replicates.
  Dashed orange: the balance value `log[(L_2 a_2^2)/(L_1 a_1^2)]` of the declining
  model. Panels D–F use `ylim = [-3, 2]`, as in Figures 3 and 5.

`'proximityCutoff'` defaults to `[0, 0.25, 0.25]` — no cutoff on SSWM, 0.25 on the
two CM panels — and accepts a scalar or a 3-vector. Rename the output with
`'outputFile'` once you settle the main-text figure number.

### Suggested caption wording for the cutoff

> Panels E and F show the phase in which both modules are still adapting
> (trajectories are truncated once either |x_i| < 0.25). Beyond that point the
> weakly supplied module has reached its deleterious-mutation load floor and the
> other is at the fitness threshold that ends the simulation, so the terminal log
> ratio is set by the stopping conditions rather than by the dynamics.

**Why this matters, with the numbers.** At the default parameters module 2 settles
at `|x_2*| = U(1-f_2)·sigma^2·a_2^2/delta = 0.072` (simulation gives 0.073) and
module 1 is caught by the `W >= 0.99` rule at `|x_1| ≈ 0.23`. Neither depends on the
initial condition, so with no cutoff all six trajectories are pinned to nearly the
same terminal value and the panel reads as convergence onto an attractor — landing
near the orange line, which makes it worse. The endpoint spread (initial spread
3.41) behaves as:

| cutoff | endpoint spread | generations retained, angle 1 |
|---|---|---|
| 0 | 0.65 | 3341 |
| 0.15 | 1.19 | 579 |
| 0.25 | 2.06 | 363 |
| 0.40 | 3.01 | 78 |

0.25 is the compromise: the trajectories stay clearly separated while each still
retains a usable stretch of its adapting phase.

SSWM needs no cutoff — there are no deleterious mutations, so there is no load
floor. Its curves nonetheless end where `x_2` reaches the last nonzero lattice
site, at `log(x_2/x_1) = log(delta/|x_1|)` exactly. That is lattice discreteness,
not dynamics, and it is why the divergence in panel D looks shallower than the
closed form predicts.

---

## Three things to check on the first run

### 1. `mybinornd` is supplied by the repository

`mybinornd` is called by `simulateModularSSWM.m`,
`simulatePleiotropicSSWM.m`, `simulateNestedSSWM.m` and by
`simulateModularSSWM_ConstantSupply.m` here. It is defined in `utils/`, so a clean
checkout can run the successive-mutations simulations without any file outside the
project. `Run_modularConstantSupply` checks for it and errors with an explicit
message rather than failing inside a `parfor`.

### 2. Does the dashed black line in panel D track the simulated drift?

It is the closed-form SSWM prediction plotted against *absolute* generations, so it
tests the prefactor of `beta_i`, not just the shape. I derived
`beta_i = N*U*f_i*delta^2/(sigma^2*a_i^2)` from `flux = N*U_i*2*s_i*delta` with
`s_i = |x_i|*delta/(sigma^2*a_i^2)`.

Applying the same derivation to the declining-supply model gives
`alpha_i = N*U*delta/(L_i*sigma^2*a_i^2)`, whereas `predictModularSSWM.m` line ~50
has `alpha_i = 4*N*U*delta/(L_i*a_i^2)` — larger by a factor of `sigma^4 = 16` at
`sigma = 2`. This does **not** affect Figure 3: `predictModularSSWM` calibrates
`t_final` from the simulated endpoint, and the trajectory shape in the `(x1, x2)`
plane depends only on `alpha_2/alpha_1`, so the prefactor cancels. But it would
matter if a predicted trajectory were ever plotted against generations. Worth
confirming against the paper's equations rather than taking my word for it — I
derived it independently and may be missing a convention.

### 3. The CM sexual warning

`simulateModularCM_ConstantSupply` prints "Pre-run population data not found for
sexual CM. Using monomorphic initialization." because `Run_modularConstantSupply`
mirrors `Run_modularFGM`, which also does not pre-run for the modular CM sexual
case. If you want standing variation, build `simParams.populationMatrices` with
`freezeParam` + `preRunSimulation` before the `cm_sexual` stage, as elsewhere in
the project.

---

## Floating-point note

Repeated `x = x + delta` on the lattice can leave a residue of order `1e-17`
instead of an exact zero. The existing code gates the supply on
`max(-x_i, 0) > 0`, which treats such a residue as a nonzero distance from the
optimum and can produce enormous waiting times and a wildly negative
`log(x_2/x_1)` at the very end of a run. The files here gate on
`x_i <= -delta/2` instead — exact on the lattice, robust to the residue — and
`makeFigure_ConstantSupply` filters `|x| > 1e-9` before taking the log. If you ever
see a stray spike at the end of a Figure 3 panel, this is the likely cause.
