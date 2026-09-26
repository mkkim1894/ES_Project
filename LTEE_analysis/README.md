# LTEE re-analysis

Diversity of gene targets over time in the six non-mutator populations of
Lenski's Long-Term Evolution Experiment, from the metagenomic time courses of
Good et al. (2017) *Nature* 551:45-50.

Everything is in one script, `ltee_analysis.py`. It replaces the earlier
Jupyter notebook, which is kept in `_relegated/LTEE_notebooks/`.

## Run

    cd LTEE_analysis
    python3 ltee_analysis.py --download      # first run only: fetch the data
    python3 ltee_analysis.py                 # afterwards

Requires `pandas`, `numpy` and `matplotlib`. Paths are resolved relative to the
script, so the working directory does not matter. `--data-dir` points at an
existing copy of the Good et al. `data_files` folder; `--out-dir` changes where
output goes.

`--download` fetches only the 12 files the analysis needs (about 8 MB) from
https://github.com/benjaminhgood/LTEE-metagenomic, rather than the whole
repository. The data folder is git-ignored.

## What it does

1. Reads the annotated time course of each population and keeps mutations with
   `PASS` status, excluding intergenic ones (no unambiguous gene assignment).
   Appearance time is the first sampled generation at which the HMM state of
   Good et al. becomes positive.
2. Builds multi-hit gene sets at multiplicity thresholds m >= 2, 3 and 4.
3. Computes Simpson's concentration index D = sum_i p_i^2 over genes in sliding
   windows 10,000 generations wide, offset by 2,500.
4. Rarefies each window to the smallest qualifying window count, 2,000
   subsamples without replacement, and reports the mean of D with 2.5/97.5
   percentiles. The manuscript plots the reciprocal of the mean D. Averaging D
   and inverting once, rather than averaging 1/D, is the estimator used
   throughout the paper, including for n_sel in the multi-module simulations.
5. Splits the time course into an early epoch (t_a <= 17,500) and a late epoch,
   the shortest terminal interval holding at least as many mutations.

## Outputs, in `results/` (git-ignored)

- `F7.pdf` main-text LTEE figure
- `FigureS_LTEE_multihit.pdf` multi-hit sensitivity
- `ltee_report.txt` every number the manuscript quotes
- `ltee_rarefied_all_genes.csv`, `ltee_rarefied_multihit_m{2,3,4}.csv` per-window
  D, its percentiles and the reciprocals actually plotted
- `ltee_windows_all_genes.csv` raw per-window counts and D
- `ltee_mutations.csv` the mutation table the analysis is built from

## Notes on the numbers

The last appearance time in the data is generation **60,500**, not 60,000. The
late epoch therefore starts at t_a > 35,000 and spans 25,500 generations, which
is why both figures appear in the manuscript.

The `min_count = 20` guard on window qualification never binds: the smallest
qualifying window already holds 104 mutations. It is kept as a guard against
degenerate windows rather than as an active filter.

The final window contains exactly 104 mutations, so every rarefaction subsample
is the whole window and its confidence band collapses to a point. This is
expected, not a plotting error.

The correspondence between the annotated file and the HMM state file is by row
index. The script asserts that the two have equal length rather than assuming
it; the earlier notebook silently dropped rows if they disagreed.
