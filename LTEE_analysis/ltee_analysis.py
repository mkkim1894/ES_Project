#!/usr/bin/env python3
"""LTEE re-analysis: diversity of gene targets over time.

Re-analysis of the six non-mutator populations of Lenski's Long-Term Evolution
Experiment, using the metagenomic time courses of Good et al. (2017) Nature
551:45-50 (https://github.com/benjaminhgood/LTEE-metagenomic).

Produces the main-text LTEE figure, the multi-hit sensitivity supplementary
figure, and a text/CSV report of every number quoted in the manuscript.

Usage
-----
    python ltee_analysis.py --download      # first run: fetch the data files
    python ltee_analysis.py                 # subsequent runs

    python ltee_analysis.py --data-dir /path/to/LTEE-metagenomic-master/data_files
    python ltee_analysis.py --out-dir results

All paths default to locations relative to this file, so the script can be run
from any working directory.
"""

from __future__ import annotations

import argparse
import sys
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")           # no display needed; figures are written to disk
import matplotlib.pyplot as plt  # noqa: E402

# Keep every glyph inside the core font: mathtext for the comparison
# operators rather than Unicode, which not every backend/font supplies.
plt.rcParams["mathtext.fontset"] = "dejavusans"
plt.rcParams["axes.unicode_minus"] = False

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent

NONMUTATOR_POPULATIONS = ["m5", "m6", "p1", "p2", "p4", "p5"]

T_STAR = 17_500          # early/late boundary, from Good et al. (2017)
WINDOW_SIZE = 10_000     # sliding window width, generations
STEP_SIZE = 2_500        # sliding window offset, generations
MIN_WINDOW_COUNT = 20    # a window must hold this many mutations to qualify
N_BOOTSTRAP = 2_000      # rarefaction subsamples per window
LATE_SEARCH_STEP = 500   # resolution of the late-epoch boundary search
SEED = 42

MULTIPLICITY_THRESHOLDS = (2, 3, 4)

GENOME_LENGTH = 4.7e6    # E. coli REL606, for the genomic-density histogram
N_POSITION_BINS = 80

C_EARLY = "#0077BB"
C_LATE = "#EE7733"
MULTIHIT_COLORS = {2: "#000000", 3: "#800020", 4: "#003153"}

DATA_BASE_URL = (
    "https://raw.githubusercontent.com/benjaminhgood/LTEE-metagenomic/master/data_files"
)
REQUIRED_SUFFIXES = ("_annotated_timecourse.txt", "_well_mixed_state_timecourse.txt")


# --------------------------------------------------------------------------
# Data acquisition
# --------------------------------------------------------------------------

def required_files() -> list[str]:
    return [f"{pop}{suffix}"
            for pop in NONMUTATOR_POPULATIONS
            for suffix in REQUIRED_SUFFIXES]


def download_data(data_dir: Path) -> None:
    """Fetch only the 12 files this analysis needs, not the whole repository."""
    data_dir.mkdir(parents=True, exist_ok=True)
    for name in required_files():
        target = data_dir / name
        if target.exists():
            print(f"  have    {name}")
            continue
        url = f"{DATA_BASE_URL}/{name}"
        print(f"  fetching {name} ...", end=" ", flush=True)
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                target.write_bytes(response.read())
        except urllib.error.URLError as exc:
            target.unlink(missing_ok=True)
            raise SystemExit(
                f"\nCould not download {url}\n  {exc}\n"
                "Download the repository manually from "
                "https://github.com/benjaminhgood/LTEE-metagenomic and point "
                "--data-dir at its data_files folder."
            ) from exc
        print(f"{target.stat().st_size / 1e6:.1f} MB")


def check_data(data_dir: Path) -> None:
    missing = [n for n in required_files() if not (data_dir / n).exists()]
    if missing:
        raise SystemExit(
            f"Missing {len(missing)} of {len(required_files())} data files in "
            f"{data_dir}\n  e.g. {missing[0]}\n\n"
            "Run this script once with --download, or download the repository "
            "from https://github.com/benjaminhgood/LTEE-metagenomic and pass "
            "--data-dir <repo>/data_files."
        )


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

def load_hmm_states(data_dir: Path, population: str):
    """Return (timepoints, states) from a well-mixed HMM state time course."""
    path = data_dir / f"{population}_well_mixed_state_timecourse.txt"
    lines = path.read_text().splitlines()
    timepoints = [int(float(x)) for x in lines[0].strip().split(", ")]
    states = []
    for line in lines[5:]:                      # first five lines are headers
        values = [int(float(x)) for x in line.strip().rstrip(",").split(", ") if x]
        states.append(values)
    return timepoints, states


def appearance_time(state_trajectory, timepoints):
    """First sampled generation at which the HMM state becomes positive."""
    for state, t in zip(state_trajectory, timepoints):
        if state > 0:
            return t
    return None


def build_mutation_table(data_dir: Path, populations=NONMUTATOR_POPULATIONS) -> pd.DataFrame:
    """One row per PASS mutation, with its gene, position and appearance time.

    The well-mixed state file has one row per PASS mutation, in the same order
    as the PASS rows of the annotated file, so row i of one indexes row i of
    the other.  That correspondence is asserted rather than assumed.
    """
    records = []
    for pop in populations:
        annotated = pd.read_csv(data_dir / f"{pop}_annotated_timecourse.txt", header=0)
        annotated.columns = annotated.columns.str.strip()
        annotated["Passed?"] = annotated["Passed?"].str.strip()
        passed = annotated[annotated["Passed?"] == "PASS"].reset_index(drop=True)

        timepoints, states = load_hmm_states(data_dir, pop)
        if len(passed) != len(states):
            raise SystemExit(
                f"{pop}: {len(passed)} PASS mutations but {len(states)} HMM state "
                "rows. The two files are out of register; the row-index "
                "correspondence this analysis relies on does not hold."
            )

        for i, row in passed.iterrows():
            gene = row["Gene"].strip() if isinstance(row["Gene"], str) else row["Gene"]
            records.append({
                "population": pop,
                "position": row["Position"],
                "gene": gene,
                "appearance_time": appearance_time(states[i], timepoints),
            })
    return pd.DataFrame(records)


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------

def simpson_d(labels) -> float:
    """Simpson's concentration index D = sum_i p_i^2."""
    proportions = pd.Series(labels).value_counts(normalize=True)
    return float((proportions ** 2).sum())


def window_bounds(max_time: int):
    for t_start in np.arange(0, max_time - WINDOW_SIZE + STEP_SIZE, STEP_SIZE):
        yield t_start, t_start + WINDOW_SIZE, t_start + WINDOW_SIZE / 2


def compute_windows(data: pd.DataFrame, min_count: int = MIN_WINDOW_COUNT) -> pd.DataFrame:
    """Raw (unrarefied) Simpson's D in sliding windows over appearance time."""
    max_time = int(data["appearance_time"].max())
    rows = []
    for t_start, t_end, t_center in window_bounds(max_time):
        w = data[(data["appearance_time"] > t_start) & (data["appearance_time"] <= t_end)]
        if len(w) >= min_count:
            rows.append({"t_center": t_center, "n": len(w), "D_gene": simpson_d(w["gene"])})
    return pd.DataFrame(rows)


def rarefied_windows(data: pd.DataFrame, n_subsample: int, label: str = "") -> pd.DataFrame:
    """Rarefied Simpson's D per window.

    Each qualifying window is subsampled without replacement to n_subsample
    mutations, N_BOOTSTRAP times.  We report the MEAN of D across subsamples
    (the manuscript then plots its reciprocal) and the 2.5/97.5 percentiles of
    D.  Averaging D and inverting once -- rather than averaging 1/D -- is the
    estimator used throughout, including in the n_sel analysis of the
    multi-module simulations.
    """
    print(f"  rarefaction [{label}]: n={n_subsample} per window, "
          f"{N_BOOTSTRAP} subsamples")
    rng = np.random.default_rng(SEED)
    max_time = int(data["appearance_time"].max())
    rows = []
    for t_start, t_end, t_center in window_bounds(max_time):
        w = data[(data["appearance_time"] > t_start) & (data["appearance_time"] <= t_end)]
        if len(w) < n_subsample:
            continue
        genes = w["gene"].to_numpy()
        d_samples = np.empty(N_BOOTSTRAP)
        for b in range(N_BOOTSTRAP):
            sub = rng.choice(genes, size=n_subsample, replace=False)
            d_samples[b] = simpson_d(sub)
        rows.append({
            "t_center": t_center,
            "n_original": len(w),
            "D_mean": float(d_samples.mean()),
            "D_lo": float(np.percentile(d_samples, 2.5)),
            "D_hi": float(np.percentile(d_samples, 97.5)),
            "inv_D_mean": 1.0 / float(d_samples.mean()),
            "inv_D_lo": 1.0 / float(np.percentile(d_samples, 97.5)),
            "inv_D_hi": 1.0 / float(np.percentile(d_samples, 2.5)),
        })
    return pd.DataFrame(rows)


def late_epoch_start(data: pd.DataFrame, n_early: int) -> int:
    """Shortest terminal interval holding at least as many mutations as the early epoch."""
    max_time = int(data["appearance_time"].max())
    for t_start in range(max_time, 0, -LATE_SEARCH_STEP):
        if len(data[data["appearance_time"] > t_start]) >= n_early:
            return t_start
    raise SystemExit("No terminal interval contains as many mutations as the early epoch.")


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

def make_main_figure(rare_all, n_sub_all, early, late, t_late_start, out_path: Path):
    fig = plt.figure(figsize=(12, 4))
    gs = fig.add_gridspec(1, 3, wspace=0.35)
    ax_a = fig.add_subplot(gs[0])
    gs_right = gs[1:].subgridspec(2, 1, hspace=0)
    ax_early = fig.add_subplot(gs_right[0])
    ax_late = fig.add_subplot(gs_right[1], sharex=ax_early, sharey=ax_early)

    t = rare_all["t_center"] / 1000
    ax_a.fill_between(t, rare_all["inv_D_lo"], rare_all["inv_D_hi"],
                      color="gray", alpha=0.3)
    ax_a.plot(t, rare_all["inv_D_mean"], "o-", color="gray", markersize=5, linewidth=1.5)
    ax_a.axvline(T_STAR / 1000, color="black", linestyle="--", alpha=0.6)
    ax_a.axvline(t_late_start / 1000, color="black", linestyle=":", alpha=0.6)
    ax_a.text(T_STAR / 1000 + 0.5, 0.02, f"t*={T_STAR // 1000}k",
              transform=ax_a.get_xaxis_transform(), fontsize=8,
              va="bottom", ha="left", alpha=0.7)
    ax_a.text(t_late_start / 1000 + 0.5, 0.02, f"late>{t_late_start // 1000}k",
              transform=ax_a.get_xaxis_transform(), fontsize=8,
              va="bottom", ha="left", alpha=0.7)
    ax_a.set_xlabel("Time (thousands of generations)", fontsize=11)
    ax_a.set_ylabel("Effective number of targets (1/D)", fontsize=11)
    ax_a.set_title(f"All genic mutations\n(rarefied, n={n_sub_all}/window)", fontsize=10)
    ax_a.text(-0.2, 1.1, "A", transform=ax_a.transAxes,
              fontsize=14, fontweight="bold", va="top")

    bins = np.linspace(0, GENOME_LENGTH, N_POSITION_BINS + 1)
    centers = (bins[:-1] + bins[1:]) / 2
    width = bins[1] - bins[0]
    scale = 1e7
    for ax, subset, color, name, boundary in (
        (ax_early, early, C_EARLY, "Early", f"$t_a \\leq$ {T_STAR:,} gen"),
        (ax_late, late, C_LATE, "Late", f"$t_a >$ {t_late_start:,} gen"),
    ):
        density, _ = np.histogram(subset["position"], bins=bins, density=True)
        ax.bar(centers, density * scale, width=width, color=color,
               alpha=0.7, edgecolor="white", linewidth=0.5)
        ax.set_ylabel(r"Density ($\times10^{-7}$)", fontsize=10)
        ax.text(0.02, 0.88, f"{name} ({boundary}, n={len(subset)})",
                transform=ax.transAxes, fontsize=9, color=color, fontweight="bold")
    ax_early.set_xlim(bins[0], bins[-1])
    ax_early.set_title("Genomic distribution of mutations", fontsize=11)
    ax_early.text(-0.12, 1.1, "B", transform=ax_early.transAxes,
                  fontsize=14, fontweight="bold", va="top")
    plt.setp(ax_early.get_xticklabels(), visible=False)
    ax_late.set_xlabel("Genomic position (Mb)", fontsize=11)
    ax_late.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x / 1e6:.1f}"))

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_multihit_figure(rare_mh, window_dfs_mh, out_path: Path):
    fig, (ax, ax_n) = plt.subplots(1, 2, figsize=(12, 4))
    for m in MULTIPLICITY_THRESHOLDS:
        rd = rare_mh[m]
        t = rd["t_center"] / 1000
        ax.fill_between(t, rd["inv_D_lo"], rd["inv_D_hi"],
                        color=MULTIHIT_COLORS[m], alpha=0.15)
        ax.plot(t, rd["inv_D_mean"], "o-", color=MULTIHIT_COLORS[m],
                markersize=4, linewidth=1.5, label=f"$m \\geq {m}$")
        wdf = window_dfs_mh[m]
        ax_n.plot(wdf["t_center"] / 1000, wdf["n"], "o-",
                  color=MULTIHIT_COLORS[m], markersize=4, linewidth=1.5,
                  label=f"$m \\geq {m}$")

    for a, ylabel, title, letter in (
        (ax, "Effective number of targets (1/D)", "Multi-hit gene sensitivity", "A"),
        (ax_n, "Mutations per window", "Mutations per window (multi-hit genes)", "B"),
    ):
        a.axvline(T_STAR / 1000, color="black", linestyle="--", alpha=0.6)
        a.set_xlabel("Time (thousands of generations)", fontsize=11)
        a.set_ylabel(ylabel, fontsize=11)
        a.set_title(title, fontsize=11)
        a.text(-0.12, 1.05, letter, transform=a.transAxes,
               fontsize=14, fontweight="bold", va="top")
    ax.legend(fontsize=9, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path,
                        default=HERE / "LTEE-metagenomic-master" / "data_files",
                        help="folder holding the Good et al. (2017) time-course files")
    parser.add_argument("--out-dir", type=Path, default=HERE / "results",
                        help="where figures and the report are written")
    parser.add_argument("--download", action="store_true",
                        help="fetch the 12 required data files before running")
    args = parser.parse_args(argv)

    if args.download:
        print(f"Downloading LTEE data into {args.data_dir}")
        download_data(args.data_dir)
    check_data(args.data_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    report: list[str] = []

    def say(line: str = "") -> None:
        print(line)
        report.append(line)

    say("=" * 70)
    say("LTEE re-analysis")
    say(f"data : {args.data_dir}")
    say(f"out  : {args.out_dir}")
    say("=" * 70)

    # -- load -------------------------------------------------------------
    master = build_mutation_table(args.data_dir)
    df = master[master["gene"] != "intergenic"].copy()
    say()
    say("1. Mutations")
    say(f"   PASS mutations (all)            : {len(master)}")
    say(f"   intergenic, excluded            : {len(master) - len(df)}")
    say(f"   genic mutations analysed        : {len(df)}")
    say(f"   populations                     : {', '.join(NONMUTATOR_POPULATIONS)}")
    say(f"   appearance times span           : {int(df['appearance_time'].min())}"
        f" to {int(df['appearance_time'].max())} generations")

    # -- multi-hit gene sets ---------------------------------------------
    gene_counts = df["gene"].value_counts()
    multihit = {}
    say()
    say("2. Multi-hit gene sets")
    for m in MULTIPLICITY_THRESHOLDS:
        genes = set(gene_counts[gene_counts >= m].index)
        sub = df[df["gene"].isin(genes)].copy()
        multihit[m] = {"genes": genes, "df": sub}
        say(f"   m >= {m}: {len(genes):4d} genes, {len(sub):5d} mutations "
            f"({len(sub) / len(df):.1%} of genic)")

    # -- sliding windows --------------------------------------------------
    window_df = compute_windows(df)
    window_dfs_mh = {m: compute_windows(multihit[m]["df"]) for m in MULTIPLICITY_THRESHOLDS}
    say()
    say("3. Sliding windows "
        f"({WINDOW_SIZE} generations wide, {STEP_SIZE} apart, "
        f"min {MIN_WINDOW_COUNT} mutations)")
    say(f"   all genes: {len(window_df)} qualifying windows, "
        f"n from {int(window_df['n'].min())} to {int(window_df['n'].max())}")
    for m in MULTIPLICITY_THRESHOLDS:
        wdf = window_dfs_mh[m]
        say(f"   m >= {m}  : {len(wdf)} windows, "
            f"n from {int(wdf['n'].min())} to {int(wdf['n'].max())}")

    # -- rarefaction ------------------------------------------------------
    say()
    say("4. Rarefaction")
    n_sub_all = int(window_df["n"].min())
    rare_all = rarefied_windows(df, n_sub_all, label="all genes")
    say(f"   all genes: subsample n={n_sub_all}, {len(rare_all)} windows retained")
    rare_mh, n_sub_mh = {}, {}
    for m in MULTIPLICITY_THRESHOLDS:
        n_sub_mh[m] = int(window_dfs_mh[m]["n"].min())
        rare_mh[m] = rarefied_windows(multihit[m]["df"], n_sub_mh[m], label=f"m>={m}")
        say(f"   m >= {m}  : subsample n={n_sub_mh[m]}, {len(rare_mh[m])} windows retained")

    # -- the numbers the manuscript quotes --------------------------------
    first = rare_all.iloc[0]
    at_star = rare_all.iloc[(rare_all["t_center"] - T_STAR).abs().argmin()]
    last = rare_all.iloc[-1]
    plateau = rare_all[rare_all["t_center"] >= T_STAR]
    say()
    say("5. Effective number of gene targets, 1/D")
    say(f"   first window  (t_center={first['t_center']:.0f}) : {first['inv_D_mean']:.1f}"
        f"  [{first['inv_D_lo']:.1f}, {first['inv_D_hi']:.1f}]")
    say(f"   nearest t*    (t_center={at_star['t_center']:.0f}) : {at_star['inv_D_mean']:.1f}"
        f"  [{at_star['inv_D_lo']:.1f}, {at_star['inv_D_hi']:.1f}]")
    say(f"   last window   (t_center={last['t_center']:.0f}) : {last['inv_D_mean']:.1f}"
        f"  [{last['inv_D_lo']:.1f}, {last['inv_D_hi']:.1f}]")
    say(f"   mean over t_center >= t*                  : {plateau['inv_D_mean'].mean():.1f}"
        f" (range {plateau['inv_D_mean'].min():.1f}-{plateau['inv_D_mean'].max():.1f})")

    # -- early / late epochs ---------------------------------------------
    early = df[df["appearance_time"] <= T_STAR].copy()
    t_late_start = late_epoch_start(df, len(early))
    late = df[df["appearance_time"] > t_late_start].copy()
    max_time = int(df["appearance_time"].max())
    say()
    say("6. Early / late epochs")
    say(f"   early : t_a <= {T_STAR}, n = {len(early)}")
    say(f"   late  : t_a >  {t_late_start}, n = {len(late)}")
    say(f"   last sampled appearance time : {max_time}")
    say(f"   late window length           : {max_time - t_late_start} generations")
    say("   -> manuscript Methods should say 't_a > "
        f"{t_late_start}' and the caption 'final {max_time - t_late_start} generations'")

    # -- outputs ----------------------------------------------------------
    make_main_figure(rare_all, n_sub_all, early, late, t_late_start,
                     args.out_dir / "F7.pdf")
    make_multihit_figure(rare_mh, window_dfs_mh, args.out_dir / "FigureS_LTEE_multihit.pdf")

    rare_all.to_csv(args.out_dir / "ltee_rarefied_all_genes.csv", index=False)
    for m in MULTIPLICITY_THRESHOLDS:
        rare_mh[m].to_csv(args.out_dir / f"ltee_rarefied_multihit_m{m}.csv", index=False)
    window_df.to_csv(args.out_dir / "ltee_windows_all_genes.csv", index=False)
    df.to_csv(args.out_dir / "ltee_mutations.csv", index=False)

    say()
    say("7. Files written")
    for name in ("F7.pdf", "FigureS_LTEE_multihit.pdf",
                 "ltee_rarefied_all_genes.csv", "ltee_windows_all_genes.csv",
                 "ltee_mutations.csv"):
        say(f"   {args.out_dir / name}")

    (args.out_dir / "ltee_report.txt").write_text("\n".join(report) + "\n")
    print(f"\nReport written to {args.out_dir / 'ltee_report.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
