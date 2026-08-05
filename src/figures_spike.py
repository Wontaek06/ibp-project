"""
Report figures for the stage-1 spike.

    python -m src.figures_spike     # after python -m src.spike

Produces:
  figures/spike_cold_adaptation.png  - habitat temperature, IBP vs controls
  figures/spike_lat_vs_sst.png       - latitude/temperature space by clade
  figures/spike_stage_coverage.png   - per-stage data return rate

Palette is Okabe-Ito (colourblind-safe). Slot order is chosen so the two
hues that sit closest under deuteranopia (green/purple) are never
adjacent, and every point is also directly labelled, so identity never
rests on colour alone.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

RESULTS = "data/spike_results.csv"
FIG_DIR = "figures"

CLADE_COLOR = {
    "Bacteria": "#0072B2",
    "Fungi": "#E69F00",
    "Chromista": "#009E73",
    "Animalia": "#56B4E9",
}
CONTROL_EDGE = "#333333"

plt.rcParams.update({
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "figure.facecolor": "white",
})


def load():
    df = pd.read_csv(RESULTS)
    df["group"] = df.expect_cold.map({"yes": "IBP-bearing", "no": "Control (no AFP)"})
    return df


def _beeswarm(values, span, step=0.075, tol=0.045):
    """
    Horizontal offsets that stop near-equal values from overplotting.

    Deterministic (no RNG) so the figure is byte-identical across runs.
    Values within `tol` of the running cluster centre fan out in a
    0, +1, -1, +2, -2 pattern.
    """
    offsets, cluster = [], []
    threshold = (span or 1) * tol
    for v in values:
        if cluster and abs(v - cluster[-1]) > threshold:
            cluster = []
        rank = len(cluster)
        direction = 1 if rank % 2 else -1
        offsets.append(direction * ((rank + 1) // 2) * step)
        cluster.append(v)
    return offsets


def _legend_handles(plotted):
    """
    Legend entries for the clades actually drawn.

    Built from the plotted rows rather than from CLADE_COLOR, because the
    two terrestrial fungi drop out of any temperature panel (Bio-ORACLE is
    marine-only) - a standing "Fungi" swatch with no matching point reads
    as a missing series rather than an excluded one.
    """
    present = [c for c in CLADE_COLOR
               if ((plotted.expect_cold == "yes") & (plotted.clade == c)).any()]
    handles = [plt.Line2D([], [], marker="o", ls="", markersize=9,
                          markerfacecolor=CLADE_COLOR[c], markeredgecolor="white",
                          label=f"IBP — {c}")
               for c in present]
    if (plotted.expect_cold == "no").any():
        handles.append(plt.Line2D([], [], marker="s", ls="", markersize=9,
                                  markerfacecolor="white", markeredgecolor=CONTROL_EDGE,
                                  label="Control (no AFP)"))
    return handles


def _label(ax, xy, text, dx=11, dy=-3):
    """Annotate with a surface-coloured halo so labels stay legible on marks."""
    ax.annotate(text, xy, textcoords="offset points", xytext=(dx, dy),
                fontsize=8, color="#444444", zorder=5,
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                          edgecolor="none", alpha=0.82))


def _footnote(fig, df):
    """State which seeds are absent from a temperature panel, and why."""
    dropped = df[df.surf_min_temp.isna() & (df.expect_cold == "yes")]
    if dropped.empty:
        return

    # Name the terrestrial values rather than just declaring the taxa
    # missing - they are colder than anything in the marine panel, and
    # omitting the numbers would read as missing data instead of a
    # different measurement.
    if "land_min_temp" in dropped and dropped.land_min_temp.notna().any():
        parts = ", ".join(f"{r.label} {r.land_min_temp:.1f} °C"
                          for _, r in dropped.iterrows() if pd.notna(r.land_min_temp))
        note = (f"Shown separately (2 m air temperature, not comparable to sea water): "
                f"{parts} — terrestrial fungi, no marine grid cell.")
    else:
        note = (f"Excluded: {', '.join(dropped.label)} — terrestrial fungi, "
                f"no Bio-ORACLE marine cell (coordinates resolved normally).")
    fig.text(0.01, 0.012, note, fontsize=7.5, color="#666666")


def _style(row):
    """IBP-bearing taxa are filled by clade; controls are hollow."""
    if row.expect_cold == "yes":
        return dict(facecolor=CLADE_COLOR.get(row.clade, "#999999"),
                    edgecolor="white", marker="o")
    return dict(facecolor="white", edgecolor=CONTROL_EDGE, marker="s")


def fig_cold_adaptation(df):
    d = df.dropna(subset=["surf_min_temp"])
    groups = ["IBP-bearing", "Control (no AFP)"]
    data = [d[d.group == g].surf_min_temp.values for g in groups]

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    bp = ax.boxplot(data, tick_labels=groups, widths=0.45, showfliers=False,
                    medianprops=dict(color="#333333", linewidth=2),
                    boxprops=dict(color="#666666"),
                    whiskerprops=dict(color="#666666"),
                    capprops=dict(color="#666666"))

    # Label only the extremes and the notable outlier. Labelling all 14
    # points collides badly - the IBP group spans barely 6 C - and the
    # per-taxon values are in data/spike_results.csv anyway.
    labelled = set()
    for g in groups:
        sub = d[d.group == g]
        labelled.add(sub.surf_min_temp.idxmin())
        labelled.add(sub.surf_min_temp.idxmax())

    for i, g in enumerate(groups, start=1):
        sub = d[d.group == g].sort_values("surf_min_temp")
        offsets = _beeswarm(sub.surf_min_temp.values, span=d.surf_min_temp.max() - d.surf_min_temp.min())
        for (idx, row), offset in zip(sub.iterrows(), offsets):
            st = _style(row)
            ax.scatter(i + offset, row.surf_min_temp, s=88, zorder=3,
                       linewidths=1.2, **st)
            if idx in labelled:
                _label(ax, (i + offset, row.surf_min_temp), row.label)

    ax.axhline(0, ls="--", lw=0.9, color="#888888", zorder=1)
    ax.annotate("0 °C", (0.52, 0.15), fontsize=8, color="#888888")

    cold = d[d.group == "IBP-bearing"].surf_min_temp
    warm = d[d.group == "Control (no AFP)"].surf_min_temp
    u, p = mannwhitneyu(cold, warm, alternative="less")
    gap = warm.median() - cold.median()

    ax.set_ylabel("Annual minimum sea surface temperature (°C)\nat distribution centre")
    ax.set_title("Ice-binding-protein taxa occupy colder habitats\n"
                 f"median separation {gap:.1f} °C · Mann–Whitney U={u:.0f}, p={p:.4f} "
                 f"(n={len(cold)} vs {len(warm)})", fontsize=11)

    handles = _legend_handles(d)
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left")
    _footnote(fig, df)

    fig.tight_layout(rect=(0, 0.045, 1, 1))
    path = f"{FIG_DIR}/spike_cold_adaptation.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def fig_lat_vs_sst(df):
    d = df.dropna(subset=["surf_min_temp", "rep_lat"])
    fig, ax = plt.subplots(figsize=(7.6, 5.6))

    for _, row in d.iterrows():
        st = _style(row)
        ax.scatter(row.rep_lat, row.surf_min_temp, s=95, zorder=3,
                   linewidths=1.2, **st)
        _label(ax, (row.rep_lat, row.surf_min_temp), row.label, dx=9, dy=4)

    ax.axhline(0, ls="--", lw=0.9, color="#888888", zorder=1)
    ax.set_xlabel("Distribution-centre latitude (|°|)")
    ax.set_ylabel("Annual minimum sea surface temperature (°C)")
    ax.set_title("Habitat temperature against latitude\n"
                 "IBP-bearing taxa cluster in the cold, high-latitude corner",
                 fontsize=11)

    handles = _legend_handles(d)
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper right")

    fig.tight_layout()
    path = f"{FIG_DIR}/spike_lat_vs_sst.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def fig_stage_coverage(df):
    """Per-stage data return rate - the feasibility claim, as a figure."""
    stages = [("ok_uniprot", "UniProt\nsequence"), ("ok_organism", "Source\norganism"),
              ("ok_gbif_taxon", "GBIF\ntaxon"), ("ok_occurrences", "Occurrence\ncoords"),
              ("ok_environment", "Habitat\ntemperature"), ("ok_alphafold", "AlphaFold\nstructure")]

    labels, passed, totals = [], [], []
    for col, label in stages:
        vals = df[col].dropna()
        applicable = vals[vals.isin([True, False])] if vals.dtype == object else vals
        n_ok = int(applicable.sum())
        labels.append(label)
        passed.append(n_ok)
        totals.append(len(applicable))

    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    y = range(len(labels))
    ax.barh(list(y), totals, color="#E8E8E8", height=0.55, zorder=2)
    bars = ax.barh(list(y), passed, color="#0072B2", height=0.55, zorder=3)

    for i, (n_ok, tot) in enumerate(zip(passed, totals)):
        full = n_ok == tot
        ax.annotate(f"{n_ok}/{tot}" + ("" if full else "  ← 2 terrestrial fungi"),
                    (n_ok, i), textcoords="offset points", xytext=(8, 0),
                    va="center", fontsize=9,
                    color="#0072B2" if full else "#D55E00", weight="bold")

    ax.set_yticks(list(y), labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Seeds returning data")
    ax.set_xlim(0, max(totals) * 1.28)
    ax.grid(axis="y", visible=False)
    ax.set_title("Every pipeline stage returns data\n"
                 "stage-1 spike, 16 seeds (13 IBP-bearing + 3 controls)", fontsize=11)
    fig.tight_layout()
    path = f"{FIG_DIR}/spike_stage_coverage.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    df = load()
    for path in (fig_cold_adaptation(df), fig_lat_vs_sst(df), fig_stage_coverage(df)):
        print(f"Saved {path}")


if __name__ == "__main__":
    main()
