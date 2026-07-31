"""
Stage-2 EDA over the DUF3494 (Pfam PF11999) ice-binding-like family.

    python -m src.fetch_pfam                       # download first
    cd-hit -i data/pfam11999_proteins.fasta \\
           -o data/pfam11999_nr90.fasta -c 0.9 -n 5
    python -m src.expand_eda

Reads the CD-HIT representatives, computes sequence features, and compares
them across clades. Clade is the grouping variable here rather than AFP
type: DUF3494 is one structural family spread across bacteria, fungi and
diatoms, so "which clade" is the question the data can actually answer at
this stage. AFP-type comparisons stay with the curated fish set.

The 13 spike accessions are flagged in the output so the hand-checked
proteins can be located inside the bulk family.
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import kruskal

from src.features import seq_features

TSV = "data/pfam11999_proteins.tsv"
NR_FASTA = "data/pfam11999_nr90.fasta"
SEED = "data/spike_seed.csv"
OUT = "data/pfam11999_features.csv"
FIG_DIR = "figures"

# Lineage token -> clade label. First match wins, so specific clades first.
CLADES = [
    ("Bacillariophyta", "Diatoms"),
    ("Fungi", "Fungi"),
    ("Metazoa", "Metazoa"),
    ("Viridiplantae", "Plants"),
    ("Archaea", "Archaea"),
    ("Bacteria", "Bacteria"),
]

# Clades with too few representatives for a meaningful distribution
MIN_GROUP = 15


def clade_of(lineage):
    for token, label in CLADES:
        if isinstance(lineage, str) and token in lineage:
            return label
    return "Other"


def read_clusters(clstr_path):
    """
    CD-HIT .clstr -> list of (members, cdhit_rep) per cluster.

    The representative is the member marked with a trailing '*', not
    necessarily the first line.
    """
    clusters, members, rep = [], [], None

    def flush():
        if members:
            clusters.append((list(members), rep or members[0]))

    with open(clstr_path) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                flush()
                members, rep = [], None
            else:
                m = re.search(r">(?:sp|tr)\|([^|]+)\|", line)
                if m:
                    members.append(m.group(1))
                    if line.rstrip().endswith("*"):
                        rep = m.group(1)
    flush()
    return clusters


def choose_representatives(clstr_path, preferred):
    """
    One accession per CD-HIT cluster, preferring curated entries.

    CD-HIT picks the longest sequence as cluster representative and knows
    nothing about annotation quality, so characterized reference proteins
    get silently replaced by anonymous TrEMBL entries: ColAFP (A5XB26,
    Swiss-Prot, the Colwellia reference IBP) lost its cluster to
    A0A5C6QDJ8. Since the whole point of the curated seed set is to anchor
    the bulk family, `preferred` accessions win their cluster.
    """
    reps, rescued = set(), []
    for members, cdhit_rep in read_clusters(clstr_path):
        override = next((a for a in members if a in preferred), None)
        if override and override != cdhit_rep:
            reps.add(override)
            rescued.append(override)
        else:
            reps.add(cdhit_rep)
    return reps, rescued


def read_representatives(path, preferred=frozenset()):
    """
    CD-HIT output -> set of representative accessions.

    Uses the .clstr file when available so representatives can be
    re-chosen; falls back to the FASTA headers otherwise.
    """
    clstr = f"{path}.clstr"
    if os.path.exists(clstr):
        reps, rescued = choose_representatives(clstr, preferred)
        if rescued:
            print(f"Re-chose {len(rescued)} cluster representative(s) to keep "
                  f"curated entries: {', '.join(sorted(rescued))}")
        return reps

    if not os.path.exists(path):
        return None
    accs = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                m = re.match(r">(?:sp|tr)\|([^|]+)\|", line)
                accs.add(m.group(1) if m else line[1:].split()[0])
    return accs


def build_table():
    df = pd.read_csv(TSV, sep="\t")
    df = df.rename(columns={
        "Entry": "accession", "Reviewed": "reviewed", "Protein names": "prot_name",
        "Organism": "organism", "Organism (ID)": "taxon_id",
        "Taxonomic lineage": "lineage", "Length": "length", "Sequence": "seq",
    })

    seed = pd.read_csv(SEED)
    spike_accs = set(seed.accession.dropna())
    curated = spike_accs | set(df[df.reviewed == "reviewed"].accession)

    reps = read_representatives(NR_FASTA, preferred=curated)
    if reps is None:
        print(f"[warn] {NR_FASTA} missing - running on all {len(df):,} entries "
              f"(near-duplicates NOT removed)")
    else:
        before = len(df)
        df = df[df.accession.isin(reps)].copy()
        print(f"CD-HIT 90%: {before:,} -> {len(df):,} representatives")

    df["clade"] = df.lineage.map(clade_of)

    feats = pd.DataFrame([seq_features(s) for s in df.seq], index=df.index)
    df = pd.concat([df.drop(columns=["length"]), feats], axis=1)

    df["in_spike"] = df.accession.isin(spike_accs)
    return df


def report(df):
    print(f"\nRepresentatives: {len(df):,}   distinct taxa: {df.taxon_id.nunique():,}")
    print(f"Swiss-Prot reviewed: {(df.reviewed == 'reviewed').sum():,}")
    print(f"Spike accessions present: {df.in_spike.sum()} / 13")
    if df.in_spike.any():
        print("  " + ", ".join(sorted(df[df.in_spike].accession)))

    print("\nClade composition:")
    for clade, n in df.clade.value_counts().items():
        print(f"  {clade:12s} {n:5,d}  ({100 * n / len(df):.1f}%)")

    cols = ["length", "ala_pct", "thr_pct", "gravy", "pI", "instability"]
    big = df[df.clade.isin(df.clade.value_counts()[lambda s: s >= MIN_GROUP].index)]
    groups = sorted(big.clade.unique())

    print(f"\nFeature medians by clade (groups with n >= {MIN_GROUP}):")
    print(big.groupby("clade")[cols].median().round(2).to_string())

    print("\nKruskal-Wallis across clades:")
    for col in cols:
        samples = [big[big.clade == g][col].dropna().values for g in groups]
        samples = [s for s in samples if len(s) >= 3]
        if len(samples) >= 2:
            h, p = kruskal(*samples)
            mark = "*" if p < 0.05 else " "
            print(f"  {col:12s} H={h:8.2f}  p={p:.3e} {mark}")

    print("""
NOTE on `length`: this is whole-protein length, not domain length. The
DUF3494 domain is ~230 aa, which is roughly what the diatom, fungal and
metazoan entries measure - those proteins are essentially the domain plus
a signal peptide. The bacterial (441) and archaeal (861) medians are
domain-fusion architectures, not larger ice-binding domains; MpAFP's RTX
repeats are the familiar example. So the length result is a statement
about protein architecture per clade, and length should not be fed to a
model as if it described the ice-binding module itself.""")
    return big, groups, cols


def make_figures(big, groups, cols):
    os.makedirs(FIG_DIR, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for ax, col in zip(axes.ravel(), cols):
        ax.boxplot([big[big.clade == g][col].dropna() for g in groups],
                   tick_labels=groups, showfliers=False)
        ax.set_title(col)
        ax.tick_params(axis="x", rotation=20)
    fig.suptitle("DUF3494 (PF11999) sequence features by clade "
                 f"- CD-HIT 90% representatives, n={len(big):,}")
    fig.tight_layout()
    path = f"{FIG_DIR}/duf3494_features_by_clade.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"\nSaved {path}")


def main():
    df = build_table()
    df.drop(columns=["seq"]).to_csv(OUT, index=False, encoding="utf-8-sig")
    print(f"Saved {OUT}")
    big, groups, cols = report(df)
    make_figures(big, groups, cols)


if __name__ == "__main__":
    main()
