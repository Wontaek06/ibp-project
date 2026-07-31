"""
Stage-2 environment mapping: habitat temperature for the DUF3494 family.

    python -m src.expand_env                # resumable; safe to re-run
    python -m src.expand_env --limit 80     # short batch first
    python -m src.expand_env --full         # all three marine layers

Maps at GENUS level, not per taxon. The 1,835 representatives span 1,117
taxa but only 441 genera, and most of the taxa are undescribed strains
("Colwellia sp. SLW05") whose species epithet GBIF has never heard of
anyway - so genus is both the cheaper and the more honest resolution.

By default only `surf_min_temp` is fetched (one layer instead of three),
because that is the variable the analysis actually uses and each extra
layer triples the request count. `--full` adds bottom temperature and sea
ice thickness.

Results append to data/pfam11999_env.csv after every genus, and genera
already present are skipped, so an interrupted run loses nothing.

READ THIS BEFORE USING THE OUTPUT
---------------------------------
Occurrence coordinates are a much weaker habitat proxy for bacteria than
for fish. A GBIF record for a marine bacterium is usually a metabarcoding
survey hit, so it marks where somebody sequenced seawater, not the
organism's thermal niche - and genus-level ranges are enormous. The
stage-1 spike measured this failure directly: Leucosporidium sp. AY30 is
an Arctic strain, but the genus is distributed across Finnish boreal
forest and lands at 61 N. Treat these numbers as a coarse screen, and
prefer strain collection coordinates (NCBI BioSample) for any claim about
a specific protein.
"""
import argparse
import os

import pandas as pd

from src.fetch_bio_oracle import bo_point, bo_species, CONFIG
from src.fetch_land_climate import land_min_temp
from src.occurrence import fetch_occurrences, core_points
from src.taxonomy import gbif_match

FEATURES = "data/pfam11999_features.csv"
OUT = "data/pfam11999_env.csv"

# features.csv clade label -> GBIF backbone kingdom
CLADE_TO_KINGDOM = {
    "Bacteria": "Bacteria", "Archaea": "Archaea", "Fungi": "Fungi",
    "Diatoms": "Chromista", "Metazoa": "Animalia", "Plants": "Plantae",
}

PRIMARY_LAYER = "surf_min_temp"

# Bio-ORACLE search radius for this module, deliberately much tighter than
# the 0.5 deg default. The wide default exists to rescue coastal MARINE
# species whose exact coordinate lands on a masked cell - but applied to a
# soil bacterium it walks inland-to-offshore and invents a sea temperature.
# Measured: Sphaerisporangium, Nocardioides and Nonomuraea (all soil
# actinomycetes) were classified "marine" at 11-15 C. DUF3494 is common in
# terrestrial microbes, so this misclassification would have hit a large
# share of the family. ~0.15 deg is roughly 15 km at these latitudes.
SEARCH_RADIUS_DEG = 0.15

# Habitat is decided on the EXACT grid cell (radius 0), with no neighbour
# search at all - a neighbour search is precisely what smuggles land taxa
# into the ocean. Validated against five genera with known habitat:
#   Colwellia 3/3 sea, Fragilariopsis 2/3 sea  -> marine   (correct)
#   Typhula 1/3, Nocardioides 0/3, Leucosporidium 0/3 -> land (correct)
# Nocardioides is the clean demonstration: null at radius 0 and 0.05, but
# 10.8 C at 0.15 - the value came entirely from an offshore neighbour.
HABITAT_RADIUS_DEG = 0.0

# A genus counts as marine only if at least this fraction of its localities
# resolve to a sea cell; otherwise it is treated as terrestrial.
MARINE_FRACTION = 0.5


def genus_table():
    """One row per genus, carrying its clade and how many proteins it holds."""
    df = pd.read_csv(FEATURES)
    df["genus"] = df.organism.astype(str).str.split().str[0]
    df = df[df.genus.str.match(r"^[A-Z][a-z]+$").fillna(False)]  # drop junk tokens
    g = (df.groupby(["genus", "clade"])
           .agg(n_proteins=("accession", "size"), n_taxa=("taxon_id", "nunique"))
           .reset_index()
           .sort_values("n_proteins", ascending=False))
    return g


def done_genera():
    if not os.path.exists(OUT):
        return set()
    try:
        return set(pd.read_csv(OUT).genus)
    except Exception:
        return set()


def append_row(row):
    header = not os.path.exists(OUT)
    pd.DataFrame([row]).to_csv(OUT, mode="a", header=header, index=False,
                               encoding="utf-8-sig")


def process(genus, clade, k, layers):
    rec = {"genus": genus, "clade": clade}
    kingdom = CLADE_TO_KINGDOM.get(clade)

    m = gbif_match(genus, kingdom=kingdom, rank="GENUS")
    rec["gbif_key"] = m.get("usageKey")
    rec["gbif_match"] = m.get("matchType")
    if not rec["gbif_key"]:
        return rec, "no GBIF taxon"

    taxon = {"usageKey": rec["gbif_key"], "name": genus}
    records, provider = fetch_occurrences(taxon)
    rec["n_occ"] = len(records)
    rec["occ_provider"] = provider
    if not records:
        return rec, "no occurrences"

    pts = core_points(records, k=k)
    abs_lats = sorted(abs(p["lat"]) for p in pts)
    rec["rep_lat"] = round(abs_lats[len(abs_lats) // 2], 2)

    # Decide marine vs terrestrial from how many localities actually sit on
    # a sea cell, rather than from "did any layer return a number" - the
    # latter is satisfied by a single coastal neighbour and mislabels soil
    # taxa as marine.
    dsid, var = CONFIG[PRIMARY_LAYER]
    on_sea = [bo_point(dsid, var, p["lat"], p["lon"], HABITAT_RADIUS_DEG) for p in pts]
    n_sea = sum(1 for v in on_sea if v is not None)
    rec["n_marine_pts"] = n_sea
    rec["n_pts"] = len(pts)

    if n_sea >= max(1, MARINE_FRACTION * len(pts)):
        rec["habitat"] = "marine"
        # Values may now use the wider radius: the taxon is established as
        # marine, so a coastal cell mask is a gap to bridge, not a hint.
        for label in layers:
            d2, v2 = CONFIG[label]
            rec[label] = bo_species(d2, v2, pts, radius_deg=SEARCH_RADIUS_DEG)
    else:
        rec["land_min_temp"] = land_min_temp(pts)
        rec["habitat"] = "terrestrial" if rec["land_min_temp"] is not None else "unresolved"

    return rec, None


def main():
    ap = argparse.ArgumentParser(description="Habitat temperature per DUF3494 genus")
    ap.add_argument("--limit", type=int, help="process at most N new genera")
    ap.add_argument("--full", action="store_true",
                    help="fetch all three marine layers instead of surface minimum only")
    ap.add_argument("-k", type=int, default=3, help="occurrence points per genus (default 3)")
    args = ap.parse_args()

    layers = list(CONFIG) if args.full else [PRIMARY_LAYER]
    table = genus_table()
    already = done_genera()
    todo = table[~table.genus.isin(already)]
    if args.limit:
        todo = todo.head(args.limit)

    print(f"{len(table)} genera total, {len(already)} done, processing {len(todo)}")
    print(f"layers={layers}  points/genus={args.k}  "
          f"~{len(todo) * args.k * len(layers):,} Bio-ORACLE requests\n")

    for i, (_, row) in enumerate(todo.iterrows(), start=1):
        rec, problem = process(row.genus, row.clade, args.k, layers)
        rec["n_proteins"] = row.n_proteins
        append_row(rec)
        status = problem or (f"{rec.get(PRIMARY_LAYER)} C" if rec.get("habitat") == "marine"
                             else f"AIR {rec.get('land_min_temp')} C")
        print(f"[{i:3d}/{len(todo)}] {row.genus:24s} {row.clade:9s} "
              f"lat={rec.get('rep_lat')} {status}", flush=True)

    print(f"\nSaved {OUT}")


if __name__ == "__main__":
    main()
