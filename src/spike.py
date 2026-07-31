"""
Stage-1 spike: drive a hand-curated IBP set through every pipeline stage
once, and report where data stops coming back.

    python -m src.spike

This is a feasibility probe, not an analysis. It answers exactly two
questions before any scale-up work is done:

  1. Does every stage return data, or is there a stage that comes back
     empty for a whole clade? (per-stage coverage table)
  2. Do organisms that are known to be cold-adapted actually land on
     cold coordinates? (cold vs. negative-control temperature contrast)

The seed set is pinned by UniProt accession rather than keyword search -
see uniprot_by_accession() for why. Negative controls carry no accession
and only exercise the occurrence -> environment path.
"""
import pandas as pd

from src.fetch_uniprot import uniprot_by_accession
from src.taxonomy import resolve_taxon
from src.occurrence import fetch_occurrences, core_points, cold_points
from src.fetch_bio_oracle import bo_species, CONFIG
from src.fetch_alphafold import alphafold_meta

SEED_PATH = "data/spike_seed.csv"
OUT_PATH = "data/spike_results.csv"

STAGES = ["uniprot", "organism", "gbif_taxon", "occurrences", "environment", "alphafold"]


def _blank(value):
    """True for None, NaN and empty/whitespace strings.

    pandas reads an empty CSV cell as NaN, which is truthy - so a plain
    `row.get("accession") or None` sends the negative controls off to
    fetch UniProt entry "nan".
    """
    return value is None or pd.isna(value) or not str(value).strip()


def run_one(row):
    """One seed row -> a flat result dict with a per-stage ok flag."""
    acc = row.get("accession")
    rec = {
        "label": row["label"], "clade": row["clade"], "afp_class": row["afp_class"],
        "expect_cold": row["expect_cold"],
        "accession": None if _blank(acc) else str(acc).strip(),
    }
    ok = dict.fromkeys(STAGES, False)

    # --- Stage 1: UniProt -------------------------------------------------
    override = row.get("organism_override")
    organism = None if _blank(override) else str(override).strip()
    lineage = []
    if rec["accession"]:
        up = uniprot_by_accession(rec["accession"])
        if up and up.get("seq"):
            ok["uniprot"] = True
            organism, lineage = up["organism"], up["lineage"]
            rec.update(prot_name=up["prot_name"], organism=organism,
                       taxon_id=up["taxon_id"], seq_len=len(up["seq"]),
                       reviewed=up["reviewed"])
    else:
        ok["uniprot"] = None  # not applicable: negative controls have no protein
        rec.update(prot_name=None, organism=organism, taxon_id=None,
                   seq_len=None, reviewed=None)

    # --- Stage 2: source organism ----------------------------------------
    if organism:
        ok["organism"] = True
        rec["organism"] = organism
    else:
        print(f"  {rec['label']:16s} STOPPED - no source organism")
        return rec, ok

    # --- Stage 3: GBIF taxon resolution ----------------------------------
    taxon = resolve_taxon(organism, lineage)
    rec.update(gbif_name=taxon["name"], gbif_rank=taxon["rank"],
               gbif_key=taxon["usageKey"], gbif_match=taxon["matchType"],
               kingdom=taxon["kingdom"])
    ok["gbif_taxon"] = taxon["usageKey"] is not None

    # --- Stage 4: occurrences --------------------------------------------
    # Habitat temperature is read at the distribution centre; the poleward
    # tail is recorded as latitude only, for the range-limit question.
    records, provider = fetch_occurrences(taxon)
    pts = core_points(records, k=5)
    tail = cold_points(records, k=5)
    rec.update(n_occ=len(records), occ_provider=provider, n_core_pts=len(pts))
    ok["occurrences"] = len(pts) > 0

    if pts:
        lats = sorted(abs(p["lat"]) for p in pts)
        rec["rep_lat"] = round(lats[len(lats) // 2], 2)
    if tail:
        rec["tail_lat"] = round(max(abs(p["lat"]) for p in tail), 2)

    # --- Stage 5: environment (Bio-ORACLE) --------------------------------
    if pts:
        for label, (dsid, var) in CONFIG.items():
            rec[label] = bo_species(dsid, var, pts)
        ok["environment"] = rec.get("surf_min_temp") is not None

    # --- Stage 6: AlphaFold ----------------------------------------------
    if rec["accession"]:
        af = alphafold_meta(rec["accession"])
        rec.update(af=af["af"], plddt=af["plddt"])
        ok["alphafold"] = bool(af["af"])
    else:
        ok["alphafold"] = None  # not applicable

    print(f"  {rec['label']:16s} {str(rec.get('organism'))[:28]:30s} "
          f"gbif={rec.get('gbif_key')} occ={rec.get('n_occ')}({rec.get('occ_provider')}) "
          f"lat={rec.get('rep_lat')} sst_min={rec.get('surf_min_temp')} af={rec.get('af')}")

    for stage in STAGES:
        rec[f"ok_{stage}"] = ok[stage]
    return rec, ok


def coverage_report(oks):
    """Per-stage pass counts, ignoring rows where the stage does not apply."""
    print("\n" + "=" * 62)
    print("STAGE COVERAGE  (n/a = stage does not apply to that row)")
    print("=" * 62)
    for stage in STAGES:
        vals = [o[stage] for o in oks]
        applicable = [v for v in vals if v is not None]
        passed = sum(1 for v in applicable if v)
        na = len(vals) - len(applicable)
        bar = "#" * passed + "." * (len(applicable) - passed)
        flag = "OK " if passed == len(applicable) else "!! "
        print(f"{flag}{stage:14s} {passed:2d}/{len(applicable):2d}  {bar}"
              + (f"   ({na} n/a)" if na else ""))


def cold_contrast(df):
    """The actual hypothesis check: do cold-adapted taxa land on cold cells?"""
    print("\n" + "=" * 62)
    print("COLD-ADAPTATION CHECK  (IBP-bearing vs. negative controls)")
    print("=" * 62)
    have_env = df[df.surf_min_temp.notna()] if "surf_min_temp" in df else df.iloc[0:0]
    if have_env.empty:
        print("no environment data - cannot evaluate")
        return

    for expect, group in have_env.groupby("expect_cold"):
        tag = "IBP-bearing (expect cold)" if expect == "yes" else "controls (expect warm)"
        print(f"\n{tag}: n={len(group)}")
        print(f"  surf_min_temp  median={group.surf_min_temp.median():6.2f} C  "
              f"range=[{group.surf_min_temp.min():.2f}, {group.surf_min_temp.max():.2f}]")
        print(f"  |lat| centre   median={group.rep_lat.median():6.2f} deg  "
              f"(poleward tail median {group.tail_lat.median():.2f} deg)")

    cold = have_env[have_env.expect_cold == "yes"].surf_min_temp
    warm = have_env[have_env.expect_cold == "no"].surf_min_temp
    if len(cold) >= 2 and len(warm) >= 2:
        from scipy.stats import mannwhitneyu
        u, p = mannwhitneyu(cold, warm, alternative="less")
        gap = warm.median() - cold.median()
        print(f"\n  separation = {gap:.2f} C   Mann-Whitney U={u:.1f}, p={p:.4f}")
        print("  -> " + ("SIGNAL PRESENT: IBP taxa map to colder cells"
                         if p < 0.05 else
                         "NO CLEAR SEPARATION - revisit point selection"))


def main():
    seed = pd.read_csv(SEED_PATH).fillna({"organism_override": ""})
    print(f"Spike over {len(seed)} seed entries "
          f"({(seed.expect_cold == 'yes').sum()} IBP-bearing, "
          f"{(seed.expect_cold == 'no').sum()} controls)\n")

    rows, oks = [], []
    for _, row in seed.iterrows():
        rec, ok = run_one(row)
        rows.append(rec)
        oks.append(ok)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_PATH, index=False, encoding="utf-8-sig")
    print(f"\nSaved {OUT_PATH}")

    coverage_report(oks)
    cold_contrast(df)


if __name__ == "__main__":
    main()
