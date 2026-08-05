"""
PAMC (KOPRI Polar and Alpine Microbial Collection) cross-check.

    python -m src.kpdc_crosscheck              # resumable; safe to re-run
    python -m src.kpdc_crosscheck --limit 40   # short batch first

WHY THIS EXISTS
---------------
The stage-2 habitat numbers rest on GBIF occurrence records mapped at
genus level, and stage 1 showed directly how weak that proxy is for
microbes: Leucosporidium sp. AY30 is an Arctic strain, but the genus is
distributed across Finnish boreal forest and its GBIF centre lands at
61 N in a forest, not on ice.

PAMC is the independent check. It is a culture collection of ~6,500
strains isolated from Arctic, Antarctic and alpine sites, and each strain
record carries the coordinate where that strain was actually collected -
strain-level ground truth rather than a genus-level range average. If
DUF3494-bearing genera really are polar organisms, they should appear in
PAMC with polar collection coordinates.

WHAT IT PRODUCES
----------------
  data/pamc_strain_table.csv    one row per PAMC strain matching a DUF3494 genus
  data/pamc_coordinates.csv     strains with usable lat/lon
  data/pamc_sites.csv           distinct collection sites
  data/pamc_genus_latitude.csv  per-genus latitude summary, GBIF vs PAMC

The counts these produce are whatever the data gives. They are not
tuned to match any previously reported figure - the point of a
cross-check is defeated if the target is known in advance.

SCRAPING NOTE
-------------
PAMC publishes no API, so this reads the public search and strain pages.
Requests are serialised with a delay, responses are cached on disk, and
completed genera are skipped on re-run, so a repeat run costs nothing.
Keep REQUEST_DELAY where it is - this is a small research institute's
server, not a CDN.
"""
import argparse
import html
import os
import re
import time

import pandas as pd
import requests

from src.cache import cached

BASE = "https://pamc.kopri.re.kr"
SEARCH = f"{BASE}/k_search"
STRAIN = f"{BASE}/strain/{{strain_id}}"

FEATURES = "data/pfam11999_features.csv"
GBIF_ENV = "data/pfam11999_env.csv"

OUT_STRAINS = "data/pamc_strain_table.csv"
OUT_COORDS = "data/pamc_coordinates.csv"
OUT_SITES = "data/pamc_sites.csv"
OUT_GENUS = "data/pamc_genus_latitude.csv"

REQUEST_DELAY = 1.0  # seconds between requests - be polite
HEADERS = {"User-Agent": "ibp-project academic cross-check (contact via GitHub)"}

# Locality strings are free text; these decide polar vs non-polar.
ARCTIC_TOKENS = ["arctic", "svalbard", "ny-ålesund", "ny-alesund", "greenland",
                 "chukchi", "beaufort", "barents", "kara sea", "alaska",
                 "spitsbergen", "franz josef", "bering"]
ANTARCTIC_TOKENS = ["antarctic", "king george island", "south shetland",
                    "ross sea", "weddell", "terra nova", "mcmurdo",
                    "victoria land", "livingston island", "barton peninsula"]
ALPINE_TOKENS = ["alps", "himalaya", "tibet", "andes", "tauern", "glacier",
                 "alpine", "mount ", "mt. "]


def _strip(fragment):
    return html.unescape(re.sub(r"<[^>]+>", "", fragment)).strip()


def _get(url, params=None, attempts=3, timeout=45):
    """GET with retries. The PAMC host intermittently times out under load."""
    for attempt in range(attempts):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            if attempt == attempts - 1:
                print(f"    [pamc] gave up on {url} {params or ''}: {e}")
                return None
            time.sleep(REQUEST_DELAY * (2 ** attempt))
    return None


def _parse_rows(text, genus):
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", text, re.S):
        cells = [_strip(c) for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)]
        if len(cells) < 5:
            continue
        strain_id = cells[1].replace("*", "").strip()
        species = cells[2]
        if not strain_id:
            continue
        # The search is a free-text match, so it also returns rows where the
        # genus appears in some other field. Keep only real genus matches.
        if not species.lower().startswith(genus.lower()):
            continue
        rows.append({"strain_id": strain_id, "species": species,
                     "locality": cells[3], "habitat": cells[4]})
    return rows


@cached("pamc_search")
def search_genus(genus, max_pages=20):
    """
    Genus name -> matching PAMC strain rows, following pagination.

    The result table shows 20 rows per page, so a genus with more strains
    than that (Flavobacterium, Pseudomonas) is silently truncated unless
    the pages are walked. Stops as soon as a page adds nothing new.
    """
    rows, seen = [], set()
    for page in range(max_pages):
        r = _get(SEARCH, params={"keyword": genus, "page": page})
        if r is None:
            break
        page_rows = _parse_rows(r.text, genus)
        fresh = [x for x in page_rows if x["strain_id"] not in seen]
        if not fresh:
            break
        seen.update(x["strain_id"] for x in fresh)
        rows.extend(fresh)
        if len(page_rows) < 20:
            break
        time.sleep(REQUEST_DELAY)
    return rows


@cached("pamc_strain")
def strain_detail(strain_id):
    """Strain ID -> {latitude, longitude, elevation, sampling_site} or {}."""
    r = _get(STRAIN.format(strain_id=requests.utils.quote(strain_id)))
    if r is None:
        return {}

    text = re.sub(r"<[^>]+>", " ", r.text)
    text = html.unescape(re.sub(r"\s+", " ", text))

    def grab(label, pattern):
        m = re.search(label + r"\s*:?\s*" + pattern, text, re.I)
        return m.group(1) if m else None

    def coord(label):
        """Coordinates are printed as '(N) 78.55' / '(W) 58.47' - sign the value."""
        m = re.search(label + r"\s*:?\s*\(([NSEW])\)\s*([0-9.]+)", text, re.I)
        if not m:
            return None
        hemi, value = m.group(1).upper(), float(m.group(2))
        return -value if hemi in ("S", "W") else value

    return {
        "latitude": coord(r"Latitude\s*\(start\)"),
        "longitude": coord(r"Longitude\s*\(start\)"),
        "elevation": grab(r"Elevation", r"([0-9.\-]+)"),
        "sampling_site": grab(r"Sampling site", r"([^:]{0,120})"),
    }


def classify_region(locality, latitude=None):
    """Locality text (and latitude when present) -> Arctic / Antarctic / Alpine / Other."""
    text = (locality or "").lower()
    if any(t in text for t in ANTARCTIC_TOKENS):
        return "Antarctic"
    if any(t in text for t in ARCTIC_TOKENS):
        return "Arctic"
    if any(t in text for t in ALPINE_TOKENS):
        return "Alpine"
    if latitude is not None:
        if latitude <= -60:
            return "Antarctic"
        if latitude >= 66.5:
            return "Arctic"
    return "Other"


def is_marine(habitat, locality):
    """Whether a strain record describes a marine collection."""
    text = f"{habitat} {locality}".lower()
    marine = ["marine", "sea", "ocean", "seawater", "sediment", "brine",
              "sea ice", "seaice", "fjord", "bay"]
    land = ["soil", "rock", "lichen", "moss", "cryoconite", "lake", "freshwater",
            "glacier ice", "permafrost", "plant", "feces", "guano"]
    if any(t in text for t in land) and not any(t in text for t in ["marine", "ocean", "sea ice"]):
        return False
    return any(t in text for t in marine)


def duf3494_genera():
    """DUF3494-bearing genera, most protein-rich first."""
    df = pd.read_csv(FEATURES)
    df["genus"] = df.organism.astype(str).str.split().str[0]
    df = df[df.genus.str.match(r"^[A-Z][a-z]+$").fillna(False)]
    return (df.groupby(["genus", "clade"]).size()
              .reset_index(name="n_proteins")
              .sort_values("n_proteins", ascending=False))


def done_genera():
    if not os.path.exists(OUT_STRAINS):
        return set()
    try:
        return set(pd.read_csv(OUT_STRAINS).query_genus)
    except Exception:
        return set()


def main():
    ap = argparse.ArgumentParser(description="Cross-check DUF3494 genera against PAMC")
    ap.add_argument("--limit", type=int, help="process at most N new genera")
    args = ap.parse_args()

    genera = duf3494_genera()
    already = done_genera()
    todo = genera[~genera.genus.isin(already)]
    if args.limit:
        todo = todo.head(args.limit)

    print(f"DUF3494 genera: {len(genera)}   already checked: {len(already)}   "
          f"to check: {len(todo)}\n")

    # Coordinates are a property of the collection site, not of the strain, so
    # one strain page per distinct locality is enough. Psychrobacter alone
    # returns 264 strains from a handful of sites; fetching every strain page
    # would be thousands of requests for the same few coordinates.
    site_coords = {}

    def coords_for(locality, strain_id):
        if locality not in site_coords:
            site_coords[locality] = strain_detail(strain_id)
            time.sleep(REQUEST_DELAY)
        return site_coords[locality]

    records = []
    for i, (_, row) in enumerate(todo.iterrows(), start=1):
        hits = search_genus(row.genus)
        time.sleep(REQUEST_DELAY)

        for h in hits:
            detail = coords_for(h["locality"], h["strain_id"])
            rec = {"query_genus": row.genus, "clade": row.clade,
                   "n_proteins": row.n_proteins, **h, **detail}
            rec["region"] = classify_region(h["locality"], detail.get("latitude"))
            rec["marine"] = is_marine(h["habitat"], h["locality"])
            records.append(rec)

        if hits:
            print(f"[{i:3d}/{len(todo)}] {row.genus:24s} {row.clade:9s} "
                  f"-> {len(hits):4d} strain(s), {len(set(h['locality'] for h in hits))} site(s)",
                  flush=True)
            # append incrementally so an interrupted run keeps its work
            pd.DataFrame(records).to_csv(
                OUT_STRAINS, mode="a", index=False, encoding="utf-8-sig",
                header=not os.path.exists(OUT_STRAINS))
            records = []
        elif i % 25 == 0:
            print(f"[{i:3d}/{len(todo)}] ... {row.genus} (no match)", flush=True)

    summarize()


def summarize():
    """Build the coordinate, site and genus-latitude tables from the strain table."""
    if not os.path.exists(OUT_STRAINS):
        print("no strain table yet")
        return

    df = pd.read_csv(OUT_STRAINS).drop_duplicates(subset=["strain_id"])
    print(f"\nPAMC strains matching DUF3494 genera: {len(df)}")
    print(f"  distinct genera matched: {df.query_genus.nunique()}")

    coords = df[df.latitude.notna() & df.longitude.notna()].copy()
    coords.to_csv(OUT_COORDS, index=False, encoding="utf-8-sig")
    print(f"  with usable coordinates: {len(coords)}  -> {OUT_COORDS}")

    if not coords.empty:
        print("\n  region breakdown (coordinate-bearing strains):")
        for region, n in coords.region.value_counts().items():
            print(f"    {region:12s} {n:4d}")
        print(f"\n  marine: {int(coords.marine.sum())}   "
              f"terrestrial/other: {int((~coords.marine).sum())}")
        print("\n  marine strains by region:")
        for region, n in coords[coords.marine].region.value_counts().items():
            print(f"    {region:12s} {n:4d}")

    sites = (coords.groupby(["locality", "latitude", "longitude"])
                   .agg(n_strains=("strain_id", "nunique"),
                        genera=("query_genus", lambda s: ", ".join(sorted(set(s)))),
                        region=("region", "first"),
                        marine=("marine", "first"))
                   .reset_index()
                   .sort_values("n_strains", ascending=False))
    sites.to_csv(OUT_SITES, index=False, encoding="utf-8-sig")
    print(f"\n  distinct collection sites: {len(sites)}  -> {OUT_SITES}")

    genus = (coords.groupby(["query_genus", "clade"])
                   .agg(n_strains=("strain_id", "nunique"),
                        pamc_lat_median=("latitude", "median"),
                        pamc_lat_min=("latitude", "min"),
                        pamc_lat_max=("latitude", "max"),
                        n_polar=("region", lambda s: int(s.isin(["Arctic", "Antarctic"]).sum())))
                   .reset_index())
    genus["pamc_abs_lat_median"] = genus.pamc_lat_median.abs().round(2)

    # Join the GBIF genus-level centre where the stage-2 batch has reached it,
    # so the two estimates can be compared directly.
    if os.path.exists(GBIF_ENV):
        gbif = pd.read_csv(GBIF_ENV)[["genus", "rep_lat", "habitat"]]
        gbif = gbif.rename(columns={"genus": "query_genus", "rep_lat": "gbif_abs_lat",
                                    "habitat": "gbif_habitat"})
        genus = genus.merge(gbif, on="query_genus", how="left")
        both = genus[genus.gbif_abs_lat.notna()]
        if not both.empty:
            genus["lat_gap"] = (genus.pamc_abs_lat_median - genus.gbif_abs_lat).round(2)
            print(f"\n  genera with both PAMC and GBIF estimates: {len(both)}")
            print(f"    median |lat| PAMC {both.pamc_abs_lat_median.median():.2f} "
                  f"vs GBIF {both.gbif_abs_lat.median():.2f}")

    genus.to_csv(OUT_GENUS, index=False, encoding="utf-8-sig")
    print(f"  genus latitude summary: {len(genus)} rows -> {OUT_GENUS}")


if __name__ == "__main__":
    main()
