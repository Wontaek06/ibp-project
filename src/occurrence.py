"""
Occurrence retrieval by resolved taxon, with an OBIS fallback.

GBIF is queried by usageKey rather than by name (see src/taxonomy.py for
why name-based matching fails on microbial taxa). OBIS is tried second
because it holds marine records - including microbial metabarcoding
surveys - that GBIF's backbone sometimes does not surface for a genus.

`cold_points` deliberately selects the most poleward records within the
local winter season: annual means wash out the cold extreme that actually
drives ice-binding protein expression.
"""
import time

import requests

from src.cache import cached

GBIF_OCC = "https://api.gbif.org/v1/occurrence/search"
OBIS_OCC = "https://api.obis.org/v3/occurrence"


@cached("gbif_occ_key")
def gbif_occurrences(usage_key, limit=300, request_pause=0.2):
    """GBIF usageKey -> list of {lat, lon, month, source} (coordinates only)."""
    records, offset = [], 0
    while offset < limit:
        try:
            occ = requests.get(
                GBIF_OCC,
                params={
                    "taxonKey": usage_key,
                    "hasCoordinate": "true",
                    "hasGeospatialIssue": "false",
                    "limit": min(300, limit - offset),
                    "offset": offset,
                },
                timeout=30,
            ).json()
        except (requests.RequestException, ValueError) as e:
            print(f"    [gbif] occurrence request failed: {e}")
            break

        for r in occ.get("results", []):
            lat, lon = r.get("decimalLatitude"), r.get("decimalLongitude")
            if lat is not None and lon is not None:
                records.append({"lat": lat, "lon": lon,
                                "month": r.get("month"), "source": "GBIF"})

        if occ.get("endOfRecords", True):
            break
        offset += 300
        time.sleep(request_pause)

    return records


@cached("obis_occ")
def obis_occurrences(scientific_name, size=300):
    """Scientific name -> OBIS records as {lat, lon, month, source}."""
    try:
        r = requests.get(
            OBIS_OCC,
            params={"scientificname": scientific_name, "size": size},
            timeout=40,
        ).json()
    except (requests.RequestException, ValueError) as e:
        print(f"    [obis] request failed: {e}")
        return []

    records = []
    for rec in r.get("results", []):
        lat, lon = rec.get("decimalLatitude"), rec.get("decimalLongitude")
        if lat is None or lon is None:
            continue
        month = rec.get("month")
        if month is None and rec.get("eventDate"):
            # eventDate is ISO-ish; month is the second field when present
            bits = str(rec["eventDate"])[:10].split("-")
            month = int(bits[1]) if len(bits) > 1 and bits[1].isdigit() else None
        records.append({"lat": lat, "lon": lon, "month": month, "source": "OBIS"})
    return records


def fetch_occurrences(taxon, min_records=5):
    """
    Resolved taxon dict -> (records, provider). Tries GBIF, then OBIS.

    OBIS is consulted whenever GBIF returns fewer than `min_records`, not
    only when it returns zero: a handful of stray points gives a
    meaningless habitat temperature.
    """
    records, provider = [], None

    if taxon.get("usageKey"):
        records = gbif_occurrences(taxon["usageKey"])
        provider = "GBIF" if records else None

    if len(records) < min_records and taxon.get("name"):
        obis = obis_occurrences(taxon["name"])
        if len(obis) > len(records):
            records, provider = obis, "OBIS"

    return records, provider


def dedupe_points(records, precision=1):
    """
    Collapse records that share a coordinate cell (default ~0.1 deg).

    Repeat sampling at one station is common - Colwellia returned the same
    coordinate four times in five selected points, which pins the habitat
    temperature to that single station.
    """
    seen, out = set(), []
    for r in records:
        cell = (round(r["lat"], precision), round(r["lon"], precision))
        if cell in seen:
            continue
        seen.add(cell)
        out.append(r)
    return out


def core_points(records, k=5):
    """
    Localities at the centre of a taxon's latitudinal distribution.

    This is the default representative locality, not cold_points(). The
    poleward tail turned out to be unusable as a habitat descriptor for
    wide-ranging taxa no matter how it was trimmed: Thunnus albacares has
    a median occurrence latitude of 7.7N but 29,118 GBIF records above
    40N, so any tail-based rule reports a tropical tuna sitting in ~2 C
    water. The tail measures range limits and sampling effort; the centre
    measures where the animal actually lives.

    The cold extreme is not lost by centring, because the environmental
    layer supplies it: thetao_min is the annual minimum for the chosen
    cell, so a 45N fish still gets its sub-zero winter value.

    Latitude is compared as |lat| so that bipolar taxa (several of the
    sea-ice bacteria) are not split across hemispheres; temperature is
    approximately symmetric about the equator, so this is safe here.

    The centre is computed on the RAW records and only the returned
    localities are deduplicated. Deduplicating before taking the median
    discards record density, which is the only evidence of where a taxon is
    concentrated: Colwellia's 228 Svalbard records collapse to a handful of
    cells, moving its centre from 78.7N to 56.6N and its habitat
    temperature from -2 C to +1.2 C. Leucosporidium moved 66.4N -> 45.8N
    the same way. Density is biased by sampling effort, but collapsing it
    is a heavier distortion than leaving it in.
    """
    if not records:
        return []
    abs_lats = sorted(abs(r["lat"]) for r in records)
    median_abs = abs_lats[len(abs_lats) // 2]
    pts = dedupe_points(records)
    return sorted(pts, key=lambda r: abs(abs(r["lat"]) - median_abs))[:k]


def cold_points(records, k=5, drop_extreme=0.02, tail_frac=0.15):
    """
    Localities representing the poleward tail of a taxon's range, with the
    most extreme records trimmed.

    Use core_points() for habitat temperature. This function describes a
    range *limit*, which is a different quantity and is much more sensitive
    to sampling effort - keep it for "how far poleward does this taxon get"
    questions, not for "how cold is its habitat".

    Returns up to k distinct localities drawn from the poleward tail after
    discarding the top `drop_extreme` fraction. The trim matters because
    the extreme tail is dominated by vagrants and sampling artefacts rather
    than habitat: Lutjanus campechanus (median 27.9N, a Gulf of Mexico
    snapper) contributed five records near 44N to a 300-record sample even
    though GBIF holds only 18 records above 40N for the species at all.
    Taking the single most poleward points therefore reported a tropical
    fish as living at 1.98 C. Polar taxa are unaffected by the trim -
    Colwellia stays above 80N either way.

    NOTE - do not reintroduce an observation-month filter here. Selecting
    records observed in winter looks right but inverts the intended effect,
    because polar fieldwork happens in summer: Colwellia has 231/300
    records in March and only 7 in Dec-Feb, and those 7 sit at 56.6N and
    41.3S, dragging a genuinely Svalbard-centred taxon (median 78.7N) down
    to a 3.5 C representative cell. Leucosporidium failed the same way,
    landing on 42S / 11.6 C.

    Seasonality belongs in the environmental layer, not the occurrence
    filter: Bio-ORACLE's thetao_min is already the annual minimum, so the
    coldest month is captured for whatever locality is chosen. Filtering
    occurrences by month double-counts season and imports sampling-effort
    bias on top.
    """
    if not records:
        return []

    pts = dedupe_points(records)
    pts.sort(key=lambda r: abs(r["lat"]), reverse=True)

    start = int(len(pts) * drop_extreme)
    tail_end = max(start + k, int(len(pts) * tail_frac))
    pool = pts[start:tail_end] or pts[:k]
    return pool[:k]
