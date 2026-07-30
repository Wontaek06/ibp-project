"""
GBIF occurrence fetching, restricted to winter / most-poleward records.

Rationale: annual-mean environmental values wash out the cold season that
actually triggers antifreeze protein expression, and temperate species need
their winter extreme captured, not their summer average.
"""
import time
import requests

from src.cache import cached

GBIF_MATCH = "https://api.gbif.org/v1/species/match"
GBIF_OCC = "https://api.gbif.org/v1/occurrence/search"


@cached("gbif")
def gbif_records(species_name, limit=300, request_pause=0.2):
    """species name -> list of {lat, lon, month} dicts (quality-filtered)."""
    try:
        m = requests.get(GBIF_MATCH, params={"name": species_name}, timeout=15).json()
    except requests.RequestException as e:
        print(f"[{species_name}] match request failed: {e}")
        return []

    key = m.get("usageKey")
    if key is None:
        print(f"[{species_name}] no taxonKey (matchType={m.get('matchType')})")
        return []

    records, offset = [], 0
    while offset < limit:
        try:
            occ = requests.get(
                GBIF_OCC,
                params={
                    "taxonKey": key,
                    "hasCoordinate": "true",
                    "hasGeospatialIssue": "false",
                    "limit": min(300, limit - offset),
                    "offset": offset,
                },
                timeout=25,
            ).json()
        except requests.RequestException as e:
            print(f"[{species_name}] occurrence request failed: {e}")
            break

        for r in occ.get("results", []):
            lat, lon = r.get("decimalLatitude"), r.get("decimalLongitude")
            if lat is not None and lon is not None:
                records.append({"lat": lat, "lon": lon, "month": r.get("month")})

        if occ.get("endOfRecords", True):
            break
        offset += 300
        time.sleep(request_pause)

    return records


def cold_points(records, k=3):
    """
    Pick winter-season records for the hemisphere inferred from median
    latitude; fall back to all records if too few winter observations exist.
    Returns the k most poleward points among the chosen subset.
    """
    if not records:
        return []
    med_lat = sorted(r["lat"] for r in records)[len(records) // 2]
    winter_months = {6, 7, 8} if med_lat < 0 else {12, 1, 2}
    winter = [r for r in records if r.get("month") in winter_months]
    pool = winter if len(winter) >= 3 else records
    return sorted(pool, key=lambda r: abs(r["lat"]), reverse=True)[:k]


def species_cold_points(species_name, k=3, occ_limit=300):
    """Convenience wrapper: species name -> up to k cold/poleward points."""
    recs = gbif_records(species_name, limit=occ_limit)
    return cold_points(recs, k=k)
