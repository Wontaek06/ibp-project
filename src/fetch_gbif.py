"""
GBIF occurrence fetching by species name, reduced to poleward localities.

Rationale: annual-mean environmental values wash out the cold season that
actually triggers antifreeze protein expression. That cold extreme is
captured by the environmental layer (Bio-ORACLE thetao_min, an annual
minimum), NOT by filtering occurrences to winter observation months - see
src/occurrence.cold_points for why the latter backfires.

Point selection is re-exported from src/occurrence.py so this module and
the stage-1 spike cannot drift apart.
"""
import time
import requests

from src.cache import cached
from src.occurrence import cold_points  # noqa: F401  (re-exported)

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




def species_cold_points(species_name, k=3, occ_limit=300):
    """Convenience wrapper: species name -> up to k cold/poleward points."""
    recs = gbif_records(species_name, limit=occ_limit)
    return cold_points(recs, k=k)
