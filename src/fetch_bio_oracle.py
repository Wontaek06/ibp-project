"""
Bio-ORACLE v3 (ERDDAP) point extraction.

Includes a radius search fallback: a point-exact query on a coastal/land
grid cell returns null across every layer simultaneously (confirmed
empirically - Pagothenia borchgrevinki, Gadus morhua, Hemitripterus
americanus, Osmerus mordax all failed this way before the fix). Searching
a small bounding box around the point and taking the nearest non-null cell
resolves it.

Run `list_variables` / `search_datasets` once to confirm dataset IDs and
variable names before trusting CONFIG below - Bio-ORACLE dataset naming
has changed across versions.
"""
import math
import requests

from src.cache import cached

ERDDAP = "https://erddap.bio-oracle.org/erddap"

# label -> (dataset_id, variable_name). Verify with list_variables() first.
CONFIG = {
    "surf_min_temp": ("thetao_baseline_2000_2019_depthsurf", "thetao_min"),
    "bottom_temp": ("thetao_baseline_2000_2019_depthmax", "thetao_mean"),
    "sea_ice_thick": ("sithick_baseline_2000_2020_depthsurf", "sithick_mean"),
}


def list_variables(dataset_id):
    r = requests.get(f"{ERDDAP}/info/{dataset_id}/index.json", timeout=30).json()
    return [row[1] for row in r["table"]["rows"] if row[0] == "variable"]


def search_datasets(keyword):
    r = requests.get(f"{ERDDAP}/search/index.json", params={"searchFor": keyword}, timeout=30).json()
    cols = r["table"]["columnNames"]
    idx = cols.index("Dataset ID")
    return [row[idx] for row in r["table"]["rows"]]


@cached("bio_oracle")
def bo_point(dataset_id, variable, lat, lon, radius_deg=0.5):
    """Value at (lat, lon); falls back to nearest non-null cell within radius_deg."""
    lat_lo, lat_hi = lat - radius_deg, lat + radius_deg
    lon_lo, lon_hi = lon - radius_deg, lon + radius_deg

    for time_part in ("[last]", ""):
        query = f"{variable}{time_part}[({lat_lo}):({lat_hi})][({lon_lo}):({lon_hi})]"
        try:
            r = requests.get(f"{ERDDAP}/griddap/{dataset_id}.json?{query}", timeout=40)
        except requests.RequestException:
            continue
        if not r.ok:
            continue
        try:
            rows = r.json()["table"]["rows"]
        except Exception:
            continue
        if not rows:
            continue

        best_val, best_dist = None, None
        for row in rows:
            val = row[-1]
            if val is None:
                continue
            row_lat, row_lon = row[-3], row[-2]
            d = math.hypot(row_lat - lat, row_lon - lon)
            if best_dist is None or d < best_dist:
                best_dist, best_val = d, val
        if best_val is not None:
            return round(float(best_val), 3)

    return None


def bo_species(dataset_id, variable, points, radius_deg=0.5):
    """Average a variable over a list of {lat, lon} points, skipping nulls."""
    vals = [bo_point(dataset_id, variable, p["lat"], p["lon"], radius_deg) for p in points]
    vals = [v for v in vals if v is not None]
    return round(sum(vals) / len(vals), 3) if vals else None


def fetch_all_env_vars(points, config=None):
    """points -> {surf_min_temp, bottom_temp, sea_ice_thick} using CONFIG."""
    cfg = config or CONFIG
    return {label: bo_species(dsid, var, points) for label, (dsid, var) in cfg.items()}
