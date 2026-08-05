"""
Terrestrial minimum air temperature from Open-Meteo's ERA5 archive.

Bio-ORACLE is marine-only, so land taxa come back null on every layer -
confirmed for Leucosporidium sp. AY30 (Arctic yeast) and Typhula
ishikariensis (snow mould), both of which resolve occurrence coordinates
normally and then hit an empty ocean grid. Those two are ice-binding
proteins from genuinely cold habitats; dropping them would bias the set
toward marine organisms.

IMPORTANT - `land_min_temp` is 2 m AIR temperature and is NOT comparable
to Bio-ORACLE's `surf_min_temp` (sea surface water). Water under ice is
bounded near -2 C by freezing point, while winter air over land routinely
reaches -20 C or below. Keep them in separate columns; never pool them
into one "temperature" variable for a statistical test.

WorldClim was the alternative but ships as a ~300 MB global raster per
variable, which is a heavy dependency for a handful of points. Open-Meteo
needs no key and answers per coordinate.
"""
import requests

from src.cache import cached

ARCHIVE = "https://archive-api.open-meteo.com/v1/archive"

# ERA5 reanalysis window. Kept short deliberately: each request returns one
# value per day, so a 20-year pull is ~7300 numbers per coordinate for no
# extra signal in a winter-minimum statistic.
START, END = "2010-01-01", "2019-12-31"

# Percentile used as the "annual minimum" proxy. The absolute daily
# minimum over a decade is a single freak night; the 1st percentile is the
# recurring cold extreme the organism actually has to survive.
COLD_PERCENTILE = 1


@cached("open_meteo")
def land_daily_min(lat, lon, start=START, end=END):
    """Coordinate -> list of daily 2 m minimum air temperatures (deg C)."""
    try:
        r = requests.get(
            ARCHIVE,
            params={
                "latitude": lat, "longitude": lon,
                "start_date": start, "end_date": end,
                "daily": "temperature_2m_min", "timezone": "UTC",
            },
            timeout=60,
        )
        r.raise_for_status()
        payload = r.json()
    except (requests.RequestException, ValueError) as e:
        print(f"    [open-meteo] {lat},{lon} failed: {e}")
        return []

    if payload.get("error"):
        print(f"    [open-meteo] {lat},{lon}: {payload.get('reason')}")
        return []
    return [v for v in payload.get("daily", {}).get("temperature_2m_min", []) if v is not None]


def land_min_temp(points, percentile=COLD_PERCENTILE):
    """
    Points -> recurring winter minimum air temperature, averaged over points.

    Returns None when no point yields data, so callers can distinguish
    "not a land taxon / no data" from a real cold value.
    """
    values = []
    for p in points:
        daily = land_daily_min(p["lat"], p["lon"])
        if not daily:
            continue
        daily = sorted(daily)
        idx = max(0, min(len(daily) - 1, len(daily) * percentile // 100))
        values.append(daily[idx])
    return round(sum(values) / len(values), 3) if values else None
