"""
UniProt organism name -> GBIF/OBIS taxon resolution.

This layer exists because UniProt and GBIF do not speak the same name
dialect, and the mismatch silently produces empty occurrence sets:

  * UniProt carries strain junk GBIF has never heard of
    ("Flavobacterium frigoris (strain PS1)", "Leucosporidium sp. (strain AY30)").
  * Many characterized microbial IBPs come from undescribed species
    ("Colwellia sp"), so only the genus is resolvable - which is fine,
    genus-level occurrence is enough for a habitat-temperature signal.
  * GBIF's /species/match returns matchType=NONE for a bare microbial
    genus because the same name exists in several backbone kingdoms
    (confirmed: "Colwellia" -> "Multiple equal matches"). Passing a
    kingdom hint resolves it to usageKey 3222745. The hint is derived
    from the UniProt lineage, so it costs nothing extra.

Without all three fixes every microbial/algal IBP drops out of the
pipeline at the occurrence step and only the fish survive.
"""
import re

import requests

from src.cache import cached

GBIF_MATCH = "https://api.gbif.org/v1/species/match"

# UniProt lineage token -> GBIF backbone kingdom. Order matters: the first
# token found in the lineage wins, so put the specific clades first.
LINEAGE_TO_KINGDOM = [
    ("Metazoa", "Animalia"),
    ("Fungi", "Fungi"),
    ("Viridiplantae", "Plantae"),
    ("Bacillariophyta", "Chromista"),   # diatoms
    ("Stramenopiles", "Chromista"),
    ("Sar", "Chromista"),
    ("Bacteria", "Bacteria"),
    ("Archaea", "Archaea"),
]

# Tokens that mean "this is not a real species epithet"
_PLACEHOLDER_EPITHETS = {"sp", "spp", "cf", "aff", "bacterium", "symbiont"}


def kingdom_from_lineage(lineage):
    """UniProt lineage list -> GBIF kingdom string (or None)."""
    if not lineage:
        return None
    for token, kingdom in LINEAGE_TO_KINGDOM:
        if token in lineage:
            return kingdom
    return None


def clean_organism(name):
    """
    UniProt organism name -> (query_name, rank).

    Strips strain/isolate annotations and falls back to the genus when the
    species epithet is a placeholder:
        "Flavobacterium frigoris (strain PS1)" -> ("Flavobacterium frigoris", "SPECIES")
        "Leucosporidium sp. (strain AY30)"     -> ("Leucosporidium", "GENUS")
        "Colwellia sp"                         -> ("Colwellia", "GENUS")
    """
    name = re.sub(r"\s*\((?:strain|isolate)[^)]*\)", "", name, flags=re.I)
    name = re.sub(r"\s+", " ", name).strip()

    parts = name.split()
    if not parts:
        return None, None
    if len(parts) == 1:
        return parts[0], "GENUS"

    epithet = parts[1].rstrip(".").lower()
    if epithet in _PLACEHOLDER_EPITHETS:
        return parts[0], "GENUS"
    return f"{parts[0]} {parts[1]}", "SPECIES"


@cached("gbif_match")
def gbif_match(name, kingdom=None, rank=None):
    """Name (+ optional hints) -> GBIF backbone match dict, or {} on failure."""
    params = {"name": name}
    if kingdom:
        params["kingdom"] = kingdom
    if rank:
        params["rank"] = rank
    try:
        return requests.get(GBIF_MATCH, params=params, timeout=20).json()
    except (requests.RequestException, ValueError):
        return {}


def resolve_taxon(organism, lineage=None):
    """
    UniProt organism (+ lineage) -> {name, rank, usageKey, matchType, kingdom}.

    Tries species-level first, then falls back to the genus. usageKey is
    None when nothing resolved - callers should treat that as "no GBIF
    occurrence path" rather than retrying.
    """
    kingdom = kingdom_from_lineage(lineage)
    query, rank = clean_organism(organism)
    if not query:
        return {"name": None, "rank": None, "usageKey": None,
                "matchType": "NONE", "kingdom": kingdom}

    attempts = [(query, rank)]
    if rank == "SPECIES":
        attempts.append((query.split()[0], "GENUS"))  # genus fallback

    for attempt_name, attempt_rank in attempts:
        m = gbif_match(attempt_name, kingdom=kingdom, rank=attempt_rank)
        if m.get("usageKey"):
            return {
                "name": attempt_name,
                "rank": m.get("rank", attempt_rank),
                "usageKey": m["usageKey"],
                "matchType": m.get("matchType"),
                "kingdom": m.get("kingdom") or kingdom,
            }

    return {"name": query, "rank": rank, "usageKey": None,
            "matchType": "NONE", "kingdom": kingdom}
