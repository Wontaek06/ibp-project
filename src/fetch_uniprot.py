"""
UniProt fetching for antifreeze proteins.

IMPORTANT: a naive full-text search for "antifreeze" pulls in unrelated
proteins whose description happens to mention the word (confirmed cases:
Elongation Factor-1a, Cytochrome c oxidase I). Restricting to
`protein_name:antifreeze` AND `reviewed:true` (Swiss-Prot only) fixes most
of this, but ALWAYS check the returned prot_name before trusting a hit -
some species genuinely have no curated AFP entry and should be excluded
rather than matched to the wrong protein.
"""
import requests

from src.cache import cached

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


@cached("uniprot")
def uniprot_afp(species_name, size=5):
    """
    species name -> {uniprot, prot_name, seq} for the first reviewed hit
    whose protein name contains "antifreeze", or None if no reviewed entry.
    """
    query = f'(protein_name:antifreeze) AND organism_name:"{species_name}" AND reviewed:true'
    try:
        r = requests.get(
            UNIPROT_SEARCH,
            params={
                "query": query,
                "format": "json",
                "fields": "accession,protein_name,organism_name,sequence",
                "size": size,
            },
            timeout=20,
        ).json()
    except requests.RequestException as e:
        print(f"[{species_name}] UniProt request failed: {e}")
        return None

    results = r.get("results", [])
    if not results:
        return None

    entry = results[0]
    pdd = entry.get("proteinDescription", {}) or {}
    try:
        name = (pdd.get("recommendedName") or pdd.get("submissionNames", [{}])[0]) \
            .get("fullName", {}).get("value", "")
    except Exception:
        name = ""

    return {
        "uniprot": entry.get("primaryAccession"),
        "prot_name": name,
        "seq": (entry.get("sequence") or {}).get("value"),
    }
