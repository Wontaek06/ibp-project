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
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}.json"


@cached("uniprot_acc")
def uniprot_by_accession(accession):
    """
    Accession -> {uniprot, prot_name, organism, taxon_id, lineage, seq, reviewed}.

    Used instead of a keyword search when the seed set is hand-curated.
    Keyword search is what produced the false positives recorded in
    data/seed_species.csv (Elongation Factor-1a, Cytochrome c oxidase);
    pinning the accession removes that failure mode entirely, at the cost
    of the curation step being manual.
    """
    try:
        r = requests.get(UNIPROT_ENTRY.format(accession=accession), timeout=25)
        r.raise_for_status()
        e = r.json()
    except (requests.RequestException, ValueError) as err:
        print(f"[{accession}] UniProt entry fetch failed: {err}")
        return None

    pdd = e.get("proteinDescription", {}) or {}
    named = pdd.get("recommendedName") or (pdd.get("submissionNames") or [{}])[0] or {}
    org = e.get("organism", {}) or {}

    return {
        "uniprot": e.get("primaryAccession"),
        "prot_name": (named.get("fullName", {}) or {}).get("value", ""),
        "organism": org.get("scientificName"),
        "taxon_id": org.get("taxonId"),
        "lineage": org.get("lineage", []),
        "seq": (e.get("sequence") or {}).get("value"),
        "reviewed": e.get("entryType", "").startswith("UniProtKB reviewed"),
    }


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
