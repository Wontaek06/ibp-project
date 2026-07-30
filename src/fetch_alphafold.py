"""
AlphaFold DB lookup: prediction metadata + per-residue pLDDT profile.

We use the *predicted* structure database (not folding from scratch) -
AFGP and other short/repetitive or heavily glycosylated sequences often
have no entry or a low-confidence model. That's expected, not a bug:
treat af=False / plddt=None as "no reliable structural feature available
for this species" rather than retrying.
"""
import io
import requests
import numpy as np

from src.cache import cached

AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"


@cached("alphafold_meta")
def alphafold_meta(accession):
    """UniProt accession -> {af, plddt, pdb_url}."""
    try:
        r = requests.get(AF_API.format(accession=accession), timeout=20)
    except requests.RequestException:
        return {"af": False, "plddt": None, "pdb_url": None}
    if not r.ok or not r.json():
        return {"af": False, "plddt": None, "pdb_url": None}
    entry = r.json()[0]
    return {
        "af": True,
        "plddt": entry.get("globalMetricValue"),
        "pdb_url": entry.get("pdbUrl"),
    }


@cached("alphafold_profile")
def fetch_plddt_profile(pdb_url):
    """
    Download a PDB file and return the per-residue pLDDT profile
    (AlphaFold stores per-residue pLDDT in the B-factor column).
    Requires biopython.
    """
    from Bio.PDB import PDBParser

    try:
        r = requests.get(pdb_url, timeout=30)
        r.raise_for_status()
    except requests.RequestException:
        return []

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", io.StringIO(r.text))

    profile = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    profile.append(residue["CA"].get_bfactor())
        break  # first model only
    return profile


def structural_regularity(plddt_profile):
    """
    Rough structural-regularity summary from a pLDDT profile:
      - plddt_std: how uniform the confidence is along the chain
      - plddt_periodicity: lag-2..7 autocorrelation peak, a heuristic proxy
        for repeat-unit structures (expected higher in AFGP / Type I).
    This is a simple heuristic for a preliminary report, not a validated
    structural-biology metric - state that plainly if used in the writeup.
    """
    if not plddt_profile or len(plddt_profile) < 10:
        return {"plddt_std": None, "plddt_periodicity": None}

    arr = np.array(plddt_profile, dtype=float)
    std = float(np.std(arr))

    centered = arr - arr.mean()
    autocorr = np.correlate(centered, centered, mode="full")
    autocorr = autocorr[len(autocorr) // 2:]
    denom = autocorr[0] if autocorr[0] != 0 else 1.0
    autocorr = autocorr / denom

    periodicity = float(np.max(autocorr[2:8])) if len(autocorr) > 8 else None

    return {
        "plddt_std": round(std, 2),
        "plddt_periodicity": round(periodicity, 3) if periodicity is not None else None,
    }
