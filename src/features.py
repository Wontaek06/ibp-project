"""
Sequence-derived physicochemical features (Biopython ProtParam).

Biological rationale for each feature relative to AFP biology:
  - ala_pct: Type I AFPs are alanine-rich amphipathic alpha-helices;
    AFGP repeat units (Thr-Ala-Ala) are also Ala-rich - this is the
    strongest, most reliable discriminating signal for AFGP.
  - thr_pct: AFGP repeat units are built on Thr-Ala-Ala / Thr-Pro-Ala.
  - gravy:   IMPORTANT - GRAVY is computed from the amino acid sequence
    only and does NOT account for post-translational glycosylation.
    AFGP is biologically hydrophilic in vivo because of its attached
    sugar chains, but its peptide backbone is Ala-rich and therefore
    reads as HYDROPHOBIC (positive GRAVY) by this metric - confirmed
    empirically on all 4 validated AFGP sequences (GRAVY +0.57 to +0.99).
    Do not treat a positive GRAVY on a real AFGP as a data error; do not
    use GRAVY as evidence of AFGP identity in the report - use ala_pct/
    thr_pct instead, and note the glycosylation blind spot as a stated
    limitation.
  - pI, instability, aromaticity: general descriptors, useful mainly for
    exploratory grouping rather than a specific AFP hypothesis.
"""
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def seq_features(seq):
    if not seq or len(seq) < 5:
        return {
            "length": None, "ala_pct": None, "thr_pct": None, "pro_pct": None,
            "pI": None, "gravy": None, "instability": None, "aromaticity": None,
        }

    length = len(seq)
    ala_pct = round(100 * seq.count("A") / length, 1)
    thr_pct = round(100 * seq.count("T") / length, 1)
    pro_pct = round(100 * seq.count("P") / length, 1)

    try:
        # ProtParam chokes on non-standard residues (X, U); strip/substitute
        pa = ProteinAnalysis(seq.replace("X", "").replace("U", "C"))
        pI = round(pa.isoelectric_point(), 2)
        gravy = round(pa.gravy(), 3)
        instability = round(pa.instability_index(), 1)
        aromaticity = round(pa.aromaticity(), 3)
    except Exception:
        pI = gravy = instability = aromaticity = None

    return {
        "length": length, "ala_pct": ala_pct, "thr_pct": thr_pct, "pro_pct": pro_pct,
        "pI": pI, "gravy": gravy, "instability": instability, "aromaticity": aromaticity,
    }
