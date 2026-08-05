"""
Stage-2 bulk download: the DUF3494 / Pfam PF11999 ice-binding-like family.

    python -m src.fetch_pfam

Downloads every UniProtKB entry carrying the PF11999 domain, with taxonomy,
as both TSV (metadata) and FASTA (sequences).

Query syntax note: UniProt's REST field is `xref:pfam-PF11999`. The
plausible-looking `xref_pfam:PF11999` is rejected with "not a valid search
field", and `database:pfam` matches every entry with any Pfam cross-
reference (~113 million) rather than this family.

Counts differ slightly between sources - InterPro reports 2322 proteins for
PF11999 while UniProtKB returns 2250 - because InterPro also counts entries
that UniProtKB has since merged or demoted. Use the UniProt figure, since
that is what the sequences come from.

The /stream endpoint returns the whole result set in one response (no
manual paging), but it is rate-limited; the files are cached on disk and
re-downloaded only with --force.
"""
import argparse
import os

import requests

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
PFAM_ID = "PF11999"
QUERY = f"xref:pfam-{PFAM_ID}"

TSV_FIELDS = ",".join([
    "accession", "reviewed", "protein_name", "organism_name", "organism_id",
    "lineage", "length", "sequence", "ft_domain", "cc_subcellular_location",
])

OUT_TSV = f"data/pfam{PFAM_ID[2:]}_proteins.tsv"
OUT_FASTA = f"data/pfam{PFAM_ID[2:]}_proteins.fasta"


def _download(fmt, out_path, force=False, timeout=600):
    if os.path.exists(out_path) and not force:
        size = os.path.getsize(out_path)
        print(f"[skip] {out_path} already exists ({size:,} bytes) - use --force to refetch")
        return out_path

    params = {"query": QUERY, "format": fmt}
    if fmt == "tsv":
        params["fields"] = TSV_FIELDS

    print(f"[get ] {fmt} stream for {QUERY} ...")
    with requests.get(STREAM, params=params, timeout=timeout, stream=True) as r:
        r.raise_for_status()
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 16):
                fh.write(chunk)

    print(f"[ok  ] {out_path} ({os.path.getsize(out_path):,} bytes)")
    return out_path


def summarize(tsv_path):
    """Print composition of the downloaded family - the EDA starting point."""
    import pandas as pd

    df = pd.read_csv(tsv_path, sep="\t")
    print(f"\nEntries: {len(df):,}")
    if "Reviewed" in df:
        print(f"  reviewed (Swiss-Prot): {(df.Reviewed == 'reviewed').sum():,}")
    if "Length" in df:
        print(f"  length: median={df.Length.median():.0f}  "
              f"range=[{df.Length.min()}, {df.Length.max()}]")

    if "Taxonomic lineage" in df:
        def top_clade(lin):
            for token in ("Bacteria", "Archaea", "Fungi", "Bacillariophyta",
                          "Viridiplantae", "Metazoa"):
                if isinstance(lin, str) and token in lin:
                    return token
            return "other/unclassified"

        counts = df["Taxonomic lineage"].map(top_clade).value_counts()
        print("\n  clade composition:")
        for clade, n in counts.items():
            print(f"    {clade:22s} {n:5,d}  ({100 * n / len(df):.1f}%)")

    if "Organism (ID)" in df:
        print(f"\n  distinct taxa: {df['Organism (ID)'].nunique():,}")
    return df


def main():
    ap = argparse.ArgumentParser(description=f"Download Pfam {PFAM_ID} (DUF3494) from UniProt")
    ap.add_argument("--force", action="store_true", help="re-download even if files exist")
    args = ap.parse_args()

    tsv = _download("tsv", OUT_TSV, force=args.force)
    _download("fasta", OUT_FASTA, force=args.force)
    summarize(tsv)

    print("\nNext: cluster to remove near-duplicates before feature extraction, e.g.")
    print(f"  cd-hit -i {OUT_FASTA} -o data/pfam{PFAM_ID[2:]}_nr90.fasta -c 0.9 -n 5")
    print("  (brew install cd-hit  /  conda install -c bioconda cd-hit)")


if __name__ == "__main__":
    main()
