"""
End-to-end pipeline: seed species -> environment + sequence + structure
features -> final_merged.csv + figures.

Run from the project root:
    python -m src.pipeline

Takes a few minutes due to API rate limits (GBIF, Bio-ORACLE, UniProt,
AlphaFold). Safe to re-run; each step prints progress per species.
"""
import time
import pandas as pd
import matplotlib.pyplot as plt

from src.fetch_gbif import species_cold_points
from src.fetch_bio_oracle import fetch_all_env_vars
from src.fetch_uniprot import uniprot_afp
from src.fetch_alphafold import alphafold_meta, fetch_plddt_profile, structural_regularity
from src.features import seq_features
from src.stats_model import group_comparison, posthoc_dunn, classify_afp_type, bootstrap_ci

SEED_PATH = "data/seed_species.csv"
OUT_PATH = "data/final_merged.csv"
FIG_DIR = "figures"

AFP_TYPES = ["AFGP", "Type I", "Type II", "Type III"]  # Type IV excluded from group plots (n=1)


def build_environment_table(seed_df):
    """All 18 seed species: GBIF winter points -> Bio-ORACLE env vars."""
    rows = []
    for _, s in seed_df.iterrows():
        name = s["scientific_name"]
        points = species_cold_points(name, k=3)
        env = fetch_all_env_vars(points) if points else {
            "surf_min_temp": None, "bottom_temp": None, "sea_ice_thick": None
        }
        rep_lat = sorted(p["lat"] for p in points)[len(points) // 2] if points else None
        rows.append({"scientific_name": name, "rep_lat": round(rep_lat, 2) if rep_lat else None, **env})
        print(f"[env] {name:32s} lat={rep_lat} {env}")
        time.sleep(0.2)
    return pd.DataFrame(rows)


def build_sequence_structure_table(seed_df):
    """
    Only species with validation_status == 'verified' or 'RELABEL' are
    queried for sequence/structure - species marked EXCLUDE in
    seed_species.csv failed UniProt validation earlier (see notes column)
    and are intentionally skipped here rather than re-fetched.
    """
    rows = []
    for _, s in seed_df.iterrows():
        name, status = s["scientific_name"], s["validation_status"]
        if status == "EXCLUDE":
            rows.append({"scientific_name": name, "uniprot": None, "seq": None})
            continue

        up = uniprot_afp(name)
        if up is None:
            rows.append({"scientific_name": name, "uniprot": None, "seq": None})
            print(f"[seq] {name:32s} no reviewed UniProt entry")
            continue

        feats = seq_features(up["seq"])
        af = alphafold_meta(up["uniprot"])
        struct = {"plddt_std": None, "plddt_periodicity": None}
        if af["af"] and af["pdb_url"]:
            profile = fetch_plddt_profile(af["pdb_url"])
            struct = structural_regularity(profile)

        rows.append({
            "scientific_name": name, "uniprot": up["uniprot"], "prot_name": up["prot_name"],
            "seq": up["seq"], **feats, "af": af["af"], "plddt": af["plddt"], **struct,
        })
        print(f"[seq] {name:32s} {up['uniprot']} gravy={feats['gravy']} af={af['af']}")
        time.sleep(0.3)
    return pd.DataFrame(rows)


def make_figures(merged, seq_cols_present):
    plt.rcParams.update({"font.size": 10})

    # Environment figures (use full 18-species environment data)
    d = merged.dropna(subset=["surf_min_temp"])
    colors = {"AFGP": "#1f77b4", "Type I": "#ff7f0e", "Type II": "#2ca02c",
              "Type III": "#d62728", "Type IV": "#9467bd"}
    plt.figure(figsize=(7, 5))
    for t, c in colors.items():
        sub = d[d.afp_type == t]
        if len(sub):
            plt.scatter(sub.rep_lat, sub.surf_min_temp, c=c, label=t, s=70, edgecolor="k")
    plt.axhline(0, ls="--", c="gray", lw=0.8)
    plt.xlabel("Representative poleward latitude (deg)")
    plt.ylabel("Surface MIN temp (C, winter)")
    plt.title("Cold-habitat winter temperature vs latitude")
    plt.legend(); plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/env_latitude_vs_sst.png", dpi=150); plt.close()

    plt.figure(figsize=(7, 5))
    plt.boxplot([merged[merged.afp_type == t].bottom_temp.dropna() for t in AFP_TYPES],
                tick_labels=AFP_TYPES)
    plt.ylabel("Bottom mean temp (C)")
    plt.title("Bottom habitat temperature by AFP type")
    plt.tight_layout(); plt.savefig(f"{FIG_DIR}/env_bottom_temp_by_type.png", dpi=150); plt.close()

    plt.figure(figsize=(7, 5))
    plt.boxplot([merged[merged.afp_type == t].sea_ice_thick.dropna() for t in AFP_TYPES],
                tick_labels=AFP_TYPES)
    plt.ylabel("Sea ice thickness (m)")
    plt.title("Sea ice thickness by AFP type")
    plt.tight_layout(); plt.savefig(f"{FIG_DIR}/env_seaice_by_type.png", dpi=150); plt.close()

    # Sequence figures (only species with validated sequence features)
    if seq_cols_present:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        for ax, col, ylab in zip(axes, ["gravy", "pI", "instability"],
                                  ["GRAVY (hydrophobicity)", "Isoelectric point (pI)", "Instability index"]):
            ax.boxplot([merged[merged.afp_type == t][col].dropna() for t in AFP_TYPES],
                       tick_labels=AFP_TYPES)
            ax.set_ylabel(ylab)
            ax.set_title(f"{ylab} by AFP type")
        plt.tight_layout(); plt.savefig(f"{FIG_DIR}/seq_physicochem_by_type.png", dpi=150); plt.close()

    # Structure figures (AlphaFold pLDDT-derived; species without a model are dropped)
    struct_cols = [c for c in ("plddt", "plddt_std", "plddt_periodicity") if c in merged]
    if struct_cols and merged[struct_cols[0]].notna().sum() >= 3:
        labels = {"plddt": "Mean pLDDT (model confidence)",
                  "plddt_std": "pLDDT std (profile uniformity)",
                  "plddt_periodicity": "pLDDT autocorr peak (repeat proxy)"}
        fig, axes = plt.subplots(1, len(struct_cols), figsize=(5 * len(struct_cols), 5))
        if len(struct_cols) == 1:
            axes = [axes]
        for ax, col in zip(axes, struct_cols):
            ax.boxplot([merged[merged.afp_type == t][col].dropna() for t in AFP_TYPES],
                       tick_labels=AFP_TYPES)
            ax.set_ylabel(labels[col])
            ax.set_title(f"{labels[col]} by AFP type")
        plt.tight_layout(); plt.savefig(f"{FIG_DIR}/struct_plddt_by_type.png", dpi=150); plt.close()

    print(f"Figures saved to {FIG_DIR}/")


def main():
    seed_df = pd.read_csv(SEED_PATH)

    env_df = build_environment_table(seed_df)
    seq_df = build_sequence_structure_table(seed_df)

    # afp_type in seed_df already reflects the Gadus morhua -> Type IV relabel
    merged = seed_df.merge(env_df, on="scientific_name", how="left") \
                     .merge(seq_df, on="scientific_name", how="left")
    merged.to_csv(OUT_PATH, index=False, encoding="utf-8-sig")
    print(f"\nSaved {OUT_PATH} ({len(merged)} species)")

    seq_cols_present = merged["gravy"].notna().sum() >= 3 if "gravy" in merged else False
    make_figures(merged, seq_cols_present)

    env_cols = ["surf_min_temp", "bottom_temp", "sea_ice_thick"]
    seq_stat_cols = ["ala_pct", "thr_pct", "gravy", "pI", "instability"] if seq_cols_present else []
    struct_stat_cols = [c for c in ("plddt", "plddt_std", "plddt_periodicity")
                        if c in merged and merged[c].notna().sum() >= 3]

    print("\n--- Kruskal-Wallis + bootstrap CI ---")
    stat_rows = []
    for label, cols in (("environment", env_cols), ("sequence", seq_stat_cols),
                        ("structure", struct_stat_cols)):
        for col in cols:
            res = group_comparison(merged, col)
            if res:
                res["block"] = label
                stat_rows.append(res)
                print(res)
                for t in AFP_TYPES:
                    vals = merged[merged.afp_type == t][col].dropna().values
                    lo, hi = bootstrap_ci(vals)
                    if lo is not None:
                        print(f"    {t:9s} median CI95 = [{lo}, {hi}]  (n={len(vals)})")

    if stat_rows:
        pd.DataFrame(stat_rows).to_csv("data/stats_summary.csv", index=False, encoding="utf-8-sig")
        print("Saved data/stats_summary.csv")

    print("\n--- Dunn's post-hoc (pairwise, Holm-corrected) ---")
    for col in env_cols + seq_stat_cols + struct_stat_cols:
        ph = posthoc_dunn(merged, col)
        if ph is not None:
            print(f"\n[{col}]")
            print(ph.round(4))

    print("\n--- Feature importance (RandomForest, LOO-CV) ---")
    feats = [c for c in ["gravy", "pI", "ala_pct", "thr_pct", "plddt", "plddt_std",
                         "surf_min_temp", "bottom_temp", "sea_ice_thick"] if c in merged]
    acc, importances = classify_afp_type(merged, feats)
    print(f"LOO accuracy: {acc}")
    print(f"Feature importances: {importances}")

    if importances:
        plt.figure(figsize=(7, 5))
        items = sorted(importances.items(), key=lambda kv: kv[1])
        plt.barh([k for k, _ in items], [v for _, v in items], color="#4c72b0")
        plt.xlabel("Feature importance")
        plt.title(f"What explains AFP type? (RandomForest, LOO acc={acc})")
        plt.tight_layout(); plt.savefig(f"{FIG_DIR}/model_feature_importance.png", dpi=150); plt.close()
        print(f"Saved {FIG_DIR}/model_feature_importance.png")


if __name__ == "__main__":
    main()
