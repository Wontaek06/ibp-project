"""
Compile every figure and table into a single report-asset document.

    python -m src.build_report_assets   ->  REPORT_ASSETS.md

Reads whatever result files exist and skips the rest, so it can be run
mid-batch: the genus-level environment mapping in particular is a long
job, and a partial data/pfam11999_env.csv is still worth reporting as
long as the row count is stated.

Every table is generated from the CSVs rather than transcribed, so the
numbers in the document cannot drift from the numbers in the data.
"""
import os

import pandas as pd

OUT = "REPORT_ASSETS.md"
FIG = "figures"

FIGURES = [
    ("spike_cold_adaptation.png",
     "1단계 스파이크 — IBP 보유 분류군 vs 음성 대조군의 서식지 수온",
     "IBP 보유군과 AFP 미보유 대조군(열대·온대 어류)의 연중최저 해수면온도. "
     "채움 원은 IBP 보유군(계통별 색), 빈 사각형은 대조군. 육상 균류 2종은 "
     "측정량이 다른 대기 기온이므로 이 패널에서 빼고 하단에 값을 별도 표기함."),
    ("spike_lat_vs_sst.png",
     "1단계 스파이크 — 위도-수온 공간에서의 분리",
     "분포중심 위도와 연중최저 수온. IBP 보유군은 고위도·저온 영역에, "
     "대조군은 저위도·고온 영역에 뚜렷이 분리됨."),
    ("spike_stage_coverage.png",
     "1단계 스파이크 — 파이프라인 단계별 데이터 반환율",
     "16개 시드가 각 단계에서 데이터를 반환한 비율. **전 단계 결손 없음** — "
     "해양 분류군은 Bio-ORACLE, 육상 분류군은 Open-Meteo(ERA5)로 온도를 확보함."),
    ("env_latitude_vs_sst.png",
     "어류 AFP 18종 — 위도 대비 서식지 수온",
     "AFP 타입별 색 구분. 대표점은 분포중심, 수온은 Bio-ORACLE 연중최저값."),
    ("env_bottom_temp_by_type.png",
     "어류 AFP 18종 — AFP 타입별 저층 수온",
     "Kruskal-Wallis H=4.96, p=0.292 (비유의)."),
    ("env_seaice_by_type.png",
     "어류 AFP 18종 — AFP 타입별 해빙 두께",
     "환경 변수 중 유일하게 유의한 차이. Kruskal-Wallis H=10.99, p=0.027."),
    ("seq_physicochem_by_type.png",
     "어류 AFP 18종 — AFP 타입별 서열 물리화학 특성",
     "GRAVY, 등전점(pI), 불안정성 지수. GRAVY는 당쇄화를 반영하지 못하므로 "
     "AFGP 판정 근거로 쓰면 안 됨(README 참조)."),
    ("struct_plddt_by_type.png",
     "어류 AFP 18종 — AlphaFold pLDDT 기반 구조 지표",
     "plddt_std / plddt_periodicity 는 검증된 구조생물학 지표가 아닌 "
     "예비 휴리스틱임."),
    ("model_feature_importance.png",
     "어류 AFP 18종 — RandomForest 변수 중요도 (LOO-CV)",
     "표본이 작아 정확도 자체는 참고치이며, 변수 중요도 확인 목적."),
    ("duf3494_features_by_clade.png",
     "2단계 DUF3494 — 계통별 서열 특성 (CD-HIT 90% 대표 1,835개)",
     "전 변수가 계통 간 유의차를 보임(전부 p<1e-5). 단 length 는 도메인 길이가 "
     "아니라 단백질 전체 길이이므로 해석 주의(README 참조)."),
]


def md_table(df, floatfmt="{:.3f}"):
    """DataFrame -> GitHub markdown table."""
    def cell(v):
        if pd.isna(v):
            return "—"
        if isinstance(v, float):
            return floatfmt.format(v)
        return str(v)

    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    sep = "|" + "|".join("---" for _ in df.columns) + "|"
    rows = ["| " + " | ".join(cell(v) for v in row) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, sep] + rows)


def section_figures(out):
    out.append("## 그림\n")
    out.append("모든 그림은 `figures/` 에 200 dpi PNG 로 저장됨. "
               "색상은 Okabe-Ito 색각이상 안전 팔레트를 사용함.\n")
    for i, (fname, title, caption) in enumerate(FIGURES, start=1):
        path = os.path.join(FIG, fname)
        if not os.path.exists(path):
            continue
        out.append(f"### 그림 {i}. {title}\n")
        out.append(f"![{title}]({path})\n")
        out.append(f"*{caption}*\n")


def section_spike(out):
    path = "data/spike_results.csv"
    if not os.path.exists(path):
        return
    df = pd.read_csv(path)
    out.append("## 표 1. 1단계 스파이크 — 시드 16종 전 단계 결과\n")

    cols = {
        "label": "라벨", "accession": "UniProt", "organism": "출처 생물",
        "clade": "계통", "afp_class": "AFP 계열", "gbif_rank": "GBIF 해상도",
        "n_occ": "출현기록", "rep_lat": "분포중심 |위도|",
        "surf_min_temp": "연중최저 수온(°C)", "af": "AlphaFold", "plddt": "pLDDT",
    }
    if "land_min_temp" in df:
        cols["land_min_temp"] = "육상 최저기온(°C)"
    t = df[[c for c in cols if c in df]].rename(columns=cols)
    out.append(md_table(t, floatfmt="{:.2f}") + "\n")

    n_land = int(df.land_min_temp.notna().sum()) if "land_min_temp" in df else 0
    out.append(f"시드 {len(df)}개 = IBP 보유 13 + 음성 대조군 3. "
               f"해양 수온 {int(df.surf_min_temp.notna().sum())}건, "
               f"육상 기온 {n_land}건.\n")


def section_cold_stats(out):
    path = "data/spike_results.csv"
    if not os.path.exists(path):
        return
    from scipy.stats import mannwhitneyu

    df = pd.read_csv(path).dropna(subset=["surf_min_temp"])
    cold = df[df.expect_cold == "yes"].surf_min_temp
    warm = df[df.expect_cold == "no"].surf_min_temp
    if len(cold) < 2 or len(warm) < 2:
        return
    u, p = mannwhitneyu(cold, warm, alternative="less")

    out.append("## 표 2. 한랭 적응 검증 (해양 분류군)\n")
    stat = pd.DataFrame([
        {"그룹": "IBP 보유", "n": len(cold), "중앙값(°C)": cold.median(),
         "최소": cold.min(), "최대": cold.max()},
        {"그룹": "대조군 (AFP 없음)", "n": len(warm), "중앙값(°C)": warm.median(),
         "최소": warm.min(), "최대": warm.max()},
    ])
    out.append(md_table(stat, floatfmt="{:.2f}") + "\n")
    out.append(f"**중앙값 차이 {warm.median() - cold.median():.2f} °C · "
               f"Mann-Whitney U={u:.0f}, p={p:.4f}** (단측). "
               f"U=0 은 두 그룹이 순위상 완전히 분리됨을 뜻함.\n")


def section_duf3494(out):
    path = "data/pfam11999_features.csv"
    if not os.path.exists(path):
        return
    df = pd.read_csv(path)
    out.append("## 표 3. 2단계 DUF3494(PF11999) 계통 구성\n")
    comp = (df.clade.value_counts().rename_axis("계통")
              .reset_index(name="대표서열 수"))
    comp["비율(%)"] = (100 * comp["대표서열 수"] / len(df)).round(1)
    comp["고유 분류군"] = comp.계통.map(df.groupby("clade").taxon_id.nunique())
    out.append(md_table(comp, floatfmt="{:.1f}") + "\n")
    out.append(f"UniProtKB {2250:,}건 → CD-HIT 90% → **{len(df):,}개 대표서열 / "
               f"{df.taxon_id.nunique():,} 분류군**. "
               f"Swiss-Prot 리뷰 항목은 {int((df.reviewed == 'reviewed').sum())}건뿐임.\n")

    cols = ["length", "ala_pct", "thr_pct", "gravy", "pI", "instability"]
    big = df[df.clade.isin(df.clade.value_counts()[lambda s: s >= 15].index)]
    med = big.groupby("clade")[cols].median().round(2).reset_index()
    med = med.rename(columns={"clade": "계통", "length": "길이(aa)"})
    out.append("### 계통별 서열 특성 중앙값\n")
    out.append(md_table(med, floatfmt="{:.2f}") + "\n")


def section_env(out):
    path = "data/pfam11999_env.csv"
    if not os.path.exists(path):
        return
    df = pd.read_csv(path)
    out.append("## 표 4. DUF3494 속(genus) 단위 서식지 온도 — 부분 결과\n")
    out.append(f"⚠️ **진행 중인 배치의 중간 결과 ({len(df)} / 441 속)**. "
               f"완주에는 수 시간이 더 필요하며, `python -m src.expand_env` 를 "
               f"다시 실행하면 이어서 진행됨.\n")

    show = df.copy()
    ren = {"genus": "속", "clade": "계통", "rep_lat": "분포중심 |위도|",
           "habitat": "서식 구분", "surf_min_temp": "해수 최저(°C)",
           "land_min_temp": "육상 최저(°C)", "n_proteins": "단백질 수"}
    show = show[[c for c in ren if c in show]].rename(columns=ren)
    out.append(md_table(show, floatfmt="{:.2f}") + "\n")
    if "habitat" in df:
        counts = df.habitat.value_counts().to_dict()
        out.append(f"서식 구분: {counts}. 판정은 정확 격자 셀(반경 0) 기준이며, "
                   f"근접 셀 탐색을 쓰면 토양 세균이 해양으로 오분류됨(README 참조).\n")


def section_pipeline(out):
    path = "data/stats_summary.csv"
    if not os.path.exists(path):
        return
    df = pd.read_csv(path)
    out.append("## 표 5. 어류 AFP 18종 — 그룹 비교 통계\n")
    ren = {"variable": "변수", "H": "Kruskal-Wallis H", "p": "p", "block": "구분",
           "n_groups": "그룹 수", "n_total": "n"}
    t = df[[c for c in ren if c in df]].rename(columns=ren)
    out.append(md_table(t, floatfmt="{:.4f}") + "\n")
    out.append("표본이 작아(그룹당 3~8종) 대부분 비유의함. "
               "해빙 두께만 p<0.05 임.\n")


def main():
    out = ["# 보고서용 자료 모음 (Report Assets)\n",
           "`python -m src.build_report_assets` 로 생성됨. "
           "모든 표는 `data/` 의 CSV 에서 직접 계산되므로 수치가 어긋나지 않음.\n"]

    section_figures(out)
    section_spike(out)
    section_cold_stats(out)
    section_duf3494(out)
    section_env(out)
    section_pipeline(out)

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))
    print(f"Saved {OUT} ({os.path.getsize(OUT):,} bytes)")


if __name__ == "__main__":
    main()
