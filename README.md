# IBP 프로젝트 — 위도·수온 구배와 항결빙단백질(AFP)

## 실행 환경
- Python 3.10+
- `pip install -r requirements.txt`
- 인터넷 필요 (GBIF, Bio-ORACLE ERDDAP, UniProt, AlphaFold DB API 호출)

## 실행 순서

**Colab에서 실행 시**: 이 폴더 전체를 압축 해제 후 업로드하거나, 각 `src/*.py` 파일 내용을
`%%writefile src/파일명.py` 셀로 만들어 실행. 이후:

```python
import sys; sys.path.insert(0, ".")   # 필요시
from src.pipeline import main
main()
```

**로컬/터미널 실행 시**:
```bash
cd ibp_project
pip install -r requirements.txt
python -m src.pipeline
```

첫 실행은 몇 분 걸림 (API rate limit 대응 sleep 포함). **2회차부터는 캐시 덕분에 훨씬 빠름.**

### 캐싱
모든 API 호출(GBIF, Bio-ORACLE, UniProt, AlphaFold)은 `.cache/` 폴더에 자동 저장됨.
분석 코드만 고치고 재실행할 때 API를 다시 안 때림.

```python
from src.cache import clear_cache
clear_cache()              # 전체 캐시 삭제
clear_cache("bio_oracle")  # 특정 API만 삭제
```
환경변수 `IBP_NO_CACHE=1` 을 주면 캐시를 완전히 우회함.

**주의**: 캐시 키는 함수 *인자*만 보고 만들어짐. fetch 함수의 *반환 형식*을 바꿨다면
반드시 `clear_cache()` 할 것 (안 그러면 옛 형식을 계속 읽음).

## 데이터 흐름

```
data/seed_species.csv (18종, afp_type + validation_status)
        │
        ├─ src/fetch_gbif.py        → 종명 → GBIF 관측좌표
        ├─ src/occurrence.py        → 좌표 → 분포중심 대표점 (core_points)
        ├─ src/fetch_bio_oracle.py  → 좌표 → 표층최저·저층수온·해빙두께
        ├─ src/fetch_uniprot.py     → 종 → 검증된 AFP 서열 (reviewed:true만)
        ├─ src/fetch_alphafold.py   → accession → 구조 신뢰도 + pLDDT 프로파일
        ├─ src/features.py          → 서열 물리화학 특징 (GRAVY, pI 등)
        └─ src/stats_model.py       → 그룹비교 통계 + 분류모델(feature importance)
                │
                ▼
        data/final_merged.csv, figures/*.png
```

## 1단계 스파이크 — 미생물·조류 IBP 포함 전 구간 관통 (2026-07)

`python -m src.spike` → `data/spike_results.csv`

기존 파이프라인은 **어류 AFP 전용**이었음. 본선의 DUF3494(미생물·균류·규조류) 확장이
가능한지 보려고, 문헌상 특성이 규명된 IBP 13개 + **음성 대조군 3종**(열대·온대 어류,
AFP 없음)을 `data/spike_seed.csv` 에 고정하고 전 단계를 한 번 관통시킴.

시드는 키워드 검색이 아니라 **UniProt accession 을 직접 지정**함 — 키워드 검색이
`seed_species.csv` 의 EXCLUDE 5종을 만든 원인이기 때문.

### 단계별 결과 (16개 시드)

| 단계 | 통과 | 비고 |
|---|---|---|
| UniProt 서열 | 13/13 | ColAFP, FfIBP, MpAFP, LeIBP, TisAFP, FcIBP1, CnAFP + 어류 6 |
| 출처 생물명 | 16/16 | |
| GBIF 분류 매칭 | 16/16 | 전부 EXACT |
| 출현 좌표 | 16/16 | 전부 GBIF, OBIS 폴백 불필요 |
| Bio-ORACLE 수온 | **14/16** | 결손 2건 모두 **육상 균류** (아래 참조) |
| AlphaFold 구조 | 13/13 | pLDDT 75.1~95.4, 전부 모델 존재 |

### 한랭 적응 검증

| 그룹 | n | 표층 연중최저수온(중앙값) | 분포중심 위도 |
|---|---|---|---|
| IBP 보유 | 11 | **-0.25 °C** (범위 -1.99 ~ 4.41) | 63.5° |
| 음성 대조군 | 3 | **14.04 °C** (범위 12.38 ~ 23.20) | 28.0° |

**분리 14.30 °C, Mann-Whitney U=0.0, p=0.0027 (완전 분리)** — "추운 생물이 실제로
저수온 좌표에 매핑되는가"에 대한 답은 예. 프로젝트 전제가 성립함.

### 스파이크가 잡아낸 결함 (기존 코드 포함)

1. **`protein_name:antifreeze` 쿼리로는 미생물 IBP 가 전멸함.** UniProt 등재명이
   "Ice-binding protein" 이라 한 건도 안 잡힘. 미생물까지 보려면 `protein_name:"ice-binding"`
   을 OR 로 넣어야 함.
2. **`reviewed:true` 필터는 남극 규조류를 통째로 날림.** *Fragilariopsis cylindrus*
   IBP1~6(D0FHA3~A8), *Chaetoceros neogracilis*(D2DLE1)는 전부 TrEMBL 임.
   어류만 볼 땐 안전한 필터였지만 확장 시엔 못 씀.
3. **GBIF `species/match` 에 이름만 넘기면 세균 속이 `matchType: NONE` 이 됨.**
   ("Colwellia" → Multiple equal matches) UniProt lineage 에서 kingdom 을 뽑아
   힌트로 주면 해결됨 → `src/taxonomy.py`.
4. **⚠️ 겨울 관측월 필터가 극지 taxa 를 저위도로 끌어내리고 있었음** (기존
   `fetch_gbif.cold_points`). 극지 현장조사는 여름에 집중되기 때문 —
   Colwellia 는 관측 300건 중 231건이 3월이고 12~2월은 7건뿐인데, 그 7건이
   56.6°N·41.3°S 라 실제 중앙값 78.7°N 인 taxon 이 3.5 °C 로 찍혔음.
   Leucosporidium 은 남위 42도·11.6 °C 로 찍혔음. **필터 제거함** — 계절 극값은
   Bio-ORACLE `thetao_min`(연중 최저)이 이미 담당하므로 이중 계산이었고,
   샘플링 노력 편향만 들여왔음.
5. **대표점을 "가장 극쪽 k개"로 잡으면 안 됨.** 분포 꼬리는 미아(vagrant) 기록과
   샘플링 편향이 지배함 — 황다랑어(분포중심 7.7°N)가 2.7 °C, 홍돔이 1.98 °C 로
   나왔음. **분포 중심(`core_points`)으로 전환**함. 꼬리는 `tail_lat` 컬럼에
   위도만 보존(서식 한계 질문용).
6. **좌표 중복 제거 필요.** 한 정점 반복 조사가 대표점 5개 중 4개를 차지하는 일이
   있었음. 단, **중복 제거는 대표점 선택에만 적용하고 분포중심 계산에는 적용하면 안 됨**
   — 관측 밀도가 사라져 Colwellia 중심이 78.7°N→56.6°N 으로 밀림.

→ 4·5·6 은 `fetch_gbif.py` 에도 있던 문제라 **기존 18종 결과도 같은 편향을 받음.
   `python -m src.pipeline` 재실행 필요.**

### 확인된 한계: 육상 IBP 는 해양 레이어로 못 잡음

`LeIBP`(*Leucosporidium* sp. AY30, 북극 효모)와 `TisAFP`(*Typhula ishikariensis*,
눈곰팡이)는 좌표는 정상적으로 나오지만(61.1°N / 47.8°N) **Bio-ORACLE 수온이 null** 임.
육상 종이라 해양 격자에 없음. 버그가 아니라 데이터 도메인 불일치.
→ 육상 taxa 는 WorldClim/CHELSA 기온을 붙이거나, 해양 taxa 로 범위를 한정해야 함.

또한 **균주 수준 서식지는 속 단위 GBIF 좌표로 대리할 수 없음** — *Leucosporidium* 속은
핀란드 북방림 등 전 세계 분포라, 북극 균주 AY30 의 실제 채집 환경과 다름.
(인수인계 5번에 적힌 우려가 실측으로 확인된 셈)

## 2단계 — DUF3494 확장 경로 (검증 완료)

`python -m src.fetch_pfam` → `data/pfam11999_proteins.tsv` / `.fasta`

- **Pfam PF11999 (DUF3494, "Ice-binding-like")**: UniProtKB 2,250건, InterPro 기준
  2,322건, AlphaFold 모델 1,582개 존재.
- UniProt 검색 필드는 **`xref:pfam-PF11999`** 임. `xref_pfam:PF11999` 은
  "not a valid search field" 로 거부되고, `database:pfam` 은 Pfam 상호참조가 있는
  전체 엔트리(~1.1억)를 반환하므로 쓰면 안 됨.
- InterPro API 의 `?format=fasta` 는 404 — 서열은 UniProt `/stream` 으로 받을 것.
- CD-HIT 는 미설치 상태: `brew install cd-hit` 또는 `conda install -c bioconda cd-hit`

## 알아둬야 할 것 (데이터 정제 이력)

1. **환경 변수는 연평균이 아니라 "분포중심 관측점의 연중 최저수온"** 을 씀
   (`thetao_min`). 연평균은 온대종의 겨울 저온 노출을 가려서 AFP 유도 조건을 못 잡음.
   계절 극값은 이 레이어가 담당하므로 **관측월로 좌표를 거르지 말 것** —
   그렇게 하면 극지 taxa 가 오히려 저위도로 끌려감(1단계 스파이크 4번 항목).
2. **Bio-ORACLE 좌표 요청은 반경 탐색 포함**: 정확한 좌표가 육지/마스킹 격자에
   걸리면 모든 레이어가 동시에 null이 됨(확인된 사례: Pagothenia borchgrevinki,
   Gadus morhua, Hemitripterus americanus, Osmerus mordax). 주변 0.5도 박스에서
   가장 가까운 유효 셀을 찾도록 처리함.
3. **UniProt 검색은 `protein_name:antifreeze AND reviewed:true`로 한정**.
   느슨한 전체텍스트 검색은 무관한 단백질(Elongation Factor-1a, Cytochrome c
   oxidase I 등)을 "antifreeze"라는 단어가 우연히 포함된 이유로 잘못 매칭함.
   **단, 이 쿼리는 어류 전용임** — 미생물 IBP 는 등재명이 "Ice-binding protein"
   이라 안 잡히고, 규조류 IBP 는 TrEMBL 이라 `reviewed:true` 에서 탈락함.
   확장 시에는 accession 고정(`uniprot_by_accession`) 또는 Pfam 도메인 기반
   수집(`src/fetch_pfam.py`)을 쓸 것.
4. **`data/seed_species.csv`의 `validation_status` 컬럼이 최종 판정**:
   - `verified`: 서열이 해당 AFP 타입의 알려진 특징과 일치함을 확인
   - `RELABEL`: Gadus morhua는 AFGP가 아니라 Type IV(ice-structuring protein)로 확인, 재분류
   - `EXCLUDE`: UniProt에 리뷰된 항목이 없거나(Liopsetta putnami), 검색이
     무관한 단백질을 반환함(Trematomus bernacchii, Boreogadus saida,
     Arctogadus glacialis, Microgadus tomcod). 이 5종은 서열/구조 분석에서
     제외하되, **환경 데이터(GBIF/Bio-ORACLE)는 계속 18종 전체로 사용** —
     종의 존재·서식지 자체는 문제없기 때문.
5. **구조 규칙성 지표(`plddt_std`, `plddt_periodicity`)는 휴리스틱**이며
   검증된 구조생물학 지표가 아님. 보고서에서 "예비적 지표"로 명시할 것.
6. **RandomForest 분류모델은 예측력 자랑용이 아니라 feature importance 확인용**
   — 표본이 작아(n=13 내외) 정확도 자체는 참고치일 뿐임.

## ⚠️ seed_species.csv 의 출처와 한계 (반드시 읽을 것)

`data/seed_species.csv` 는 **코드가 생성하는 파일이 아니라 사람이 작성한 입력 파일**임.
나머지 두 CSV(`final_merged.csv`, `stats_summary.csv`)만 파이프라인이 생성함.

이 파일의 `afp_type` 라벨(어느 종이 AFGP인지, Type I인지 등)은 초안 작성 시
**문헌 기억에 의존해 작성된 것**이며, **모든 그룹 비교 분석이 이 라벨을 기준으로 함**.
라벨이 틀리면 그림·통계가 전부 틀림.

실제로 이미 오류가 하나 발견됨: **Gadus morhua 를 AFGP 로 적었으나 UniProt 검증 결과
Type IV(ice-structuring protein)** 였음 → `RELABEL` 처리.

### ✅ 라벨 검증 완료 (2026-07)
13종 전부 UniProt 등재명과 교차확인 완료. 전부 일치함.
- UniProt은 "antifreeze protein" 대신 **"ice-structuring protein(ISP)"** 명명을 쓰는 경우가 많음 (동의어)
- 이름에 `glyco` 포함 → AFGP / `Type-N` 명시 → 해당 타입 / 번호 없으면 서열로 판정
  (Ala 40%↑=Type I, 130aa↑&Cys 다수=Type II, 100aa↓=Type III)
- Gadus morhua 는 UniProt 이 "Type-4 ice-structuring protein" 으로 명시 → Type IV 재분류 확정

### ⚠️ 검증 중 발견된 추가 한계: 성숙 펩타이드 vs 전구체 혼재
AFGP 그룹의 서열 길이가 종마다 크게 다름:
Chaenocephalus 17aa / Pagothenia 31aa / Dissostichus 33aa / **Notothenia 790aa**

이는 오류가 아니라 AFGP 생물학 때문임 — AFGP 는 큰 **폴리단백질 전구체**로 합성된 뒤
절단되어 다수의 짧은 AFGP 분자가 됨. UniProt 등록이 종에 따라 성숙 펩타이드 또는
전구체 전체로 되어 있음.

영향: 길이 의존적 지표(**pI, instability index**)는 서로 다른 대상을 비교하게 됨.
조성 비율 지표(**ala_pct, thr_pct**)는 상대적으로 견고함(790aa 에서도 Ala 49.7%).
또한 Chaenocephalus(17aa)는 ProtParam 통계 자체가 노이즈에 가까움.

→ 보고서 한계에 명시할 것. 개선하려면 UniProt `features` 의 CHAIN/PEPTIDE 위치로
   전구체에서 성숙 영역만 잘라내 재계산 (팀원 개선 항목 후보).

## 팀원 인수인계 — 발전시킬 만한 지점

현 상태로도 실행되며 결과가 나옴. 아래는 개선 여지:

1. **표본 확대** — 현재 환경 n=18 / 서열·구조 n=13. UniProt 에서 AFP 보유 어류를
   추가 확보하거나, `EXCLUDE` 5종을 NCBI/문헌에서 정확한 accession 으로 복구하면
   통계 검정력이 올라감. (현재 그룹당 2~4종이라 CI 가 매우 넓음)
2. **구조 지표 고도화** — `plddt_std`/`plddt_periodicity` 는 임시 휴리스틱임.
   실제 2차구조 비율(DSSP), 회전반경, 표면 소수성 패치, 추정 ice-binding face
   평탄도 등이 훨씬 의미 있음.
3. **환경 변수 확장** — 현재 수온·해빙 3종. 염분, 용존산소, 수심, 계절 진폭 추가 가능.
4. **모델링** — 현재 RandomForest + LOO-CV(표본이 작아 정확도는 참고치).
   계통(phylogeny) 통제 교차검증(held-out-by-clade)을 넣으면 "환경 신호 vs 계통 암기"
   분리가 가능해짐 — 이게 심사에서 가장 강한 rigor 포인트가 될 수 있음.
5. **미생물 IBP 확장(본선용)** — Pfam DUF3494 도메인 서열(수천 개) + Ocean Gene Atlas
   메타게놈으로 규모 확대. 단, 미생물은 GBIF 좌표가 서식지 대리값으로 부정확하므로
   균주 채집 좌표(NCBI BioSample) 또는 메타게놈 샘플 환경값을 써야 함.

## LLM 사용 내역 (규정상 명시 필요)
분석 설계, 데이터 정제 로직, 코드 작성 및 디버깅 과정에서 Claude(Anthropic)를
보조 도구로 활용함. 생물학적 판단(종 선정, 검증 여부 확정) 및 결과 해석은
팀원이 직접 수행함.
