# IBP 보유 분류군의 서식 환경 모델링
### GBIF 좌표 편향성 분석 및 PAMC 기반 교정

2026 극지 빅데이터-인공지능 활용 경진대회 · 데이터 분석 부문 · 팀 ColdFold

본 저장소는 예선 제출 분석보고서에 수록된 모든 수치·표·그림을 재현하는 코드와 결과 데이터를 포함한다.

---

## 1. 제출 요건 확인

| 항목 | 내용 |
|---|---|
| 사용 언어 | Python 3.10 이상 |
| 코드 인코딩 | UTF-8 |
| 개발 환경(OS) | Ubuntu 22.04 LTS |
| 실행 확인 환경 | Google Colab (Python 3.12) 및 로컬 |
| 재현성 | `python -m src.pipeline` 단일 명령으로 보고서 수록 결과 전체 재생성 |
| 외부 데이터 | 전량 공개·무상, API 키 불요 (§3 참조) |

### 주요 라이브러리 버전

```
requests        2.31.0
pandas          2.0.3
numpy           1.24.4
scipy           1.11.4
scikit-learn    1.3.2
scikit-posthocs 0.7.0
biopython       1.81
matplotlib      3.7.5
```

전체 의존성 목록은 `requirements.txt`, 실행 환경의 설치 버전 전량은 `environment.txt` 에 기록하였다.

---

## 2. AI 및 사전학습 모델 사용 내역

대회 규칙 3항(AI·사전학습 모델 사용)에 따라 아래와 같이 명시한다.

### 사전학습 모델

| 모델 | 제공 | 활용 범위 | 라이선스 |
|---|---|---|---|
| AlphaFold Protein Structure Database | EMBL-EBI / Google DeepMind | 시드 단백질의 예측 구조를 내려받아 잔기별 pLDDT 프로파일 및 파생 지표(평균 pLDDT, 프로파일 표준편차, 자기상관 기반 규칙성) 산출 | CC-BY-4.0 |

구조 예측을 직접 수행하지 않았으며, 공개된 예측 결과를 조회하여 사용하였다.

### 머신러닝 모델

RandomForest 분류기(scikit-learn)를 Leave-One-Out 교차검증과 함께 사용하였다.
표본 규모가 작아 예측 성능을 주장하지 않으며, 어떤 변수가 AFP 타입 구분에
기여하는지 확인하는 용도로만 활용하였다.

### 생성형 AI

분석 설계 검토, 데이터 정제 로직 수립, 코드 작성 및 디버깅, 보고서 문안 정리에
Claude(Anthropic)를 보조 도구로 사용하였다. 대상 종 선정, 라벨 검증, 결과 해석 및
보고서의 판단은 참가자가 직접 수행하였다.

---

## 3. 활용 데이터 출처 및 사용 범위

대회 규칙 2항(데이터 활용)에 따라 명시한다. 전 데이터는 공개·무상 접근이며 API 키를 요구하지 않는다.

### 극지연구소(KPDC) 데이터 — 주 분석 대상

| 자원 | 사용 내용 | 접근 |
|---|---|---|
| PAMC 극지미생물자원은행 | 극지·고산 유래 균주의 채집지·서식지·좌표 정보. 전지구 출현 좌표의 타당성 검증 및 환경 매핑 기준 좌표로 사용 (576 균주 / 479 좌표 / 95 채집 지점) | 공개 검색 |

### 외부 공개 데이터 — 참조 및 검증용

| 데이터 | 제공 | 사용 내용 | 접근 경로 |
|---|---|---|---|
| GBIF Occurrence | GBIF | 종·속별 출현 좌표, 분류 매칭 | REST API |
| Bio-ORACLE v3 | ERDDAP | 표층 연중최저수온, 저층수온, 해빙 두께 | griddap |
| ERA5 재분석 | Copernicus / Open-Meteo | 육상 2 m 일최저기온 | REST API |
| UniProtKB | UniProt Consortium | IBP/AFP 서열, 등재명, 계통 | REST / stream |
| Pfam PF11999 (DUF3494) | InterPro / UniProt | IBP-like 도메인 보유 엔트리 | `xref:pfam-PF11999` |
| AlphaFold DB | EMBL-EBI / DeepMind | 예측 구조 pLDDT | REST API |

접근 일자 및 질의 조건은 보고서 7장에 상세히 기재하였다.

---

## 4. 재현 방법

### 실행

```bash
git clone https://github.com/Wontaek06/ibp-project.git
cd ibp-project
pip install -r requirements.txt
mkdir -p figures report
python -m src.pipeline
```

Google Colab에서 실행하는 경우:

```python
!git clone https://github.com/Wontaek06/ibp-project.git
%cd ibp-project
!pip install -q -r requirements.txt
!mkdir -p figures report

from src.pipeline import main
main()
```

### 산출물

실행이 완료되면 다음이 생성된다.

- `data/` — 통합 결과 및 통계 요약 CSV
- `figures/` — 보고서 수록 그림

### 캐싱

외부 API 응답은 `.cache/` 에 저장되어 재실행 시 재호출하지 않는다.
첫 실행은 API 호출로 수 분이 소요되며, 이후에는 수 초 내에 완료된다.

```python
from src.cache import clear_cache
clear_cache()               # 전체 초기화
clear_cache("bio_oracle")   # 특정 원천만 초기화
```

환경변수 `IBP_NO_CACHE=1` 을 지정하면 캐시를 사용하지 않는다.
수집 함수의 반환 구조를 변경한 경우에는 캐시를 초기화한 뒤 실행한다.

---

## 5. 저장소 구성

```
ibp-project/
├── README.md
├── requirements.txt
├── data/
│   ├── seed_species.csv          # 시드 목록 (사람이 작성한 입력)
│   └── (그 외)                   # 파이프라인이 생성하는 결과 CSV
├── src/
│   ├── pipeline.py               # 전체 실행 진입점
│   ├── fetch_uniprot.py          # 서열 수집
│   ├── fetch_pfam.py             # DUF3494 계열 수집
│   ├── taxonomy.py               # 계통 정보 → GBIF 분류 매칭 보정
│   ├── fetch_gbif.py             # 출현 좌표 수집
│   ├── occurrence.py             # 분포중심 기반 대표점 선정
│   ├── fetch_bio_oracle.py       # 해양 환경 매핑
│   ├── fetch_land_climate.py     # 육상 기온 매핑
│   ├── fetch_alphafold.py        # 예측 구조 및 pLDDT
│   ├── features.py               # 서열 물리화학 특성
│   ├── stats_model.py            # 통계 검정 및 변수 중요도
│   ├── kpdc_crosscheck.py        # PAMC 교차검증
│   └── cache.py                  # API 응답 캐시
└── figures/
```

---

## 6. 분석 파이프라인

```
시드 목록 / PF11999 도메인 질의
        │
        ├─ 서열 수집 (UniProt · 수탁번호 고정)
        ├─ 분류 매칭 (계통 힌트 기반)
        ├─ 출현 좌표 수집 (GBIF)
        ├─ 대표점 선정 (분포중심 기준)
        ├─ 환경 매핑 (해양: Bio-ORACLE / 육상: ERA5)
        ├─ 구조 조회 (AlphaFold DB)
        └─ 특성 산출 → 통계 검정 → 시각화
                │
                ▼
        결과 CSV · 그림
```

### 주요 설계 결정

파이프라인의 좌표 처리 규칙은 보고서 3.2~3.3절에 서술한 네 건의 편향 검토 결과를
반영한 것이다.

1. **관측월 필터를 적용하지 않는다.** 계절 극값은 환경 레이어(연중 최저수온)가
   담당하므로 좌표를 월로 거르면 극지 현장조사의 계절 편중이 유입된다.
2. **대표점은 분포 중심으로 산출한다.** 분포 꼬리는 미아 개체 기록에 지배되어
   상시 서식 환경을 대표하지 못한다. 서식 한계 분석을 위해 꼬리 위도는 별도 보존한다.
3. **경도 ±180 구간을 이어붙인다.** 남극 환극분포종은 날짜변경선 부근에 위치하며,
   경도 범위를 자르면 환경값이 전 레이어에서 결손된다.
4. **서식지 판정에는 반경 탐색을 쓰지 않는다.** 반경 탐색을 허용하면 육상 좌표가
   인근 해양 셀의 수온을 받아 서식지가 오분류된다.
5. **해양 수온과 육상 기온을 하나의 변수로 합치지 않는다.** 해수는 결빙점에서
   하한이 막히지만 육상 기온은 그렇지 않다. 결과 테이블에서 `habitat` 컬럼으로 구분한다.

---

## 7. 데이터 및 방법상의 한계

보고서 6.1절에 상세히 기술하였으며, 코드 사용자가 알아야 할 사항은 다음과 같다.

- **`data/seed_species.csv` 는 파이프라인이 생성하는 파일이 아니라 문헌 검토로
  작성한 입력 파일이다.** 모든 그룹 비교가 이 라벨에 의존한다. 13종 전부에 대해
  UniProt 등재명과 교차확인을 완료하였으며, 그 과정에서 *Gadus morhua* 1건을
  AFGP에서 Type IV로 재분류하였다.

- **전지구 출현 좌표는 미생물의 극지성을 반영하지 못한다.** 보고서 4.5절에서
  576 균주 규모로 정량화한 결과이며, 해당 좌표에 기반한 환경값은 거친 스크리닝
  지표로만 제시하고 정량적 주장에는 사용하지 않는다.

- **AFGP는 UniProt 등재 형태가 성숙 펩타이드와 전구체로 혼재한다**(17 aa ~ 790 aa).
  길이 의존 지표(pI, 불안정지수)는 서로 다른 대상을 비교하게 되며, 조성 비율
  지표(Ala·Thr 비율)가 상대적으로 견고하다.

- **pLDDT 기반 규칙성 지표는 검증된 구조생물학 지표가 아닌 예비적 휴리스틱이다.**

- **RandomForest 결과는 변수 중요도 확인 용도이며 예측 성능을 주장하지 않는다.**
  어류 패널은 그룹당 1~8종으로 신뢰구간이 넓다.

---

## 8. 라이선스 및 이용 조건

본 저장소의 코드는 2026 극지 빅데이터-인공지능 활용 경진대회 제출을 위해 작성되었다.
KPDC 데이터는 대회 참가 목적 외의 복제·재배포·상업적 이용을 하지 않는다.
외부 데이터는 각 제공기관의 이용 조건을 따르며, AlphaFold DB 예측 구조는 CC-BY-4.0 하에 사용하였다.
