# RT(maze routing) 검증 결과 (2026-09-20)

`core/router.py`(RT, maze routing) 구현 후 6칩 전체에 대해 측정한 결과와, 그 과정에서
나온 설계 결정 근거를 정리한다. `docs/20260917_gp_qc_findings.md`와 같은 형식 — 전부
측정 기반이고, 재현 커맨드를 하단에 남긴다.

## 목차
- 배경 — RT의 범위와 GP/DP와의 인터페이스
- 발견 1 — 라우팅 순서(rt_coupler_order)는 결과에 전혀 영향이 없다 (구조적 이유)
- 발견 2 — 꺾임 페널티(rt_turn_penalty)는 이 6칩 범위에서 효과가 작다
- 발견 3 — 리드 leg를 미앤더 후보에서 빼자 교차가 6칩 중 5칩에서 0이 됐다
- 발견 4 — 길이 미달(len_err)이 크고 구조적이다 (DP 압축과의 상호작용)
- 최종 검증표 (6칩)
- 재현

## 배경 — RT의 범위와 GP/DP와의 인터페이스

RT는 GP(`core/globalplacement.py`)가 이미 확정한 `Coupler.segments`(어느 셀이 이
커플러 차지인가) 안에서, 실제로 선을 어떻게 그을지(꺾임 최소화 + 목표 길이 채우기)만
정한다. GP의 `segment_occupied`가 이미 전 칩에 걸쳐 셀을 커플러당 정확히 하나씩만
배정했으므로, "자기 세그먼트 셀만 통과 가능"이라는 필터 하나로 "다른 배선과 비교차"가
구조적으로 따라온다 — 리드인/아웃(포트↔첫/마지막 셀 직선 구간)은 예외이고, 발견 3이
그 예외를 다룬다.

DP(`core/detailedplacement.py`)가 전 좌표를 alpha배 닮음변환하므로 GP가 만든 전역
원점(0,0) 기준 격자는 RT 시점엔 무효하다 — `_build_cell_index()`가 커플러별로 자기
세그먼트의 최소 x/y를 로컬 원점 삼아 격자 인덱스를 역산한다(같은 커플러의 세그먼트끼리는
닮음변환 후에도 정확히 `segment_size_um`(post-DP) 간격으로 정렬돼 있으므로 안전).

## 발견 1 — 라우팅 순서(rt_coupler_order)는 결과에 전혀 영향이 없다 (구조적 이유)

6칩 × 3전략(shortest_first/longest_first/random) 전수 비교 — routed/failed 수,
교차 수, turns_avg/max, 길이 오차, 이탈 길이가 **모든 지표에서 소수점까지 완전히
동일**했다(부동소수 잔차조차 없음 — 재현 커맨드 참고).

이유: GP는 `segment_occupied`라는 커플러 루프 전체가 공유하는 가변 상태가 있어서
순서가 "누가 경합 셀을 먼저 차지하는가"를 바꾼다. RT의 `_route_coupler`는 그런 공유
상태가 전혀 없다 — 자기 커플러의 (이미 GP가 확정한) `Coupler.segments`와 자기
`port_assignment` 항목만 읽는다. 그래서 "어느 커플러가 먼저 라우팅되는가"는 다른
어떤 커플러의 결과에도 영향을 줄 수 없다 — 순서 무관은 튜닝 결과가 아니라 GP/RT
분리 구조 자체에서 나오는 성질이다. 01_mainref의 단일 패스 조인트 라우터(넷 사이에
점유/충돌 상태를 실시간 공유하는 라우터)에서 shortest-first가 실제로 교차를 줄였던
것(2/624 -> 0/624)과는 다른 상황 — 그 라우터는 이 저장소처럼 GP가 미리 셀 소유권을
분리해두지 않는다.

`gp_coupler_order`의 값 그대로 `rt_coupler_order=shortest_first`를 기본값으로 뒀다 —
`_coupler_order`(core/globalplacement.py) 재사용이라 비용이 0이고, 향후 RT가 리드
leg 충돌 회피처럼 커플러 간 공유 상태를 갖는 쪽으로 재설계되면(발견 3의 남은 과제)
다시 의미를 가질 수 있어 인터페이스를 유지한다.

## 발견 2 — 꺾임 페널티(rt_turn_penalty)는 이 6칩 범위에서 효과가 작다

grid_25/aspen_11에 대해 turn_penalty ∈ {0, 1, 3, 5, 10, 20}을 스윕했다.

- **grid_25**: turns_avg/turns_max/len_err 전부 0~20 구간에서 완전히 동일. 이 칩의
  커플러당 own-cell 블록이 평균 ~9.6칸으로 작아서(각주: `num_segments`가 GP 채우기
  당시 값이며 DP 압축 후에도 개수는 불변) A*가 애초에 고를 수 있는 대안 경로가 없다 —
  꺾임 비용을 아무리 올려도 바꿀 경로 자체가 없다.
- **aspen_11**: 0->1에서만 작은 개선(turns_avg 20.46->20.28, turns_max 61->59), 1
  이후 20까지는 평탄. penalty>=1이면 이 벤치마크 풀에서 얻을 수 있는 꺾임 감소 여지는
  이미 거의 다 얻는다는 뜻.

기본값 3.0을 선택했다: 측정상 효과가 포화되는 지점(1)보다 확실히 위라 이 6칩 풀에
없는 더 넓은 회랑(대안 경로가 실제로 있는 경우)에도 여유가 있고, 그러면서도 셀 이동
기본비용(1.0)의 ~4배 수준으로 낮게 유지해 "꺾임 하나 피하자고 훨씬 긴 직진 우회를
택해 긴 leg를 다 써버리는"(미앤더가 접을 leg 자체가 줄어드는) 부작용을 피한다.

## 발견 3 — 리드 leg를 미앤더 후보에서 빼자 교차가 6칩 중 5칩에서 0이 됐다

**최초 구현(리드 leg 포함)**: 미앤더를 "가장 긴 leg 하나"에만 접다가 → "모든 leg에
욕심 있게 분배"로 바꾼 뒤 측정했더니 xtree_53 22건, eagle 2건, aspen_m 2건 교차가
나왔다. `core/router.py`의 `find_crossing_legs()`로 교차 쌍을 leg 인덱스까지 추적한
결과, 전부 **포트→첫 셀(리드인) 또는 마지막 셀→포트(리드아웃) leg가 미앤더로 접히면서**
생긴 것이었다.

원인: `_select_endpoint_cell()`은 포트를 포함하는 자기 셀이 있으면 그걸 쓰지만
(리드 길이 0), 없으면 가장 가까운 자기 셀로 대체한다 — 이 경우 리드 직선은
own_cells 안에 있다는 보장이 전혀 없다(포트에서 코리도까지 임의 방향/거리일 수
있음). 미앤더 진폭(`_MEANDER_AMPLITUDE_FACTOR * cell`)은 "leg가 자기 셀들로 이뤄진
한 줄(row/column)을 따라간다"는 전제에 기대는데, 리드 leg는 이 전제가 성립하지
않는다 — 그 상태로 접었더니 실제로 이웃 커플러 영역을 침범했다.

**수정**: `_route_coupler`가 `foldable` 마스크(첫/마지막 leg=False, 나머지 내부
leg만 True)를 계산해 `_apply_meander`에 넘긴다 — 리드 leg는 항상 미접힌 채로 남는다.
수정 후 재측정:

| 칩 | 교차(수정 전) | 교차(수정 후) |
|---|---|---|
| grid_25 | 0 | 0 |
| xtree_53 | 22 | 6 |
| falcon | 0 | 0 |
| eagle | 2 | 0 |
| aspen_11 | 0 | 0 |
| aspen_m | 2 | 0 |
| **합계** | **26** | **6** |

**xtree_53에 남은 6건**은 리드 leg 자체(미접힌 원래 직선)와 다른 커플러의 (정상적으로
접힌) 내부 leg 사이 교차다 — `find_crossing_legs`로 확인: 전부 `lead-interior` 유형,
`interior-interior`는 0건. 이건 미앤더 접기의 버그가 아니라 **더 근본적인 GP/RT
인터페이스 한계**다: GP의 boustrophedon 채우기는 박스의 p1쪽 모서리부터 필요한
개수(`num_segments`)만 채우고 멈춘다 — 박스가 두 포트를 포함할 만큼 큰데 필요
세그먼트 수가 박스 용량보다 훨씬 적으면, 세그먼트 뭉치가 p1 근처에만 몰리고 p2
근처엔 자기 셀이 전혀 없을 수 있다(실측: aspen_11 최악 사례에서 포트→가장 가까운
자기 셀 거리가 p1쪽 77um인데 p2쪽은 1667um). 이 경우 p2쪽 리드는 필연적으로 own_cells
밖을 길게 지나가고, 그 직선이 다른 커플러의 (정상 배선된) 영역과 마주치면 교차가
생긴다. RT가 이걸 완전히 없애려면 "다른 모든 커플러의 소유 셀을 피해 리드 경로 자체를
다시 미로 탐색"해야 하는데, 이는 커플러 간 공유 상태가 필요한 **조인트 문제**라
(발견 1이 설명하는 RT의 "커플러 독립" 설계와 정면으로 충돌한다) 이번 범위에서는
손대지 않았다 — 5절 원칙(막힌 걸 무시하는 fallback을 만들지 않는다)을 지키면서 이
문제까지 풀려면 GP 쪽에서 "박스 안에서 두 포트 모두에 가깝게 채운다"처럼 채우기
전략을 바꾸는 게 맞는 방향으로 보인다(별도 후속 과제로 남긴다). 지금은 **정직하게
측정해 보고**한다 — `route_report()`의 `crossings`/`oob_length_um`이 그 측정값이다.

## 발견 4 — 길이 미달(len_err)이 크고 구조적이다 (DP 압축과의 상호작용)

6칩 전체에서 `route_length - l`(달성 길이 - 목표 반파장 길이)의 평균이 -6600~-8800um
로 전부 큰 음수다 — 목표 길이에 크게 못 미친다. 원인은 RT의 버그가 아니라 상류
설계의 이미 문서화된 간극이 DP를 거치며 더 벌어진 것이다:

- GP가 커플러당 셀 개수(`num_segments`)를 정할 때 쓰는 면적 예산은
  `wire_area = l * meander_spacing_um`(40um)이고, 이걸 (pre-DP) `segment_size_um`
  (200um) 격자 칸 수로 나눈다.
- DP(`core/detailedplacement.py`)는 큐빗+세그먼트 전 좌표를 alpha배 축소하면서
  **세그먼트 "크기"도 같은 alpha배로 줄인다**(모듈 docstring: "segment_size_um 자체가
  GP의 배치 이산화 단위일 뿐 고정된 물리 치수가 아니다"). 문제는 `meander_spacing_um`
  (실제 CPW 폴드 피치, 물리량)은 alpha와 무관하게 고정이라는 것 — 그래서 압축된
  커플러 영역의 실제 면적은 alpha² 배로 줄어드는데, 그 안에 필요한 배선 길이(l)는
  그대로다. grid_25는 alpha_min=0.4739였으므로 이론상 최대로 채워도(면적을 100%
  meander로 빈틈없이 쓴다고 가정해도) 목표 길이의 alpha² ≈ 22.5%밖에 못 채운다 —
  이건 DP 모듈 docstring이 이미 "이 간극은 파이프라인의 segment 모델 자체가 처음부터
  갖고 있던 근사이고, DP가 새로 만든 문제가 아니다"라고 명시한 부분이다.

RT가 할 수 있는 건 "주어진 회랑을 최대한 효율적으로 쓰는 것"뿐이다 — 미앤더를 leg
하나에 몰아넣는 대신 긴 leg부터 순서대로 각 leg의 물리적 한도(진폭 `0.45*cell`,
최소 피치 `min(meander_spacing_um, amplitude)`)까지 채우는 탐욕적 분배로 바꾼 뒤
(`_apply_meander`) 그나마 나아졌지만(예: grid_25 turns_avg가 초기 단일-leg 버전
대비 안정적으로 유지되면서 oob_length_um도 줄었다), 근본적인 면적 부족 자체는
RT 범위 밖의 문제라 해소되지 않는다. **이 오차는 버그가 아니라 측정치다** — 목표
길이(l)에 물리적으로 맞는 배선이 필요하다면 DP의 alpha 하한을 조정하거나
GP/DP의 면적 예산 모델을 재검토해야 한다(후속 과제).

## 최종 검증표 (6칩, rt_coupler_order 무관 — 발견 1)

| 칩 | 커플러 | 세그먼트있음 | 라우팅성공 | 실패 | 교차 | 꺾임avg/max | 길이오차avg(um) | 길이오차 max\|abs\|(um) | 이탈길이(um) |
|---|---|---|---|---|---|---|---|---|---|
| grid_25 | 40 | 40 | 40 | 0 | 0 | 3.18 / 12 | -8845.9 | 9660.1 | 1669.8 |
| xtree_53 | 52 | 49 | 49 | 0 | 6 | 17.73 / 28 | -6648.0 | 8305.4 | 8977.0 |
| falcon | 28 | 27 | 26 | 1 (no_path) | 0 | 9.92 / 26 | -7977.4 | 9584.4 | 2394.5 |
| eagle | 144 | 140 | 140 | 0 | 0 | 13.19 / 24 | -7333.5 | 9420.9 | 17469.9 |
| aspen_11 | 48 | 46 | 46 | 0 | 0 | 11.26 / 41 | -7200.8 | 9189.9 | 15655.0 |
| aspen_m | 106 | 105 | 105 | 0 | 0 | 13.53 / 22 | -7269.2 | 9367.7 | 23784.1 |
| **합계** | **418** | **407** | **406** | **1** | **6** | | | | |

- **세그먼트있음 = 418-11**: GP가 세그먼트를 못 놓은 11개(`docs/20260917_gp_qc_findings.md`
  발견 4)는 RT가 다룰 대상이 아니라 건너뛴다(`route_failures`에 집계 안 됨 —
  `num_unplaced_couplers`로 이미 보고됨).
- **실패 1건 (falcon (1,2), no_path)**: GP는 성공했지만(9개 세그먼트, `num_segments`
  프로퍼티 기준 13개 필요) own_cells가 4-연결 기준 두 조각(5칸+4칸)으로 쪼개져
  있었다 — 다른 커플러가 중간 열(column 1~2)을 먼저 가져간 결과. GP의 박스-박스
  경합(`gp_coupler_order` 설명 참고)이 원인 — RT는 이 경우를 "막힌 셀을 억지로
  뚫는" fallback 없이 정직하게 실패로 보고한다(5절 원칙).
- 실행 시간: RT 자체는 6칩 합계 20~30ms 수준(A*/미앤더가 커플러당 셀 수십 개 규모라
  거의 즉시 끝남) — FP~DP 전체(6칩, ~8.5초)에 비하면 무시할 수준.
- 렌더링(`output/4_RT/`)은 실제 배선(초록)을 그리고, 교차가 있는 leg만 경고색(빨강,
  굵게)으로 강조한다 — `utils/rendering.py`의 `_draw_coupler_route`.

## 재현

```bash
python3 - <<'PY'
import sys, time
sys.path.insert(0, "src")
from core.state import ChipState
from core.floorplan import Floorplan
from core.globalplacement import GlobalPlacement
from core.legalization import Legalization
from core.detailedplacement import DetailedPlacement
from core.router import Router, route_report
from utils.params import Params
from utils.parsing import Parser

params = Params()
bench = Parser().load_configs(["grid_25", "xtree_53", "falcon", "eagle", "aspen_11", "aspen_m"])
states = [ChipState.processor_config(c, params) for c in bench]
states = Floorplan(params).run(states)
states = GlobalPlacement(params).run(states)
states = Legalization(params).run(states)
states = DetailedPlacement(params).run(states)

for strategy in ("shortest_first", "longest_first", "random"):
    params.rt_coupler_order = strategy
    rt = Router(params)
    out = rt.run(states)
    for s in out:
        print(strategy, route_report(s))
    print(strategy, "failures:", rt.route_failures)
PY
```
