# 세그먼트 모델 재검토 — QPlacer 원본 대조 (2026-09-20, 코드 수정 없음)

`/home/LabMember/ngchoi/research/ref/80_QPlacer/Qplacer`의 `qplacer_bm/design_format.py`
(`DesignFormator.create_edges_def_data`/`create_pos_matrix`/`create_metal_resonator`)를
직접 읽고, 우리 파이프라인의 세그먼트 개수·면적·배치 가정을 대조했다. 결론부터: **개수
공식(면적 기준) 자체는 QPlacer를 정확히 그대로 옮긴 게 맞다** — `num_poly` 공식의
차원분석이 그걸 수학적으로 강제한다. 다만 **QPlacer는 그 블록들을 "각자 안에서 꽉 채워
접는 통"으로 쓰지 않는다** — 체인으로 이은 앵커 점일 뿐이고, 실제 접기(meander)는 죽은
코드가 보여주는 대로 앵커 사이 경로의 일부(M 구간)에서, 전역 목표 길이에 맞춰 유연하게
일어난다. 우리 RT(v1도 v2도)는 "블록 하나 = 그 안에서 꽉 채워 접는 하드 컨테이너"로
구현했는데, 이건 QPlacer의 실제 메커니즘보다 더 엄격한 근사다. 그리고 완전히 별개로,
**자기 큐빗 겹침은 실측 결과 실제로 심각하다** — 세그먼트가 큐빗 몸체 안쪽으로 최대
182um(큐빗 반폭 200um 대비 91%)까지 파고든다.

## 목차
- 발견 1 — QPlacer의 세그먼트(wireblk)는 배선 한 토막(체인 링크)이지 접기 통이 아니다
- 발견 2 — num_segments 재산정: 길이 기준으로 바꾸면 박스가 칩보다 커진다
- 발견 3 — 자기 큐빗 겹침이 실제로 심각하다 (신규 발견, 수정 필요)
- 발견 4 — 박스-큐빗 겹침: 자기 큐빗은 100% 겹침(정상), 제3자는 실배치 0(기존 결론 유지)
- 결론과 권고

## 발견 1 — QPlacer의 세그먼트(wireblk)는 배선 한 토막(체인 링크)이지 접기 통이 아니다

### 공식 자체는 면적 기준이 맞다

`design_format.py:449`:
```python
num_poly = math.ceil(edge_wirelength * self.padding_size / (self.partition_size**2))
```
차원분석: `edge_wirelength[길이] × padding_size[길이] / partition_size²[길이²]` = 무차원(개수).
뒤집으면 `edge_wirelength ≈ num_poly × partition_size²/padding_size` — **"블록 하나가
담당하는 길이" = `partition_size²/padding_size`라는 뜻이고, 이게 우리 저장소의
`wire_area = l × meander_spacing_um`, `num_segments = ceil(wire_area/segment_size_um²)`
와 정확히 같은 식이다**(대응: `partition_size`↔`segment_size_um`, `padding_size`↔
`meander_spacing_um`). 즉 **우리 공식은 QPlacer를 오독한 게 아니라 정확히 옮긴 것이다.**
(QPlacer 기본값은 `partition_size=0.2mm=200um`(우리와 동일), `padding_size=0.1mm=100um`
— 우리는 01_mainref 전례를 따라 40um을 쓴다, `config/params.json`의 `meander_spacing_um`
설명 참고. 이 차이는 이번 조사 범위 밖.)

### 그런데 블록의 실제 배치는 "접기 통"이 아니라 "체인을 따라 성기게 뿌린 앵커"다

`create_pos_matrix`(`design_format.py:277`)가 `num_poly`개 위치를 어떻게 잡는지 보면:
큐빗 q1-q2를 잇는 방향으로 폭 `rect_width=partition_size`(대략 큐빗 하나 폭, `num_rect_w
= floor(qubit_size/partition_size)`), 길이는 **두 큐빗 사이 거리 전체**에 걸친 직사각형
안에, `num_poly`개 점을 `rect_l × num_rect_w` 격자로 **균등 보간**해서 뿌린다(보스트로페돈
연결 — 이 부분은 우리 GP와 원리가 같다). 즉 블록들은 좁은 한 지점에 빽빽하게 모여
있는 게 아니라 **두 큐빗 사이 전체 거리에 걸쳐 성기게 퍼진 체인**이다. `num_poly`가
작으면(`rect_l==1`) 아예 `interm_pos`로 **p1-p2 직선 위에 균등 간격으로 점을 찍는다**
(`design_format.py:328`) — 이건 "포트에서 포트로 가는 직선을 num_poly등분한 체인"이지
"한쪽 구석에 접기 공간을 몰아준다"가 전혀 아니다.

### 죽은 코드 `create_metal_resonator`가 실제 의도를 보여준다

`design_format.py:83-158`(사용자가 언급한 대로 호출부가 없는 죽은 코드, 하지만 의도는
분명하다). 각 wireblk 중심을 `anchors`로 넘기고, `between_anchors`를 **"PF", "M"을
번갈아** 채운다:
```python
for i in range(len(poly_list)):
    between_anchors[2*i] = "PF"
    between_anchors[2*i+1] = "M"
between_anchors[2*len(poly_list)] = "PF"
```
`assert len(between_anchors) == len(anchors) + 1`. 이건 qiskit-metal `RouteMixed`의
구간 타입 지정이다 — **PF(PathFinder, 직선/충돌회피 경로)와 M(Meander, 접기)이 앵커
사이사이 번갈아 나온다.** 그리고 전체 목표 길이는:
```python
options['total_length'] = f'{len(poly_list)*4}mm'
```
(이 `4mm` 상수 자체는 죽은 코드의 디버그/자리표시자 값으로 보인다 — `partition_size²/
padding_size`=0.4mm과 10배 차이 나고, 살아있는 코드 어디서도 이 상수를 쓰지 않는다.
**상수 값은 신뢰하지 않지만 구조는 신뢰할 수 있다.**) 핵심은 구조: **`total_length`는
경로 전체(PF+M 전부 합)의 목표이고, 실제 접기(M)는 qiskit-metal의 라우터가 그 목표에
맞춰 필요한 만큼만 채운다** — "블록 하나가 반드시 `partition²/padding`만큼을 자기
안에서 자체 완결적으로 접어야 한다"는 요구가 아니다. `advanced: {avoid_collision: true}`,
`snap: true`도 이 경로가 다른 컴포넌트를 피해 유연하게 휘어질 수 있음을 뜻한다 — 블록
경계에 하드하게 갇혀 있지 않다.

### 정리 — 사용자의 가설과 우리 구현, 어느 쪽이 맞나

| | 주장 | 판정 |
|---|---|---|
| "세그먼트 하나 = 배선 한 토막" | 체인 연결(`new_pin_distribution`의 IN/OUT)이 이를 뒷받침 — **맞다** |
| "세그먼트들이 이어져 미앤더를 이룬다" | 부분적으로 맞다 — 이어지는 건 맞지만(체인), 미앤더(M 구간)는 앵커 사이 경로 일부일 뿐 전체가 미앤더는 아니다(PF 구간은 그냥 직선) |
| 우리 RT(v1/v2): "블록 하나 = 그 안에서 꽉 채워 접는 하드 컨테이너" | QPlacer 원본보다 **더 엄격한 근사** — 원본은 블록을 앵커로만 쓰고 실제 접기는 전역 목표 길이에 맞춰 유연하게 배분하는데, 우리는 각 셀을 독립적으로 최대 용량까지 채우려 든다 |

이번 조사로 새로 확정된 건 없고("면적 공식은 맞다"는 이미 알던 사실을 재확인), **RT의
"셀 = 하드 컨테이너" 가정이 QPlacer 원본의 유연한 체인+PF/M 모델보다 더 제약적**이라는
점이 새로 드러났다 — 이게 RT v2가 리드(port→cell) 구간에서 교차를 겪는 근본 이유와도
통한다(포트가 own_cells 안에 없으면 그 구간은 하드 컨테이너 모델이 아예 표현 못 하는
경로다 — QPlacer식 PF 구간이 정확히 이 역할을 하도록 설계돼 있었다).

## 발견 2 — num_segments 재산정: 길이 기준으로 바꾸면 박스가 칩보다 커진다

발견 1의 결론(면적 공식이 맞다)에도 불구하고, 요청대로 "길이 기준"(`num_segments =
ceil(l / segment_size_um)`, 즉 "블록 하나 = 체인 한 칸 = `segment_size_um`만큼의
직선 진행거리"로 보는 대안)을 계산했다.

| 칩 | 커플러 | 면적기준 합계(현재) | 길이기준 합계(대안) | 배수 |
|---|---|---|---|---|
| grid_25 | 40 | 384 | 1,842 | 4.80x |
| xtree_53 | 52 | 504 | 2,420 | 4.80x |
| falcon | 28 | 269 | 1,293 | 4.81x |
| eagle | 144 | 1,381 | 6,634 | 4.80x |
| aspen_11 | 48 | 460 | 2,216 | 4.82x |
| aspen_m | 106 | 1,012 | 4,871 | 4.81x |

**배수가 정확히 `segment_size_um/meander_spacing_um = 200/40 = 5`에 근접**한다(올림
때문에 4.80~4.82x) — 사용자가 예시로 든 "10 → 46"(4.6x)과 같은 규모다.

**파급(박스 면적)**: 길이 기준 세그먼트 하나도 여전히 자기 몫의 `segment_size_um²`
칸이 필요하다고 가정하면(현재 GP 격자 모델을 그대로 쓸 경우), 박스에 필요한 총 면적이
정확히 5배가 된다(세그먼트 칸 크기는 그대로인데 칸 개수만 5배 늘어나므로):

| 칩 | 칩 면적 | 현재 박스 총면적(칩 대비) | 길이기준 박스 총면적(칩 대비) | 배수 |
|---|---|---|---|---|
| grid_25 | 6.08e7 | 1.46e7 (24.0%) | 7.30e7 (**120.0%**) | 5.00x |
| xtree_53 | 8.65e7 | 1.91e7 (22.1%) | 9.57e7 (**110.7%**) | 5.00x |
| falcon | 4.76e7 | 1.02e7 (21.5%) | 5.11e7 (**107.3%**) | 5.00x |
| eagle | 2.22e8 | 5.25e7 (23.6%) | 2.63e8 (**118.2%**) | 5.00x |
| aspen_11 | 9.99e7 | 1.76e7 (17.6%) | 8.77e7 (87.8%) | 5.00x |
| aspen_m | 1.66e8 | 3.86e7 (23.2%) | 1.93e8 (**115.8%**) | 5.00x |

**6칩 중 5칩에서 필요 박스 총면적이 칩 전체 면적을 넘는다(107.3~120.0%)** — 겹침 없이
다 배치하는 게 산수로 불가능하다(다른 커플러 박스는커녕 큐빗 자리조차 남지 않는다).
**길이 기준으로 세그먼트 개수를 5배 늘리는 건, 지금의 "세그먼트 칸 = 고정 면적" GP
모델과는 근본적으로 양립 불가능하다** — 발견 1이 보여준 QPlacer 원본의 실제 설계(블록을
성기게 체인으로 뿌리고 접기는 별도 M 구간에 맡김)로 전면 재설계하지 않는 한, 길이 기준
공식으로의 단순 교체는 이 파이프라인을 깬다.

## 발견 3 — 자기 큐빗 겹침이 실제로 심각하다 (신규 발견)

렌더링에서 세그먼트가 큐빗 위에 얹힌 게 보인다는 지적을 실측했다(GP 직후, DP 이전 —
좌표가 가장 이해하기 쉬운 시점).

| 칩 | 자기 큐빗과 겹치는 세그먼트 | 큐빗 몸체 안에 완전히 들어간 세그먼트 | 침투 깊이(중심 기준) avg/max |
|---|---|---|---|
| grid_25 | 85 | 2 | 39.1 / **121.1**um |
| xtree_53 | 103 | 1 | 45.1 / **103.0**um |
| falcon | 102 | 6 | 52.3 / **127.2**um |
| eagle | 506 | 28 | 53.8 / **177.8**um |
| aspen_11 | 167 | 10 | 58.5 / **165.1**um |
| aspen_m | 281 | 14 | 48.7 / **182.3**um |

큐빗 반폭은 200um이다 — **eagle/aspen_m에서 세그먼트 중심이 큐빗 경계에서 최대
177~182um(반폭의 88~91%) 안쪽까지 들어간다**, 사실상 큐빗 중심 코앞이다. 6칩 전부에서
"큐빗 몸체에 완전히 파묻힌" 세그먼트가 1~28개씩 나온다.

**원인**: `core/globalplacement.py`의 `_cell_blocked()`가 `Coupler.is_own_qubit()`
로 자기 큐빗을 장애물에서 뺄 때, **포트 근처인지 몸체 중앙인지를 구분하지 않는다** —
자기 q1/q2가 걸리는 셀이면 어디든(코너든 중앙이든) 무조건 통과시킨다. 물리적으로는
커플러의 CPW 트레이스가 자기 큐빗의 **포트(패드)** 한 점에만 붙으면 되고, 큐빗 포켓
안쪽(트랜스몬 커패시터/조셉슨 접합이 있는 자리)을 가로질러 지나갈 이유가 없다 —
`docs/20260917_gp_qc_findings.md` 발견 4가 "own-qubit 겹침은 설계상 정상"이라고 판단한
근거는 "포트에서 커플러가 시작하니 그 지점 겹침은 당연하다"였는데, 그 판단이 "포트
지점"을 "큐빗 몸체 전체"로 지나치게 넓게 적용된 것으로 보인다. **이건 실제 버그/설계
결함으로 봐야 한다** — 포트에 닿기 위해 큐빗 경계선 근처 셀 하나 정도가 걸리는 것과,
세그먼트가 큐빗 중심 코앞까지 파고드는 것은 전혀 다른 문제다.

## 발견 4 — 박스-큐빗 겹침: 자기 큐빗은 100% 겹침(정상), 제3자는 실배치 0(기존 결론 유지)

| 칩 | 박스 수 | 자기 큐빗과 겹침 | 제3자 큐빗과 겹침(박스 기준) |
|---|---|---|---|
| grid_25 | 40 | 40 (100%) | 0 (0%) |
| xtree_53 | 52 | 52 (100%) | 7 (13%) |
| falcon | 28 | 28 (100%) | 7 (25%) |
| eagle | 144 | 144 (100%) | 13 (9%) |
| aspen_11 | 48 | 48 (100%) | 5 (10%) |
| aspen_m | 106 | 106 (100%) | 4 (4%) |

**두 판정 다 맞다 — 서로 다른 걸 쟀을 뿐이다.** 이전 진단(`docs/20260920_rt_segment_
breakdown.md` 발견 4)이 "실제 겹침 0"이라 한 건 **제3자 큐빗 기준**이었고, 그건 지금도
맞다(박스는 13~25%가 겹쳐 보여도 실제 배치된 세그먼트가 제3자 큐빗을 밟은 사례는
0이었다 — GP의 `_cell_blocked`가 셀 단위로 걸러내므로). **렌더링에서 보이는 "눈에 띄는"
겹침은 대부분 자기 큐빗과의 겹침(100%, 발견 3에서 확인했듯 실제로 문제)이다** — 박스
자체가 포트에서 시작해 자기 큐빗 영역을 포함하도록 정의돼 있으니 박스 레벨 겹침은
당연하지만, 그 안에서 실제로 세그먼트가 큐빗 중심까지 파고드는 건(발견 3) 별개의,
새로 확인된 문제다.

## 결론과 권고

1. **면적 기준 num_segments 공식은 QPlacer를 정확히 반영한 것으로 확인됐다** — 길이
   기준으로 바꾸지 말 것. 바꾸면 5칩에서 박스 총면적이 칩 면적을 넘어(최대 120%) 현재
   GP 모델(세그먼트=고정 면적 칸) 자체가 성립하지 않는다.
2. **RT(v1/v2 둘 다)가 "셀 = 하드 컨테이너"로 접는 건 QPlacer 원본보다 엄격한 근사다**
   — 원본은 블록을 체인 앵커로만 쓰고 실제 접기(M)는 PF 구간과 섞어 전역 목표 길이에
   맞춰 유연하게 배분한다. 이 차이가 RT v2의 리드 구간 교차 문제와 직접 연결된다(포트가
   own_cells 밖에 있을 때 하드 컨테이너 모델은 그 구간을 표현할 방법이 없다 — QPlacer의
   PF 구간이 정확히 이 역할을 하도록 설계됐었다). 이번 조사 범위는 진단까지다 — RT를
   "체인 + PF/M" 모델로 재설계할지는 별도 결정 사항으로 남긴다.
3. **자기 큐빗 겹침은 실제 문제다 — 규모부터 다시 봐야 한다.** 포트 근처 셀 하나 정도
   걸치는 것과, 세그먼트가 큐빗 중심 코앞(반폭의 91%)까지 들어가는 것은 다르다. GP의
   own-qubit 예외를 "포트 근접 셀만" 허용하도록 좁히는 방향이 물리적으로 타당해
   보인다 — 다만 그러면 GP의 박스 용량이 줄어(own-qubit 셀 중 상당수가 다시 장애물이
   되므로) 실패율이 올라갈 수 있다(`docs/20260917_gp_qc_findings.md` 발견 4가 이미
   "own-qubit 겹침 중 95.8%가 실패 원인이었다"는 반대 방향 실측을 남겨뒀다 — 이번엔
   반대로 "너무 관대했다"는 쪽의 재조정이 필요해 보인다). 수정 방향과 폭은 별도 결정
   사항으로 남긴다.
4. **박스-큐빗 겹침의 두 판정(자기 100%, 제3자 0)은 둘 다 정확하다** — 모순이 아니라
   서로 다른 질문의 답이었다.

## 재현

```bash
# 발견 1: QPlacer 소스 직접 확인
sed -n '440,472p' /home/LabMember/ngchoi/research/ref/80_QPlacer/Qplacer/qplacer/qplacer_bm/design_format.py   # num_poly 공식
sed -n '277,330p' /home/LabMember/ngchoi/research/ref/80_QPlacer/Qplacer/qplacer/qplacer_bm/design_format.py   # create_pos_matrix (블록 배치)
sed -n '83,158p'  /home/LabMember/ngchoi/research/ref/80_QPlacer/Qplacer/qplacer/qplacer_bm/design_format.py   # create_metal_resonator (죽은 코드, PF/M)

# 발견 2/3/4: 6칩 실측 (코드 수정 없음)
python3 - <<'PY'
import sys, math
sys.path.insert(0, "src")
from core.state import ChipState
from core.floorplan import Floorplan
from core.globalplacement import GlobalPlacement
from core.legalization import Legalization
from utils.params import Params
from utils.parsing import Parser

params = Params()
bench = Parser().load_configs(["grid_25", "xtree_53", "falcon", "eagle", "aspen_11", "aspen_m"])
states = [ChipState.processor_config(c, params) for c in bench]
states = Floorplan(params).run(states)
states = GlobalPlacement(params).run(states)

for s in states:
    total_area = sum(c.l * params.meander_spacing_um for c in s.couplers.values())
    total_length_based_area = sum(c.l * params.segment_size_um for c in s.couplers.values())
    print(s.processor_name, total_area, total_length_based_area,
          total_length_based_area / (s.chip_width * s.chip_height))
PY
```
