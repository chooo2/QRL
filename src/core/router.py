"""Router(RT) v3: 세그먼트 체인 연결 — 셀 단위 comb(v2, _cell_meander)을 버린다.

v2는 "세그먼트 하나 = 그 안에서 미앤더를 접는 셀"이라는 잘못된 전제로 만들어졌다. QPlacer
원본(Zhang et al., ISCA'25) Figure 8과 그 소스(design_format.py, docs/20260920_segment_
model_review.md 발견 1에서 이미 직접 대조 확인함)가 보여주는 실제 모델은 다르다:

  (c) 공진기를 l_b 길이 조각(wireblk)으로 자른다 — 세그먼트 하나 = 배선 한 토막(체인 링크).
  (d) 그 조각들을 인접하게 배치한다(직선/L자/블록 — create_pos_matrix가 q1-q2 사이 전체
      거리에 걸쳐 성기게 뿌리는 앵커, 본 저장소의 GP boustrophedon 채우기와 같은 원리).
  (e) 인접한 조각들을 그대로 이으면 미앤더 모양이 "저절로" 나온다 — 각 조각 내부에서
      따로 접는 게 아니라, 조각들의 지그재그 배치 자체가 미앤더를 만든다.

v2의 _cell_meander(셀 하나 안에서 진입/이탈 변 사이를 pitch 간격 comb으로 채움)는 이
오해의 산물이었다 — 실측 결과 교차가 0~6건에서 0~51건으로 늘었다(셀 내부를 촘촘히 채우며
꺾다 보니 인접 셀의 comb과 만나는 경우가 잦아짐). 이번 버전은 (e)를 있는 그대로 구현한다:
세그먼트 중심을 idx 순서로 그냥 잇는다. 폴딩(접기)이 전혀 없다.

## 연결점: 세그먼트 중심 (셀 경계 중점 아님)

두 후보를 검토했다:
  - 셀 경계 중점(진입변 중점 -> 이탈변 중점): 직선 통과 셀은 어차피 중심을 지나는 직선과
    거의 같은 길이지만, 꺾이는 셀에서는 L자로 짧아져 "세그먼트 개수 x segment_size_um"라는
    예측과 어긋나는 셀별 가변 길이가 생긴다.
  - 세그먼트 중심(채택): 4-인접한 두 세그먼트의 중심 간 거리는 항상 정확히 segment_size_um
    이다(격자 간격이 세그먼트 크기와 같으므로) — 그래서 체인 구간 전체 길이가 정확히
    "(세그먼트 개수-1) x segment_size_um"로 나온다. 2절(길이 검증)이 요구하는 예측과
    오차 없이 맞아떨어지는 쪽을 골랐다 — 셀 경계 중점은 이 예측 자체를 흐린다.

## idx 점프는 실패다 (폴백 없음)

GP가 boustrophedon으로 채운 순서(Segment.idx)를 그대로 체인 순서로 쓴다. idx가 연속인
두 세그먼트가 격자상 4-인접(정확히 한 셀, segment_size_um)하지 않으면(점프 — 이전 조사
기준 전체 idx 쌍의 0~3.5%, docs/20260920_rt_v2_redesign.md) 그 자리에서 즉시 실패
처리한다 — BFS로 다리를 놓거나 다른 경로를 찾지 않는다. v1이 blocked 상태를 무시하고
포트-포트 직선으로 폴백하던 게 교차의 원인이었다(docs/20260920_rt_findings.md) — 그
교훈을 그대로 지킨다: 막히면 막힌 대로 보고한다.

## 길이는 억지로 안 맞춘다 — 그리고 실제로 안 맞는다

경로 길이는 세그먼트 개수만으로 결정된다: 체인 구간 (k-1) x segment_size_um + 리드인/
리드아웃(포트<->첫/마지막 세그먼트, 가변). num_segments = ceil(l x meander_spacing_um /
segment_size_um^2)의 차원분석이 강제하는 "세그먼트 하나가 담당하는 길이"는
segment_size_um^2 / meander_spacing_um = 200^2/40 = 1000um인데(docs/20260920_segment_
model_review.md 발견 1), 체인 연결은 세그먼트 하나당 segment_size_um=200um만 기여한다 —
5배 차이. 실측(아래 6절)으로 확인되는 대로: 목표 l 대비 대략 -80%대 미달이 나온다. 이건
버그가 아니라 모델 불일치다 — num_segments 공식은 "세그먼트가 자기 면적을 접어서 채운다"
(2D 인필)는 QPlacer 원본과 다른 가정(2026-09-20 이전 이 저장소의 실제 구현, v1/v2)을
전제로 도출됐는데, QPlacer의 실제 배치(create_pos_matrix)는 세그먼트를 "성긴 체인 앵커"로
만 쓰고, 목표 길이를 채우는 진짜 접기는 앵커 사이 별도 구간(죽은 코드 create_metal_
resonator가 보여주는 PF/M 교대 — PathFinder 직선 구간과 Meander 접기 구간)에서 일어난다.
이 저장소 파이프라인엔 그 PF/M 구간에 해당하는 단계 자체가 없다 — 그래서 세그먼트 개수를
QPlacer 공식 그대로 쓰면서 폴딩만 빼면 구조적으로 미달할 수밖에 없다. 억지로 접어 채우지
않는다(그러면 v2로 되돌아가 다시 교차가 는다) — 미달을 그대로 보고하고 원인(모델 불일치)
을 지목하는 것이 이 버전의 역할이다. 채우려면 PF/M 구간 자체를 다음 단계로 도입해야 한다.

## GP/RT 인터페이스 한계는 여전하다

리드인/리드아웃(포트 <-> 체인 첫/마지막 세그먼트)은 own_cells 밖을 직선으로 지나갈 수
있다(docs/20260920_rt_findings.md 발견 3) — 이 재설계로도 해소되지 않는다. 다만 체인
구간(리드가 아닌 leg)은 항상 자기 세그먼트 두 개의 중심을 잇는 직선이라, 그 두 세그먼트가
차지하는 정사각형 합집합 안에서만 움직인다 — 다른 커플러의 세그먼트와 겹치지 않으므로
(GP가 이미 보장) 리드가 아닌 leg끼리는 서로 다른 커플러 사이에서 교차할 수 없다
(find_crossing_legs의 최적화 근거, v2와 동일한 증명이 그대로 성립한다).
"""
import logging
import math
from dataclasses import replace

from core.floorplan import segments_cross
from core.globalplacement import _coupler_order
from core.state import ChipState, Coupler, Qubit, Segment


class Router:
    def __init__(self, params):
        self.params = params
        # RT의 실패 단위는 "커플러 하나"이지 "칩 전체"가 아니다 — FP/GP/LG/DP처럼 칩
        # 통째로 버릴 조건이 RT엔 없다. 그래도 인터페이스는 형제 단계들과 맞춘다 — 이
        # 필드는 항상 비어 있는 게 정상이다.
        self.skipped: list[tuple[str, str]] = []
        # 칩 이름 -> [(커플러 키, 실패 사유), ...].
        self.route_failures: dict[str, list[tuple[tuple[int, int], str]]] = {}

    def run(self, states: list[ChipState]) -> list[ChipState]:
        out: list[ChipState] = []
        for state in states:
            out.append(self._route(state))
        return out

    def _route(self, state: ChipState) -> ChipState:
        order_strategy = getattr(self.params, "rt_coupler_order", "shortest_first")
        # GlobalPlacement의 순서 로직을 재사용한다 — RT 커플러 사이엔 공유 상태가 없어
        # 순서 자체는 결과에 영향이 없지만, 인터페이스 일관성과 로그 순서를 위해 유지한다.
        order = _coupler_order(state.couplers, state.qubits, state.port_assignment, order_strategy)

        new_couplers: dict[tuple[int, int], Coupler] = {}
        failures: list[tuple[tuple[int, int], str]] = []
        n_with_segments = 0
        for key in order:
            coupler = state.couplers[key]
            if not coupler.segments:
                # GP가 세그먼트를 하나도 못 놓은 커플러 — waypoints를 비워둔 채 그대로
                # 넘긴다(4절: fallback 없음, 세그먼트가 없으면 라우팅할 대상 자체가 없다).
                new_couplers[key] = coupler
                continue
            n_with_segments += 1
            wp, reason = _route_coupler(coupler, state.qubits, state.port_assignment.get(key))
            if wp is None:
                failures.append((key, reason))
                new_couplers[key] = coupler
            else:
                new_couplers[key] = replace(coupler, waypoints=wp)

        if failures:
            self.route_failures[state.processor_name] = failures

        n_routed = n_with_segments - len(failures)
        logging.info(
            "[RT] %s: couplers=%d 세그먼트있음=%d 라우팅성공=%d 실패=%d",
            state.processor_name, len(state.couplers), n_with_segments, n_routed, len(failures),
        )
        if failures:
            by_reason: dict[str, int] = {}
            for _key, reason in failures:
                by_reason[reason] = by_reason.get(reason, 0) + 1
            logging.info("[RT] %s: 실패 사유별 건수 %s", state.processor_name, by_reason)

        return replace(state, couplers=new_couplers)


# ---------------------------------------------------------------------------
# 커플러 하나 라우팅 — 세그먼트 체인을 idx 순서로 그냥 잇는다. 접지 않는다.
# ---------------------------------------------------------------------------

def _route_coupler(
    coupler: Coupler, qubits: dict[int, Qubit], assignment: tuple[str, str] | None,
) -> tuple[list[tuple[float, float]] | None, str | None]:
    if assignment is None:
        return None, "no_port_assignment"

    p1_name, p2_name = assignment
    q1, q2 = qubits[coupler.q1], qubits[coupler.q2]
    port1, port2 = q1.ports[p1_name], q2.ports[p2_name]

    cell = coupler.segment_size_um
    segs_sorted = sorted(coupler.segments, key=lambda s: s.idx)

    for a, b in zip(segs_sorted, segs_sorted[1:]):
        if not _grid_adjacent(a, b, cell):
            # idx 순서상 다음 세그먼트가 격자상 이웃이 아니다(GP boustrophedon의 점프) —
            # 다리를 놓지 않고 그 자리에서 실패 처리한다(모듈 docstring "idx 점프는
            # 실패다" 참고).
            return None, "chain_broken"

    centers = [(s.x, s.y) for s in segs_sorted]

    # 체인의 두 끝(idx 순서상 처음/마지막) 중 어느 쪽을 포트1에 붙일지: 총 리드 거리
    # (포트-끝점 직선 두 개의 합)가 더 짧은 방향을 고른다 — idx 순서 자체(=방문하는
    # 세그먼트 집합)는 뒤집어도 그대로이므로 길이/커버리지에 영향 없이 리드만 최소화한다.
    c0, c1 = centers[0], centers[-1]
    cost_fwd = math.hypot(port1[0] - c0[0], port1[1] - c0[1]) + math.hypot(port2[0] - c1[0], port2[1] - c1[1])
    cost_rev = math.hypot(port1[0] - c1[0], port1[1] - c1[1]) + math.hypot(port2[0] - c0[0], port2[1] - c0[1])
    if cost_rev < cost_fwd:
        centers = list(reversed(centers))

    wp = _dedupe_close([port1, *centers, port2])
    if len(wp) < 2:
        return None, "degenerate_path"
    return wp, None


# 두 세그먼트가 전역 격자상 정확히 한 칸(segment_size_um) 이웃인지 — 대각선/2칸 이상
# 떨어짐은 이웃 아님. eps는 부동소수 표현 오차만 흡수한다(다른 물리량 대비 1e-6 상대
# 오차라 실제 어긋남을 가릴 걱정이 없다, core/state.py의 AABB_EPS_UM과 같은 원칙).
def _grid_adjacent(a: Segment, b: Segment, cell: float) -> bool:
    dx, dy = abs(a.x - b.x), abs(a.y - b.y)
    eps = cell * 1e-6
    return (abs(dx - cell) < eps and dy < eps) or (abs(dy - cell) < eps and dx < eps)


def _polyline_length(points: list[tuple[float, float]]) -> float:
    return sum(
        math.hypot(points[i + 1][0] - points[i][0], points[i + 1][1] - points[i][1])
        for i in range(len(points) - 1)
    )


def _dedupe_close(points: list[tuple[float, float]], eps: float = 1e-6) -> list[tuple[float, float]]:
    out = [points[0]]
    for p in points[1:]:
        if math.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) > eps:
            out.append(p)
    return out


# 폴리라인의 실제 방향 전환 횟수 — v2는 "꺾임"이 탐색 목적이 아니라 순수 보고 지표였고
# (셀마다 강제로 여러 번 꺾었으므로 "최소화"가 의미 없었다), v3도 마찬가지로 설계 목표가
# 아니라 보고용이다. 다만 이제 각 leg가 실제 직선 구간이라(더 이상 셀 내부 comb의 일부가
# 아니라) 방향이 실제로 바뀌는 지점만 센다 — 연속 leg의 단위 방향벡터가 같지 않으면 꺾임.
def _count_turns(wp: list[tuple[float, float]]) -> int:
    prev_dir = None
    turns = 0
    for i in range(len(wp) - 1):
        dx, dy = wp[i + 1][0] - wp[i][0], wp[i + 1][1] - wp[i][1]
        norm = math.hypot(dx, dy)
        if norm < 1e-9:
            continue
        d = (dx / norm, dy / norm)
        if prev_dir is not None and (prev_dir[0] * d[0] + prev_dir[1] * d[1]) < 1.0 - 1e-6:
            turns += 1
        prev_dir = d
    return turns


# ---------------------------------------------------------------------------
# 검증 — 교차 개수 / 세그먼트 이용률 (5절 보고용, main.py나 검증 스크립트가 호출)
# ---------------------------------------------------------------------------

# 실제 배선 leg(꺾은선의 각 구간)끼리 교차하는 쌍을 비교한다. 정점(포트) 공유 쌍도 빼지
# 않고 전부 비교한다(core/floorplan.py의 assign_ports 모듈 주석에 적힌 전례와 같은 원칙).
#
# 전수(O(전체 leg 수^2)) 대신 "리드 leg(각 커플러의 첫/마지막 leg) vs 전체"만 비교한다 —
# v3에서도 이게 완전하다: 리드가 아닌 체인 leg는 언제나 자기 세그먼트 두 개(4-인접, 둘 다
# 이 커플러 소유)의 중심을 잇는 직선이라 그 두 세그먼트 정사각형의 합집합 안에서만
# 움직이고, GP가 서로 다른 커플러의 세그먼트를 겹치지 않게 배치하므로 그 합집합은 다른
# 커플러의 세그먼트와 겹치지 않는다 — 그래서 "리드가 아닌 leg 대 리드가 아닌 leg" 쌍은
# 서로 다른 커플러 사이에서 수학적으로 교차할 수 없다(모듈 docstring 참고). 리드는
# 커플러당 정확히 2개뿐이라 이 최적화로 O(리드 수 x 전체 leg 수)로 줄인다.
def find_crossing_legs(state: ChipState) -> tuple[int, dict[tuple[int, int], set[int]]]:
    all_legs: list[tuple[tuple[int, int], int, tuple[float, float], tuple[float, float]]] = []
    lead_legs: list[tuple[tuple[int, int], int, tuple[float, float], tuple[float, float]]] = []
    for key in sorted(state.couplers):
        wp = state.couplers[key].waypoints
        n_legs = len(wp) - 1
        for i in range(n_legs):
            entry = (key, i, wp[i], wp[i + 1])
            all_legs.append(entry)
            if i == 0 or i == n_legs - 1:
                lead_legs.append(entry)

    flagged: dict[tuple[int, int], set[int]] = {}
    n_cross = 0
    counted: set[frozenset] = set()
    for ka, ia, pa1, pa2 in lead_legs:
        for kb, ib, pb1, pb2 in all_legs:
            if ka == kb:
                continue
            pair_id = frozenset(((ka, ia), (kb, ib)))
            if pair_id in counted:
                continue
            if segments_cross(pa1, pa2, pb1, pb2):
                counted.add(pair_id)
                n_cross += 1
                flagged.setdefault(ka, set()).add(ia)
                flagged.setdefault(kb, set()).add(ib)
    return n_cross, flagged


# find_crossing_legs가 찾은 교차 쌍 각각을 (lead, lead)/(lead, chain)/(chain, chain)으로
# 분류한다(3절: "lead 구간인지 체인 구간인지"). chain-chain 교차는 위 최적화 증명이 서로
# 다른 커플러 사이에서 불가능함을 보장하므로, 나온다면 leg 인덱싱 버그를 뜻한다 — 그래서
# 이 분류 자체가 그 불변식의 실행 시점 점검이기도 하다.
def classify_crossings(state: ChipState) -> dict[str, int]:
    counts = {"lead_lead": 0, "lead_chain": 0, "chain_chain": 0}
    legs_by_key: dict[tuple[int, int], list[tuple[tuple[float, float], tuple[float, float]]]] = {}
    for key in sorted(state.couplers):
        wp = state.couplers[key].waypoints
        legs_by_key[key] = [(wp[i], wp[i + 1]) for i in range(len(wp) - 1)]

    keys = sorted(legs_by_key)
    for ai in range(len(keys)):
        ka = keys[ai]
        legs_a = legs_by_key[ka]
        for bi in range(ai + 1, len(keys)):
            kb = keys[bi]
            legs_b = legs_by_key[kb]
            for ia, (pa1, pa2) in enumerate(legs_a):
                a_is_lead = ia == 0 or ia == len(legs_a) - 1
                for ib, (pb1, pb2) in enumerate(legs_b):
                    if not segments_cross(pa1, pa2, pb1, pb2):
                        continue
                    b_is_lead = ib == 0 or ib == len(legs_b) - 1
                    if a_is_lead and b_is_lead:
                        counts["lead_lead"] += 1
                    elif a_is_lead or b_is_lead:
                        counts["lead_chain"] += 1
                    else:
                        counts["chain_chain"] += 1
    return counts


# 세그먼트 이용률 — 이 커플러의 세그먼트 중 실제로 경로의 꼭짓점으로 쓰인 비율. v3(순수
# 체인 연결)에서는 라우팅에 성공하면 구조적으로 항상 1.0이다(모든 세그먼트가 정확히 하나의
# 꺾은선 꼭짓점이 되고, v2처럼 "지나가기만 하고 안 접은 셀"이나 "안 지나간 own_cell"이 있을
# 수 없다) — 그래도 가정이 아니라 실측으로 확인한다(5절 원칙: 측정 없이 단정하지 않는다).
def cell_utilization(coupler: Coupler) -> float | None:
    if not coupler.segments:
        return None
    if not coupler.waypoints:
        return 0.0
    seg_centers = {(round(s.x, 6), round(s.y, 6)) for s in coupler.segments}
    wp_points = {(round(x, 6), round(y, 6)) for x, y in coupler.waypoints}
    return len(seg_centers & wp_points) / len(seg_centers)


# 칩 하나의 라우팅 결과를 요약한다 — 5절 보고용. 실행 시간은 여기 포함하지 않는다 —
# main.py의 run_stage()가 이미 스테이지별 칩별 시간을 재고 있다.
def route_report(state: ChipState) -> dict:
    couplers = state.couplers
    with_segments = [c for c in couplers.values() if c.segments]
    routed = [c for c in with_segments if c.waypoints]
    turns = [_count_turns(c.waypoints) for c in routed]
    length_err_pct = [(c.route_length - c.l) / c.l * 100.0 for c in routed]
    n_cross, _ = find_crossing_legs(state)
    utils = [cell_utilization(c) for c in with_segments]
    utils = [u for u in utils if u is not None]

    return {
        "processor": state.processor_name,
        "couplers": len(couplers),
        "with_segments": len(with_segments),
        "routed": len(routed),
        "failed": len(with_segments) - len(routed),
        "crossings": n_cross,
        "crossings_by_kind": classify_crossings(state),
        "turns_avg": (sum(turns) / len(turns)) if turns else 0.0,
        "turns_max": max(turns) if turns else 0,
        "len_err_avg_pct": (sum(length_err_pct) / len(length_err_pct)) if length_err_pct else 0.0,
        "len_err_max_abs_pct": max((abs(e) for e in length_err_pct), default=0.0),
        "within_5pct_frac": (sum(1 for e in length_err_pct if abs(e) <= 5.0) / len(length_err_pct)) if length_err_pct else 0.0,
        "cell_utilization_avg": (sum(utils) / len(utils)) if utils else 0.0,
    }
