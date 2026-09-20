"""Router(RT) v2: 커플러별 실제 배선 경로(Coupler.waypoints) 확정 — 전수 순회 + 셀 단위 미앤더.

v1(꺾임-최소 A* + leg 단위 지그재그)은 목표 길이의 21~33%밖에 못 채웠다
(docs/20260920_dp_alpha_diagnosis.md, docs/20260920_rt_segment_breakdown.md) — 정확히
곱해지는 두 손실 때문이었다: (A) 꺾임-최소 경로가 own_cells의 평균 48%만 지나감(나머지는
아예 안 씀), (B) leg를 따라가는 대각선 지그재그 한 줄이 이론 용량(셀당 cell²/pitch)의
60.6~69.2%만 뽑아냄. 이 v2는 두 손실을 각각 정면으로 제거한다:

  (A') 꺾임 최소화를 버리고 **자기 세그먼트 전부를 지나는 순회 경로**를 만든다(_build_tour).
       GP가 boustrophedon으로 셀을 채운 순서(Segment.idx)를 그대로 순회 순서로 쓴다 —
       측정 결과(재현 커맨드 참고) idx가 인접하지 않는 "점프"는 전체 idx 쌍의 0~3.5%뿐이고
       (6칩 평균), 그 점프는 own_cells 안에서 BFS로 짧게(최대 3칸) 잇는다. own_cells가
       4-연결로 안 이어지면(BFS 실패) 그대로 실패 처리한다 — 5절 원칙(막힌 걸 억지로 안
       뚫는다) 유지. 꺾임 수는 더 이상 탐색 목적이 아니라 route_report()의 보고 지표로만
       남는다(실측: 미앤더가 셀 하나당 수십 번 꺾으므로 "최소 꺾임"이라는 목적 자체가
       무의미했다 — segment_breakdown.md 발견 5).
  (B') leg를 따라가는 지그재그 대신, **셀 하나를 지날 때마다 그 안에서 접는다**(_cell_meander).
       진입/이탈 방향(4방향)에 따라 "직진 통과"(대변으로 빠짐)와 "꺾임 통과"(인접 변으로
       빠짐) 두 형태를 자동으로 만든다 — 두 경우 다 폭 방향으로 여러 줄(comb)을 채우는
       진짜 미앤더라, leg-지그재그(대각선 한 줄)보다 구조적으로 효율적이다(실측:
       cell=200/pitch=40에서 이론값 1,000um 대비 128~131%까지 나옴 — 재현 커맨드 참고).

길이 맞추기(3절)는 그대로 별도 단계다: 셀별 "0접기(그냥 통과) 길이"와 "최대접기 길이"를
전부 계산해 두고, 목표 l에 못 미치면(총 최대접기 길이 < l) 전부 최대로 접어 그 미달을
정직하게 보고한다. 남으면(총 최대접기 길이 >= l) 순회 순서대로 필요한 만큼만 접고 나머지는
통과만 시킨다 — 남는 셀에 억지로 더 접어 넣지 않는다.

GP/RT 인터페이스의 근본적 한계(리드인/아웃이 own_cells 밖을 지날 수 있는 문제,
docs/20260920_rt_findings.md 발견 3)는 이 재설계로도 해소되지 않는다 — 포트가 own_cells
"안"에 없으면 그 구간은 여전히 직선 보간이다. 다만 이제 own_cells 거의 전부가 실제로
쓰이므로(발견 A' 해소) 전체 길이에서 리드 구간이 차지하는 비중은 상대적으로 작아진다.
"""
import logging
import math
from collections import deque
from dataclasses import replace

from core.floorplan import segments_cross
from core.globalplacement import _coupler_order
from core.state import ChipState, Coupler, Qubit, Segment


class Router:
    def __init__(self, params):
        self.params = params
        # RT의 실패 단위는 "커플러 하나"이지 "칩 전체"가 아니다(아래 route_failures가 그
        # 단위로 기록한다) — FP/GP/LG/DP처럼 칩 통째로 버릴 조건이 RT엔 없다. 그래도
        # 인터페이스는 형제 단계들과 맞춘다 — 이 필드는 항상 비어 있는 게 정상이다.
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
        # (docs/20260920_rt_findings.md 발견 1) 순서 자체는 결과에 영향이 없지만, 인터페이스
        # 일관성과 로그 순서를 위해 유지한다.
        order = _coupler_order(state.couplers, state.qubits, state.port_assignment, order_strategy)

        new_couplers: dict[tuple[int, int], Coupler] = {}
        failures: list[tuple[tuple[int, int], str]] = []
        n_with_segments = 0
        for key in order:
            coupler = state.couplers[key]
            if not coupler.segments:
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
# 커플러 하나 라우팅
# ---------------------------------------------------------------------------

# 셀 크기(segment_size_um, post-DP) 대비 미앤더 스트로크가 셀 경계에서 남기는 여유 —
# 정확히 경계에 닿으면 부동소수 오차로 넘어갈 수 있어(core/state.py의 AABB_EPS_UM과 같은
# 취지, 다만 여긴 설계상 필요한 물리적 여유다) 2.5%를 둔다.
_CELL_MARGIN_FACTOR = 0.025

_OPPOSITE = {"E": "W", "W": "E", "N": "S", "S": "N"}
_ROTATE = {"E": "N", "N": "W", "W": "S", "S": "E"}  # 방어적 폴백(입구==출구일 때)용 90도 회전
_DIR_VEC = {"E": (1, 0), "W": (-1, 0), "N": (0, 1), "S": (0, -1)}


def _route_coupler(
    coupler: Coupler, qubits: dict[int, Qubit], assignment: tuple[str, str] | None,
) -> tuple[list[tuple[float, float]] | None, str | None]:
    if assignment is None:
        return None, "no_port_assignment"

    p1_name, p2_name = assignment
    q1, q2 = qubits[coupler.q1], qubits[coupler.q2]
    port1, port2 = q1.ports[p1_name], q2.ports[p2_name]

    cell = coupler.segment_size_um
    pitch = coupler.meander_spacing_um
    margin = cell * _CELL_MARGIN_FACTOR

    idx_cell, cell_center = _build_idx_cell(coupler.segments, cell)
    own_cells = set(cell_center)

    ordered_by_idx = [idx_cell[i] for i in sorted(idx_cell)]
    tour = _build_tour(ordered_by_idx, own_cells)
    if tour is None:
        # own_cells가 4-연결로 안 이어짐(GP의 idx 점프를 못 이음) — GP는 성공했지만(segments
        # 가 있으므로) RT가 그 안에서 전수 순회 경로를 못 찾은 것.
        return None, "no_path"

    tour = _orient_tour(tour, cell_center, port1, port2)
    tour = _snap_tour_ends(tour, cell_center, own_cells, cell, port1, port2)

    wp = _assemble_route(tour, cell_center, cell, pitch, margin, port1, port2, coupler.l)
    if wp is None or len(wp) < 2:
        return None, "degenerate_path"
    return wp, None


# ---------------------------------------------------------------------------
# 셀 격자 재구성 (idx 포함) — DP의 alpha 닮음변환 이후에도 같은 커플러 세그먼트끼리는
# 정확한 간격으로 정렬돼 있으므로(core/router.py 모듈 docstring), 이 커플러 자신의
# 최소 x/y를 로컬 원점으로 잡아 격자 인덱스를 역산한다. idx도 함께 돌려준다 — 순회
# 순서(Segment.idx, GP의 boustrophedon 채우기 순서)를 그대로 쓰기 위함이다.
# ---------------------------------------------------------------------------

def _build_idx_cell(
    segments: list[Segment], cell: float,
) -> tuple[dict[int, tuple[int, int]], dict[tuple[int, int], tuple[float, float]]]:
    ox = min(s.x for s in segments)
    oy = min(s.y for s in segments)
    idx_cell: dict[int, tuple[int, int]] = {}
    cell_center: dict[tuple[int, int], tuple[float, float]] = {}
    for s in segments:
        i = int(round((s.x - ox) / cell))
        j = int(round((s.y - oy) / cell))
        idx_cell[s.idx] = (i, j)
        cell_center[(i, j)] = (s.x, s.y)
    return idx_cell, cell_center


# ---------------------------------------------------------------------------
# 전수 순회 경로 — GP의 idx 순서를 기본으로, 인접하지 않는 점프만 BFS로 잇는다.
# ---------------------------------------------------------------------------

def _adjacent(a: tuple[int, int], b: tuple[int, int]) -> bool:
    return abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1


def _bfs_bridge(
    start: tuple[int, int], goal: tuple[int, int], own_cells: set[tuple[int, int]],
) -> list[tuple[int, int]] | None:
    if start == goal:
        return [start]
    visited = {start}
    prev: dict[tuple[int, int], tuple[int, int]] = {}
    q = deque([start])
    while q:
        cur = q.popleft()
        for dx, dy in _DIR_VEC.values():
            nxt = (cur[0] + dx, cur[1] + dy)
            if nxt not in own_cells or nxt in visited:
                continue
            visited.add(nxt)
            prev[nxt] = cur
            if nxt == goal:
                path = [goal]
                while path[-1] != start:
                    path.append(prev[path[-1]])
                path.reverse()
                return path
            q.append(nxt)
    return None


# idx 순서(ordered)를 순회 경로로 만든다 — 연속 셀이 이미 인접하면 그대로 잇고, 아니면
# own_cells 안에서 BFS로 최단 다리를 놓는다(실측: 전체 idx 쌍의 0~3.5%만 해당, 6칩 평균
# — docs/20260920_rt_v2_redesign.md 참고). 다리가 이미 방문한 셀을 다시 지나가도
# 된다(그 셀은 두 번째 통과에서는 접지 않는다, _assemble_route의 first_visit 판정) —
# own_cells가 여러 조각으로 쪼개져 다리 자체가 불가능하면 None(호출부가 no_path로 보고).
def _build_tour(
    ordered: list[tuple[int, int]], own_cells: set[tuple[int, int]],
) -> list[tuple[int, int]] | None:
    tour = [ordered[0]]
    for nxt in ordered[1:]:
        cur = tour[-1]
        if _adjacent(cur, nxt):
            tour.append(nxt)
        else:
            bridge = _bfs_bridge(cur, nxt, own_cells)
            if bridge is None:
                return None
            tour.extend(bridge[1:])
    return tour


# 순회 경로의 두 끝(idx 순서상 처음/마지막 셀) 중 어느 쪽을 포트1에, 어느 쪽을 포트2에
# 붙일지 정한다 — "시작/끝이 각자의 포트와 이어져야 한다"(요청 1절)는 요구를 총 리드
# 거리(포트-끝셀 직선 두 개의 합)가 더 짧은 방향으로 만족시킨다. idx 순서 자체(다리 잇기
# 포함)는 두 경우 모두 동일한 셀 집합을 방문하므로 뒤집어도 커버리지는 그대로다.
def _orient_tour(
    tour: list[tuple[int, int]], cell_center: dict[tuple[int, int], tuple[float, float]],
    port1: tuple[float, float], port2: tuple[float, float],
) -> list[tuple[int, int]]:
    c0, c1 = cell_center[tour[0]], cell_center[tour[-1]]

    def d(a, b):
        return math.hypot(a[0] - b[0], a[1] - b[1])

    cost_fwd = d(port1, c0) + d(port2, c1)
    cost_rev = d(port1, c1) + d(port2, c0)
    return tour if cost_fwd <= cost_rev else list(reversed(tour))


# 포트를 담는 자기 셀이 있으면 그 셀(리드 길이 0), 없으면 가장 가까운 자기 셀을 고른다.
def _select_nearest_cell(
    port: tuple[float, float], cell_center: dict[tuple[int, int], tuple[float, float]], cell: float,
) -> tuple[int, int]:
    half = cell / 2.0
    best, best_d2 = None, None
    for c, (cx, cy) in cell_center.items():
        if abs(port[0] - cx) <= half + 1e-6 and abs(port[1] - cy) <= half + 1e-6:
            return c
        d2 = (port[0] - cx) ** 2 + (port[1] - cy) ** 2
        if best_d2 is None or d2 < best_d2:
            best, best_d2 = c, d2
    return best


# idx 순서(+ orient)로 정해진 순회 경로의 두 끝은 own_cells 안에서 "구조적으로 정해진"
# 셀일 뿐, 포트에 제일 가까운 셀이라는 보장이 없다 — GP의 boustrophedon 채우기가 박스
# 한쪽 구석에서만 시작하므로(docs/20260920_rt_findings.md 발견 3) idx=0은 대개 포트1
# 근처지만 idx=k-1은 포트2에서 멀 수 있다. 그 경우 리드 직선이 own_cells 밖을 길게
# 지나가고, 실측 결과 그게 다른 커플러의 내부(접힌) leg와 교차하는 유일한 원인이었다
# (route_report의 crossings 대부분이 "리드 vs 내부" 쌍 — find_crossing_legs 주석 참고).
#
# 고친다: 각 끝을, own_cells 전체 중 그 포트에 실제로 제일 가까운 셀로 "스냅"한다. 그
# 셀은 이미 순회 경로 어딘가에 들어 있다(idx 순서가 own_cells 전부를 방문하므로) — 그냥
# 순회의 끝에서 그 셀까지 BFS로 짧은 다리를 놓아 끝을 옮기는 것뿐이라 커버리지를 잃지
# 않는다(다리가 지나가는 셀은 이미 방문한 셀이라 재방문 처리된다, _assemble_route의
# is_first_visit). own_cells가 끊겨 있어 다리를 못 놓으면(이론상 안 나옴 — own_cells는
# _build_tour가 이미 하나로 연결됨을 확인했다) 원래 끝을 그대로 둔다.
def _snap_tour_ends(
    tour: list[tuple[int, int]], cell_center: dict[tuple[int, int], tuple[float, float]],
    own_cells: set[tuple[int, int]], cell: float,
    port1: tuple[float, float], port2: tuple[float, float],
) -> list[tuple[int, int]]:
    best1 = _select_nearest_cell(port1, cell_center, cell)
    best2 = _select_nearest_cell(port2, cell_center, cell)

    if tour[0] != best1:
        bridge = _bfs_bridge(best1, tour[0], own_cells)
        if bridge is not None:
            tour = bridge[:-1] + tour

    if tour[-1] != best2:
        bridge = _bfs_bridge(tour[-1], best2, own_cells)
        if bridge is not None:
            tour = tour + bridge[1:]

    return tour


def _dir_between(a: tuple[int, int], b: tuple[int, int]) -> str:
    dx, dy = b[0] - a[0], b[1] - a[1]
    for name, (ddx, ddy) in _DIR_VEC.items():
        if (dx, dy) == (ddx, ddy):
            return name
    raise ValueError(f"{a} -> {b} not 4-adjacent")


def _side_point(cx: float, cy: float, half: float, side: str) -> tuple[float, float]:
    dx, dy = _DIR_VEC[side]
    return (cx + dx * half, cy + dy * half)


# 셀 중심에서 target(포트)까지 방향을 4방위 중 가장 가까운 쪽으로 근사 — 순회의 첫/마지막
# 셀이 포트로 나가는 변을 고르는 데 쓴다(own_cells 안에는 이웃 셀이 없으니 dir_between을
# 못 쓴다).
def _nearest_side(cell_xy: tuple[float, float], target: tuple[float, float]) -> str:
    dx = target[0] - cell_xy[0]
    dy = target[1] - cell_xy[1]
    if abs(dx) >= abs(dy):
        return "E" if dx >= 0 else "W"
    return "N" if dy >= 0 else "S"


# ---------------------------------------------------------------------------
# 셀 단위 미앤더(comb) — 직진 통과(entry/exit가 반대 변)와 꺾임 통과(인접 변)를 같은
# 코드로 처리한다. entry_side를 기준으로 쓸어담는 축(가로/세로)을 정하고, n_rows개의
# 변-끝-to-변-끝 스트로크를 pitch 간격으로 쌓은 뒤 exit_side로 빠진다. 진입점/이탈점은
# 항상 그 변의 중점(_side_point)이다 — 그래서 인접한 두 셀은 공유 경계에서 정확히 같은
# 점을 만들어(둘 다 같은 공식으로 그 점을 계산하므로) 이음매 없이 이어진다(검증:
# proto_v2.py의 adjacency continuity check).
#
# n_rows의 홀/짝에 따라 시작 스윕 방향을 바꾼다(near_is_lo, ends_far_if_start_plus1) —
# 그러지 않으면 마지막 스트로크가 진입 변 쪽에서 끝나 출구까지 폭 전체를 되돌아가는
# 낭비가 생긴다(처음 구현에서 실측: n=2일 때 660um -> 방향 보정 후 280um, 같은 진폭/피치
# 로 40% 단축 — proto_v2.py 재현 참고). n_rows=0이면 접지 않고 센터를 경유해 최단으로
# 나간다(재방문 셀이나 예산이 다 떨어진 셀에 쓴다).
def _cell_meander(
    cx: float, cy: float, cell: float, pitch: float, margin: float,
    entry_side: str, exit_side: str, n_rows: int,
) -> list[tuple[float, float]]:
    half = cell / 2.0
    entry = _side_point(cx, cy, half, entry_side)
    exit_ = _side_point(cx, cy, half, exit_side)

    if entry_side == exit_side:
        exit_side = _ROTATE[exit_side]
        exit_ = _side_point(cx, cy, half, exit_side)

    if n_rows <= 0:
        pts = [entry]
        if exit_side != _OPPOSITE.get(entry_side):
            pts.append((cx, cy))
        pts.append(exit_)
        return pts

    horizontal_sweep = entry_side in ("W", "E")
    if horizontal_sweep:
        lo, hi = cx - half + margin, cx + half - margin
        cross_lo, cross_hi = cy - half + margin, cy + half - margin
        entry_along, entry_cross = entry[0], entry[1]
    else:
        lo, hi = cy - half + margin, cy + half - margin
        cross_lo, cross_hi = cx - half + margin, cx + half - margin
        entry_along, entry_cross = entry[1], entry[0]

    def mk(along: float, cross: float) -> tuple[float, float]:
        return (along, cross) if horizontal_sweep else (cross, along)

    span = (n_rows - 1) * pitch
    start = entry_cross - span / 2.0
    if cross_hi - span >= cross_lo:
        start = max(cross_lo, min(start, cross_hi - span))
    else:
        start = cross_lo
    row_positions = [start + k * pitch for k in range(n_rows)]

    near_is_lo = entry_along <= (lo + hi) / 2.0
    ends_far_if_start_plus1 = (n_rows % 2 == 1)
    if near_is_lo:
        direction = 1 if ends_far_if_start_plus1 else -1
    else:
        direction = -1 if ends_far_if_start_plus1 else 1

    pts = [entry, mk(entry_along, row_positions[0])]
    cur_along = entry_along
    for k, rp in enumerate(row_positions):
        target_along = hi if direction == 1 else lo
        pts.append(mk(target_along, rp))
        cur_along = target_along
        if k < n_rows - 1:
            pts.append(mk(target_along, row_positions[k + 1]))
        direction *= -1

    if horizontal_sweep:
        exit_along, exit_cross = exit_[0], exit_[1]
    else:
        exit_along, exit_cross = exit_[1], exit_[0]

    if exit_side == _OPPOSITE.get(entry_side):
        if abs(cur_along - exit_along) > 1e-9:
            pts.append(mk(exit_along, row_positions[-1]))
        pts.append(mk(exit_along, exit_cross))
    else:
        pts.append(mk(cur_along, exit_cross))
        pts.append(exit_)
    return pts


def _n_rows_max(cell: float, pitch: float, margin: float) -> int:
    cross_span = cell - 2.0 * margin
    if cross_span <= 0.0 or pitch <= 0.0:
        return 0
    return max(1, int(cross_span // pitch) + 1)


def _polyline_length(points: list[tuple[float, float]]) -> float:
    return sum(
        math.hypot(points[i + 1][0] - points[i][0], points[i + 1][1] - points[i][1])
        for i in range(len(points) - 1)
    )


# ---------------------------------------------------------------------------
# 조립 — 순회 경로 전체를 실제 waypoints로 바꾸고, 목표 길이(l)에 맞춰 셀별로 몇 줄
# 접을지 정한다.
# ---------------------------------------------------------------------------

def _assemble_route(
    tour: list[tuple[int, int]], cell_center: dict[tuple[int, int], tuple[float, float]],
    cell: float, pitch: float, margin: float,
    port1: tuple[float, float], port2: tuple[float, float], target_l: float,
) -> list[tuple[float, float]] | None:
    n = len(tour)
    entry_sides = [""] * n
    exit_sides = [""] * n
    entry_sides[0] = _nearest_side(cell_center[tour[0]], port1)
    for t in range(1, n):
        d = _dir_between(tour[t - 1], tour[t])
        entry_sides[t] = _OPPOSITE[d]
        exit_sides[t - 1] = d
    exit_sides[-1] = _nearest_side(cell_center[tour[-1]], port2)

    # idx 점프를 잇는 다리가 이미 방문한 셀을 다시 지나갈 수 있다 — 재방문은 접지 않는다
    # (같은 자리에 물리적으로 두 번 겹쳐 접으면 자기교차가 생긴다). "처음 방문"만 folding
    # 예산을 받는다.
    is_first_visit = [False] * n
    seen: set[tuple[int, int]] = set()
    for t, cxy in enumerate(tour):
        if cxy not in seen:
            is_first_visit[t] = True
            seen.add(cxy)

    n_max = _n_rows_max(cell, pitch, margin)
    cache: dict[tuple[int, int], list[tuple[float, float]]] = {}

    def comb(t: int, k: int) -> list[tuple[float, float]]:
        key = (t, k)
        if key not in cache:
            cx, cy = cell_center[tour[t]]
            cache[key] = _cell_meander(cx, cy, cell, pitch, margin, entry_sides[t], exit_sides[t], k)
        return cache[key]

    base_len = [0.0] * n
    max_len = [0.0] * n
    for t in range(n):
        base_len[t] = _polyline_length(comb(t, 0))
        max_len[t] = _polyline_length(comb(t, n_max)) if is_first_visit[t] else base_len[t]

    lead_in = math.hypot(port1[0] - comb(0, 0)[0][0], port1[1] - comb(0, 0)[0][1])
    lead_out = math.hypot(comb(n - 1, 0)[-1][0] - port2[0], comb(n - 1, 0)[-1][1] - port2[1])

    total_base = lead_in + lead_out + sum(base_len)
    total_max = lead_in + lead_out + sum(max_len)

    n_choice = [0] * n
    if total_max <= target_l + 1e-6:
        # 최대로 접어도 목표 미달 — 전부 최대로 접고 미달은 route_report()의 길이오차로
        # 정직하게 드러낸다(5절 원칙: 모자라면 그대로 보고, 세그먼트 밖으로 나가지 않는다).
        for t in range(n):
            n_choice[t] = n_max if is_first_visit[t] else 0
    else:
        # 순회 순서대로, 남은 deficit이 더 필요로 하는 만큼만 접는다 — 목표에 도달한 뒤의
        # 셀은 그냥 통과만 시킨다(요청 3절: "총 수용 용량이 l보다 크면 일부 셀에서 덜 접는다").
        remaining = target_l - total_base
        for t in range(n):
            if not is_first_visit[t] or remaining <= 1e-6:
                continue
            best_k, best_added = 0, 0.0
            for k in range(1, n_max + 1):
                added = _polyline_length(comb(t, k)) - base_len[t]
                if added <= remaining + 1e-6:
                    best_k, best_added = k, added
                else:
                    break
            n_choice[t] = best_k
            remaining -= best_added

    wp: list[tuple[float, float]] = [port1]
    for t in range(n):
        pts = comb(t, n_choice[t])
        wp.extend(pts if t == 0 else pts[1:])
    wp.append(port2)
    return _dedupe_close(wp)


def _dedupe_close(points: list[tuple[float, float]], eps: float = 1e-6) -> list[tuple[float, float]]:
    out = [points[0]]
    for p in points[1:]:
        if math.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) > eps:
            out.append(p)
    return out


# ---------------------------------------------------------------------------
# 검증 — 교차 개수 / 셀 이탈 길이 / 셀 이용률 (7절 보고용, main.py나 검증 스크립트가 호출)
# ---------------------------------------------------------------------------

# 실제 배선 leg(꺾은선의 각 구간)끼리 교차하는 쌍을 비교한다. 정점(포트) 공유 쌍도 빼지
# 않고 전부 비교한다(core/floorplan.py의 assign_ports 모듈 주석에 적힌 전례와 같은 원칙).
#
# 전수(O(전체 leg 수²)) 대신 "리드 leg(각 커플러의 첫/마지막 leg) vs 전체"만 비교한다 —
# v2에서 이게 완전하다: 리드가 아닌 내부 leg(_cell_meander)는 항상 자기 셀 사각형
# 안에서만 움직이도록 설계돼 있고(proto_v2.py의 exhaustive containment sweep으로 검증),
# GP가 셀을 커플러당 정확히 하나씩만 배정하므로(core/globalplacement.py) 서로 다른
# 커플러의 두 사각형은 겹치지 않는다 — 그래서 "리드가 아닌 leg 대 리드가 아닌 leg" 쌍은
# 서로 다른 커플러 사이에서 수학적으로 교차할 수 없다(둘 다 닫힌 정사각형 내부에
# 갇혀 있고, 그 정사각형끼리 겹치지 않으므로). v2는 셀 전부를 도는 전수 순회라 커플러당
# leg 수가 v1보다 훨씬 많아져서(대표 커플러 기준 4개 -> 86개) 전수 비교가 실측으로
# O(분) 단위까지 느려졌다 — 이 최적화로 O(리드 수 × 전체 leg 수)로 줄인다(리드는 커플러당
# 정확히 2개뿐).
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


def _point_in_any_cell(
    pt: tuple[float, float], cell_center: dict[tuple[int, int], tuple[float, float]], half: float,
) -> bool:
    x, y = pt
    for cx, cy in cell_center.values():
        if abs(x - cx) <= half + 1e-6 and abs(y - cy) <= half + 1e-6:
            return True
    return False


# 각 커플러 경로가 자기 세그먼트 셀 합집합 밖으로 나간 총 길이(um)의 근사값(leg마다
# samples_per_leg+1개 표본점으로 셀 포함 여부를 검사). v2에서는 순회 경로 자체(모든 셀을
# 지나는 comb 구간)는 셀 안에서만 움직이도록 설계돼 있으므로 항상 0이어야 하고, 이 값이
# 실제로 잡아내는 건 주로 포트->첫/마지막 셀 리드 구간(own_cells 필터링 대상이 아닌 직선
# 보간)이다.
def containment_violation_length(state: ChipState, samples_per_leg: int = 10) -> float:
    total = 0.0
    for c in state.couplers.values():
        if not c.waypoints or not c.segments:
            continue
        half = c.segment_size_um / 2.0
        _idx_cell, cell_center = _build_idx_cell(c.segments, c.segment_size_um)
        wp = c.waypoints
        for i in range(len(wp) - 1):
            x1, y1 = wp[i]
            x2, y2 = wp[i + 1]
            seg_len = math.hypot(x2 - x1, y2 - y1)
            if seg_len < 1e-9:
                continue
            outside = 0
            for s in range(samples_per_leg + 1):
                t = s / samples_per_leg
                pt = (x1 + (x2 - x1) * t, y1 + (y2 - y1) * t)
                if not _point_in_any_cell(pt, cell_center, half):
                    outside += 1
            total += seg_len * (outside / (samples_per_leg + 1))
    return total


# 커플러 하나의 own_cells 중 순회 경로가 실제로 방문한 셀의 비율 — v1의 핵심 손실원(A,
# 꺾임-최소 경로가 own_cells의 평균 48%만 지나가던 문제, docs/20260920_rt_segment_breakdown.md
# 발견 1)이 v2(전수 순회)에서 해소됐는지 직접 보여주는 지표(요청 4절). tour 자체를 다시
# 들고 있지 않으므로(route_report가 Coupler만 받음) waypoints가 실제로 어느 셀을
# 지나는지 좌표로 역산한다 — "tour가 이론상 전부 방문하니 항상 1.0"이라고 가정하지 않고
# 실측한다(5절 원칙과 같은 태도: 가정이 아니라 측정).
def cell_utilization(coupler: Coupler) -> float | None:
    if not coupler.segments:
        return None
    if not coupler.waypoints:
        return 0.0
    _idx_cell, cell_center = _build_idx_cell(coupler.segments, coupler.segment_size_um)
    half = coupler.segment_size_um / 2.0
    visited: set[tuple[int, int]] = set()
    for x, y in coupler.waypoints:
        for cxy, (cx, cy) in cell_center.items():
            if abs(x - cx) <= half + 1e-6 and abs(y - cy) <= half + 1e-6:
                visited.add(cxy)
    return len(visited) / len(cell_center)


# 칩 하나의 라우팅 결과를 요약한다 — 6절/7절 보고용. 실행 시간은 여기 포함하지 않는다 —
# main.py의 run_stage()가 이미 스테이지별 칩별 시간을 재고 있다.
def route_report(state: ChipState) -> dict:
    couplers = state.couplers
    with_segments = [c for c in couplers.values() if c.segments]
    routed = [c for c in with_segments if c.waypoints]
    turns = [len(c.waypoints) - 2 for c in routed]
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
        "turns_avg": (sum(turns) / len(turns)) if turns else 0.0,
        "turns_max": max(turns) if turns else 0,
        "len_err_avg_pct": (sum(length_err_pct) / len(length_err_pct)) if length_err_pct else 0.0,
        "len_err_max_abs_pct": max((abs(e) for e in length_err_pct), default=0.0),
        "within_5pct_frac": (sum(1 for e in length_err_pct if abs(e) <= 5.0) / len(length_err_pct)) if length_err_pct else 0.0,
        "cell_utilization_avg": (sum(utils) / len(utils)) if utils else 0.0,
        "oob_length_um": containment_violation_length(state),
    }
