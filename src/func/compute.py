import math

from core.state import ChipState, Qubit
from func.metric import Metric

# QPlacer(Zhang et al., ISCA'25) Eq.15와 동일한 디튜닝 임계값(GHz).
# |f_i - f_j| <= 이 값이면 두 큐빗을 "근접 주파수 쌍(crosstalk 후보)"로 본다.
DEFAULT_DETUNING_GHZ = 0.1


# ─── 기하 헬퍼 ────────────────────────────────────────────────────────────
# 큐빗의 경계 상자를 (x0, x1, y0, y1)로 반환
def _qubit_aabb(q: Qubit) -> tuple[float, float, float, float]:
    return (q.x - q.w / 2.0, q.x + q.w / 2.0, q.y - q.h / 2.0, q.y + q.h / 2.0)


# 커플러 양 끝 큐빗의 현재 좌표로 Coupler.boundary를 계산 (끝점이 없으면 None)
def _coupler_aabb(state: ChipState, coupler) -> tuple[float, float, float, float] | None:
    q1 = state.qubits.get(coupler.q1)
    q2 = state.qubits.get(coupler.q2)
    if q1 is None or q2 is None:
        return None
    return coupler.boundary(q1, q2)


# 두 경계 상자가 겹치는지 여부
def _aabb_overlap(a: tuple, b: tuple) -> bool:
    ax0, ax1, ay0, ay1 = a
    bx0, bx1, by0, by1 = b
    return ax0 < bx1 and bx0 < ax1 and ay0 < by1 and by0 < ay1


# 두 선분 (p1,p2), (p3,p4)이 교차하는지 여부 (공선 접촉 포함, 표준 오리엔테이션 판정)
def _segments_cross(p1, p2, p3, p4) -> bool:
    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def on_segment(a, b, c):
        return min(a[0], b[0]) <= c[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= c[1] <= max(a[1], b[1])

    d1, d2 = orient(p3, p4, p1), orient(p3, p4, p2)
    d3, d4 = orient(p1, p2, p3), orient(p1, p2, p4)

    if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and d1 != 0 and d2 != 0 and d3 != 0 and d4 != 0:
        return True
    if d1 == 0 and on_segment(p3, p4, p1):
        return True
    if d2 == 0 and on_segment(p3, p4, p2):
        return True
    if d3 == 0 and on_segment(p1, p2, p3):
        return True
    if d4 == 0 and on_segment(p1, p2, p4):
        return True
    return False


def _active_qubit_indices(state: ChipState, active: set[int] | None) -> list[int]:
    active = active if active is not None else state.active_qubits
    return [i for i in state.qubits if active is None or i in active]


# ─── 개별 metric 계산 ─────────────────────────────────────────────────────

# QPlacer Eq.15 형태의 hotspot proportion을 계산.
#   P_h = Σ_{i,j} (overlap_x + overlap_y) · centroid_dist · τ(Δf) / Σ_qubit area
# 직접 커플러로 이어진 쌍(의도된 결합)은 crosstalk 대상에서 제외한다.
# 반환: (proportion, 겹침에 관여한 큐빗 인덱스 집합)
def compute_freq_hotspot_proportion(
    state: ChipState,
    active: set[int] | None = None,
    detuning_ghz: float = DEFAULT_DETUNING_GHZ,
) -> tuple[float, set[int]]:
    idxs = _active_qubit_indices(state, active)
    if len(idxs) < 2:
        return 0.0, set()

    qubits = state.qubits
    a_poly = sum(qubits[i].w * qubits[i].h for i in idxs)
    if a_poly <= 0.0:
        return 0.0, set()

    coupled = {frozenset(e) for e in state.cmap}

    numerator = 0.0
    hotspot: set[int] = set()
    for a in range(len(idxs)):
        i = idxs[a]
        qi = qubits[i]
        for b in range(a + 1, len(idxs)):
            j = idxs[b]
            if frozenset((i, j)) in coupled:
                continue
            qj = qubits[j]
            if abs(qi.f - qj.f) > detuning_ghz:
                continue
            x0i, x1i, y0i, y1i = _qubit_aabb(qi)
            x0j, x1j, y0j, y1j = _qubit_aabb(qj)
            overlap_x = min(x1i, x1j) - max(x0i, x0j)
            overlap_y = min(y1i, y1j) - max(y0i, y0j)
            if overlap_x <= 0.0 or overlap_y <= 0.0:
                continue
            dist = math.hypot(qi.x - qj.x, qi.y - qj.y)
            numerator += (overlap_x + overlap_y) * dist
            hotspot.add(i)
            hotspot.add(j)

    return numerator / a_poly, hotspot


# 주파수 hotspot에 관여한 큐빗 수 (compute_freq_hotspot_proportion과 같은 정의를 공유)
def compute_num_hotspot_qubits(
    state: ChipState,
    active: set[int] | None = None,
    detuning_ghz: float = DEFAULT_DETUNING_GHZ,
) -> int:
    _, hotspot = compute_freq_hotspot_proportion(state, active=active, detuning_ghz=detuning_ghz)
    return len(hotspot)


# 칩 면적 대비 배치 효율: (칩 면적 - 모든 큐빗을 감싸는 최소 경계상자 면적) / 칩 면적
def compute_area_efficiency(state: ChipState) -> float:
    qubits = state.qubits
    chip_area = state.chip_width * state.chip_height
    if not qubits or chip_area <= 0.0:
        return 0.0

    boxes = [_qubit_aabb(q) for q in qubits.values()]
    x0 = min(b[0] for b in boxes)
    x1 = max(b[1] for b in boxes)
    y0 = min(b[2] for b in boxes)
    y1 = max(b[3] for b in boxes)
    min_area = (x1 - x0) * (y1 - y0)

    return max(0.0, (chip_area - min_area) / chip_area)


# 커플링맵 논리 엣지(큐빗 중심을 잇는 직선 근사)끼리 교차하는 쌍의 개수.
# 끝점을 공유하는 엣지(같은 큐빗에 연결된 두 엣지)는 교차로 치지 않는다.
def compute_edge_cross_point(state: ChipState) -> int:
    qubits = state.qubits
    edges = [(u, v) for u, v in state.cmap if u in qubits and v in qubits]

    count = 0
    for a in range(len(edges)):
        u1, v1 = edges[a]
        p1, p2 = (qubits[u1].x, qubits[u1].y), (qubits[v1].x, qubits[v1].y)
        for b in range(a + 1, len(edges)):
            u2, v2 = edges[b]
            if u1 in (u2, v2) or v1 in (u2, v2):
                continue
            p3, p4 = (qubits[u2].x, qubits[u2].y), (qubits[v2].x, qubits[v2].y)
            if _segments_cross(p1, p2, p3, p4):
                count += 1
    return count


# 커플러 경계 상자(AABB)끼리 겹치는 쌍의 개수
def compute_coupler_cross_point(state: ChipState) -> int:
    boxes = [box for c in state.couplers.values() if (box := _coupler_aabb(state, c)) is not None]

    count = 0
    for a in range(len(boxes)):
        for b in range(a + 1, len(boxes)):
            if _aabb_overlap(boxes[a], boxes[b]):
                count += 1
    return count


# DRC 4종(큐빗-큐빗/큐빗-커플러/커플러-커플러 겹침, 논리 엣지 교차)을 한 번에 계산.
# edge_cross_point/coupler_cross_point를 이미 계산해뒀다면 넘겨서 재계산을 피할 수 있다.
def compute_drc(
    state: ChipState,
    edge_cross_point: int | None = None,
    coupler_cross_point: int | None = None,
) -> dict[str, bool]:
    qubits = list(state.qubits.items())

    qq_overlap = 0
    for a in range(len(qubits)):
        box_a = _qubit_aabb(qubits[a][1])
        for b in range(a + 1, len(qubits)):
            if _aabb_overlap(box_a, _qubit_aabb(qubits[b][1])):
                qq_overlap += 1

    # 커플러는 자기 자신의 두 큐빗에 물리적으로 붙어야 하므로(패드 연결), 그 겹침은
    # 설계상 정상 — DRC 위반 신호에서 자기 자신의 두 큐빗은 제외한다.
    qc_overlap = 0
    for coupler in state.couplers.values():
        cbox = _coupler_aabb(state, coupler)
        if cbox is None:
            continue
        for i, q in state.qubits.items():
            if i in (coupler.q1, coupler.q2):
                continue
            if _aabb_overlap(cbox, _qubit_aabb(q)):
                qc_overlap += 1

    if edge_cross_point is None:
        edge_cross_point = compute_edge_cross_point(state)
    if coupler_cross_point is None:
        coupler_cross_point = compute_coupler_cross_point(state)

    return {
        "drc_qq_overlap": qq_overlap == 0,
        "drc_qc_overlap": qc_overlap == 0,
        "drc_cc_overlap": coupler_cross_point == 0,
        "drc_edge_cross": edge_cross_point == 0,
    }


# TODO: fidelity 계산식 미정 (01_mainref에도 대응 지표가 없음).
# 공식이 정해지기 전까지는 이 함수를 직접 호출하지 않도록 한다 — compute_metrics()도 호출하지 않는다.
def compute_fidelity(_state: ChipState) -> float:
    raise NotImplementedError("compute_fidelity: 계산식이 아직 정의되지 않았습니다 (TODO)")


# 배치된 칩 상태(ChipState)로부터 Metric 전체를 계산해 조립.
# fidelity는 계산식 미정이라 None으로 남긴다(compute_fidelity 참고).
def compute_metrics(
    state: ChipState,
    active: set[int] | None = None,
    detuning_ghz: float = DEFAULT_DETUNING_GHZ,
) -> Metric:
    freq_hotspot_proportion, hotspot_qubits = compute_freq_hotspot_proportion(
        state, active=active, detuning_ghz=detuning_ghz,
    )
    edge_cross_point = compute_edge_cross_point(state)
    coupler_cross_point = compute_coupler_cross_point(state)
    drc = compute_drc(state, edge_cross_point=edge_cross_point, coupler_cross_point=coupler_cross_point)

    return Metric(
        freq_hotspot_proportion=freq_hotspot_proportion,
        area_utilization=compute_area_efficiency(state),
        edge_cross_point=edge_cross_point,
        coupler_cross_point=coupler_cross_point,
        num_hotspot_qubits=len(hotspot_qubits),
        drc_qq_overlap=drc["drc_qq_overlap"],
        drc_qc_overlap=drc["drc_qc_overlap"],
        drc_cc_overlap=drc["drc_cc_overlap"],
        drc_edge_cross=drc["drc_edge_cross"],
    )
