import math

from core.state import AABB_EPS_UM, ChipState, Coupler, Qubit, Segment, coupler_own_port_cell
from func.metric import Metric

# QPlacer(Zhang et al., ISCA'25) Eq.15와 동일한 디튜닝 임계값(GHz).
# |f_i - f_j| <= 이 값이면 두 큐빗을 "근접 주파수 쌍(crosstalk 후보)"로 본다.
DEFAULT_DETUNING_GHZ = 0.1


# ─── 기하 헬퍼 ────────────────────────────────────────────────────────────
# 큐빗의 경계 상자를 (x0, x1, y0, y1)로 반환
def _qubit_aabb(q: Qubit) -> tuple[float, float, float, float]:
    return (q.x - q.w / 2.0, q.x + q.w / 2.0, q.y - q.h / 2.0, q.y + q.h / 2.0)


# 세그먼트가 필요한데(num_segments > 0) 하나도 못 놓은(segments가 빈) 커플러 수 — GP가
# 실패 처리한 커플러(core/globalplacement.py의 SegmentPlacementInfeasibleError 대상은
# 아니고, 개별 실패). num_segments == 0(애초에 세그먼트가 필요 없는 짧은 커플러)은 정상
# 완료이지 미배치가 아니므로 여기서 제외한다. 이 값이 0이 아니면 qc/cc 위반 "0"은 검사가
# 그만큼 스킵됐다는 뜻이지 실제로 통과했다는 뜻이 아니다 — 반드시 같이 읽어야 한다.
def _num_unplaced_couplers(state: ChipState) -> int:
    return sum(1 for c in state.couplers.values() if c.num_segments > 0 and not c.segments)


# 세그먼트 하나의 경계 상자
def _segment_aabb(coupler: Coupler, seg: Segment) -> tuple[float, float, float, float]:
    half = coupler.segment_size_um / 2.0
    return (seg.x - half, seg.x + half, seg.y - half, seg.y + half)


# (coupler, segment) 쌍 전체 — qc/cc 둘 다 이 목록 위에서 개별 세그먼트 기준으로 판정한다.
def _all_segments(state: ChipState) -> list[tuple[Coupler, Segment]]:
    return [(c, s) for c in state.couplers.values() for s in c.segments]


# 두 경계 상자가 겹치는지 여부. AABB_EPS_UM(core/state.py)만큼의 부동소수 허용 오차를
# 둔다 — 물리적 clearance가 아니라 순수 반올림 잡음 흡수용(왜 필요한지는 그 상수 정의
# 주석 참고).
def _aabb_overlap(a: tuple, b: tuple) -> bool:
    ax0, ax1, ay0, ay1 = a
    bx0, bx1, by0, by1 = b
    return ax0 < bx1 - AABB_EPS_UM and bx0 < ax1 - AABB_EPS_UM \
        and ay0 < by1 - AABB_EPS_UM and by0 < ay1 - AABB_EPS_UM


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


# active=None이면 전체 큐빗이 활성이다(필터 없음). 예전엔 여기서 ChipState.active_qubits로
# 대체 폴백을 시도했는데, 그 필드가 ChipState(core/state.py)에 존재한 적이 없어 compute_metrics
# 호출 즉시 AttributeError로 죽었다(active를 안 넘기는 모든 호출부 — Evaluation.compute 등 —
# 가 여기 걸림). 2026-09-17 기준 이 저장소 어디에도(입력 JSON, params.json, ChipState) 결함
# 큐빗/수율 마스킹 개념이 실제로 존재하지 않는다 — params.json의 chip_width/height 설명이
# 언급하는 src/legacy_rl/tests/test_defect_mask.py도 이 저장소엔 없는 파일이라(다른 스냅샷의
# 낡은 설명으로 보임) "쓸 계획이 있다"는 근거가 못 된다. 그래서 존재하지 않는 필드를 참조하는
# 폴백을 지웠다 — active 매개변수 자체(호출부가 명시적으로 부분집합을 넘기는 경로)는 멀쩡히
# 동작하므로 그대로 둔다. 결함 마스킹이 실제로 필요해지면 그때 ChipState에 필드를 추가하고
# 여기서 다시 폴백하면 된다.
def _active_qubit_indices(state: ChipState, active: set[int] | None) -> list[int]:
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


# 서로 다른 커플러에 속한 세그먼트끼리 실제로 겹치는 쌍의 개수. 2026-09-17 이전엔
# Coupler.bbox()(그 커플러의 세그먼트 전체를 감싸는 사각형)끼리 비교했는데, bbox는 대리
# 지표일 뿐이다 — 세그먼트가 장애물을 피해 이쪽저쪽에 놓이면 실제로는 안 겹쳐도 그 전체를
# 감싸는 bbox끼리는 겹칠 수 있다(01_mainref에서 이 대리 지표 때문에 위반의 66.6%가
# false positive였다 — docs/20260917_gp_qc_findings.md 참고). 이제 개별 세그먼트
# AABB로 직접 비교한다 — GP가 격자로 세그먼트 간 비중첩을 이미 보장하므로(core/
# globalplacement.py) 정상 동작 중엔 0이 나와야 하고, 0이 아니면 그 자체가 실제 버그 신호다.
def compute_coupler_cross_point(state: ChipState) -> int:
    segs = _all_segments(state)
    count = 0
    for a in range(len(segs)):
        ca, sa = segs[a]
        box_a = _segment_aabb(ca, sa)
        for b in range(a + 1, len(segs)):
            cb, sb = segs[b]
            if ca.id == cb.id:
                continue  # 같은 커플러 내부는 cc(커플러 간) 대상이 아니다
            if _aabb_overlap(box_a, _segment_aabb(cb, sb)):
                count += 1
    return count


# DRC 4종(큐빗-큐빗/큐빗-커플러/커플러-커플러 겹침, 논리 엣지 교차)을 한 번에 계산.
# edge_cross_point/coupler_cross_point를 이미 계산해뒀다면 넘겨서 재계산을 피할 수 있다.
# num_unplaced_couplers: 세그먼트가 필요한데 못 놓은 커플러 수. 이 값이 0이 아니면
# drc_qc_overlap/drc_cc_overlap의 "위반 0"은 실제 통과가 아니라 그만큼 검사가 스킵됐다는
# 뜻이므로 반드시 함께 읽어야 한다.
def compute_drc(
    state: ChipState,
    edge_cross_point: int | None = None,
    coupler_cross_point: int | None = None,
) -> dict[str, bool | int]:
    qubits = list(state.qubits.items())

    qq_overlap = 0
    for a in range(len(qubits)):
        box_a = _qubit_aabb(qubits[a][1])
        for b in range(a + 1, len(qubits)):
            if _aabb_overlap(box_a, _qubit_aabb(qubits[b][1])):
                qq_overlap += 1

    num_unplaced_couplers = _num_unplaced_couplers(state)

    # 큐빗-세그먼트 겹침(자기 큐빗의 배정 포트 근처만 제외 — core.state.coupler_own_port_cell()
    # 로 통일. core/globalplacement.py의 GP 장애물 판정도 같은 함수를 쓴다, 두 곳에 따로
    # 구현하면 또 어긋난다). 2026-09-17 이전엔 bbox() 기준이었다 — compute_coupler_cross_point
    # 위 주석과 같은 이유로 개별 세그먼트 기준으로 바꿨다: bbox가 큐빗 위를 가로지르면 실제
    # 세그먼트는 안 건드려도 위반으로 잡히는 헛경보가 났다. 2026-09-20: 자기 큐빗이면 몸체
    # 전체를 봐주던 예외를 "배정 포트 근처 셀만"으로 좁혔다 — 실측 결과 세그먼트가 포트
    # 근처가 아니라 큐빗 중심 코앞(반폭의 최대 91%)까지 파고들었다
    # (docs/20260920_segment_model_review.md 발견 3). 이제 그 정도의 몸체 내부 겹침은
    # qc_overlap 위반으로 잡힌다.
    qc_overlap = 0
    for coupler, seg in _all_segments(state):
        sbox = _segment_aabb(coupler, seg)
        assignment = state.port_assignment.get((coupler.q1, coupler.q2))
        for i, q in state.qubits.items():
            if not _aabb_overlap(sbox, _qubit_aabb(q)):
                continue
            if coupler_own_port_cell(coupler, i, q, assignment, sbox):
                continue
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
        "num_unplaced_couplers": num_unplaced_couplers,
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
        num_unplaced_couplers=drc["num_unplaced_couplers"],
    )
