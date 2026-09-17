"""Legalization: GP 결과의 DRC 위반(qq/qc/cc/die 밖)을 최소 변위로 해소.

결정할 것은 큐빗/세그먼트 좌표뿐이다(포트 배정, 커플러 박스, 세그먼트 개수는 FP/GP가
이미 확정했고 LG는 건드리지 않는다). 하드 제약은 두 가지: (1) DRC를 만족할 것,
(2) FP가 확정한 평면 임베딩이 깨지지 않을 것(큐빗을 옮긴 뒤 proper crossing이 늘면 그
이동은 버린다 — 01_mainref의 BFS DRC 마스킹과 같은 원칙: 페널티를 주는 게 아니라 후보
자체를 배제한다). 목적은 변위 최소화 — 위반을 없애는 데 필요한 최소한만 움직인다.

DRC 정의는 func/compute.py의 compute_drc를 따르되(2026-09-17 결정 사항):
  - qc(큐빗-세그먼트)는 개별 세그먼트 AABB 기준이다. bbox()(세그먼트 전체를 감싸는 사각형)
    를 쓰면 실제로는 안 겹쳐도 우회 경로 때문에 bbox가 큐빗 위를 가로질러 헛경보가 난다
    (01_mainref에서 qc 위반의 66.6%가 이런 false positive였다 —
    docs/20260917_gp_qc_findings.md 참고). cc(커플러 간 세그먼트)도 같은 이유로 개별
    세그먼트 기준이다.
  - clearance는 zero-margin이다. 01_mainref의 qc_min_gap_um=500(QPlacer d_q=400+d_r=100)
    은 그 원문이 스스로 "aggressively set"이라고 밝힌 값이라 도입하지 않는다 — 물리적
    근거가 따로 생기면 그때 추가한다.
  - 자기 큐빗 예외는 Coupler.is_own_qubit()로 판정한다(GP의 장애물 판정, compute_drc의
    qc 판정과 동일한 정의 — 세 곳이 각자 구현하면 또 어긋난다).

GP에서 세그먼트를 하나도 못 놓은 커플러(segments=[])는 LG가 다시 배치하지 않는다 — 세그먼트
배치는 GP의 일이고, LG는 "이미 놓인 것들의 위반을 최소 변위로 없애는" 역할로 한정한다.
그런 커플러는 func/compute.py의 num_unplaced_couplers로 계속 보고된다.

2026-09-17 실측: GP 산출물(6칩)에서 qq/qc/cc/die밖 위반이 전부 0이었다(GP가 격자 기반
순차 배치로 이미 비중첩을 보장하기 때문). 그래서 이 모듈은 현재 벤치마크에 대해 사실상
no-op이다 — 위반이 없으면 즉시 원본 state를 그대로 반환한다. 그렇다고 이 모듈이 불필요한
건 아니다: GP가 실패하거나 입력이 달라지면 위반이 생길 수 있고, 그때 실제로 동작해야 한다
(합성 위반으로 이 로직 자체는 별도 검증했다).
"""
import logging
import math
from dataclasses import replace

from core.state import AABB_EPS_UM, ChipState, Coupler, Qubit, Segment
from core.floorplan import count_crossings, sort_edges
from core.globalplacement import (
    _cell_blocked, _die_max_index, _qubit_owner_cells, _region_index_range,
)


# 위반을 최소 변위로 못 없앤 커플러/큐빗이 있어도(부분 실패) 칩 전체를 버리지 않는다 —
# Floorplan/GlobalPlacement와 같은 패턴. 이 예외는 지금 어디서도 올라오지 않는다(부분
# 실패는 unresolved 카운트로 로그만 남기고 계속 진행) — 향후 "위반이 너무 많아 이 칩
# 자체가 legalize 불가능"으로 판단할 기준이 생기면 여기서 쓰면 된다.
class LegalizationInfeasibleError(Exception):
    pass


class Legalization:
    def __init__(self, params):
        self.params = params
        self.skipped: list[tuple[str, str]] = []

    def run(self, states: list[ChipState]) -> list[ChipState]:
        out: list[ChipState] = []
        for state in states:
            try:
                out.append(self._place(state))
            except LegalizationInfeasibleError as e:
                self.skipped.append((state.processor_name, str(e)))
        return out

    def _place(self, state: ChipState) -> ChipState:
        cell = float(self.params.segment_size_um)
        chip_w, chip_h = state.chip_width, state.chip_height

        oob_q = _find_qubit_oob(state.qubits, chip_w, chip_h)
        qq = _find_qq_overlaps(state.qubits)
        oob_s = _find_seg_oob(state.couplers, chip_w, chip_h)
        qc = _find_qc_overlaps(state.qubits, state.couplers)
        cc = _find_cc_overlaps(state.couplers)
        n_unplaced = sum(1 for c in state.couplers.values() if c.num_segments > 0 and not c.segments)

        if not (oob_q or qq or oob_s or qc or cc):
            logging.info(
                "[LG] %s: DRC 위반 0건(qq=0 qc=0 cc=0 die밖=0) — no-op (num_unplaced_couplers=%d)",
                state.processor_name, n_unplaced,
            )
            return state

        edges = sort_edges(state.cmap)
        qubits = dict(state.qubits)
        couplers = dict(state.couplers)

        qubit_owner = _qubit_owner_cells(qubits, cell, chip_w, chip_h)
        segment_occupied: set[tuple[int, int]] = set()
        for c in couplers.values():
            for s in c.segments:
                segment_occupied.add(_seg_cell(s, cell))
        i_max = _die_max_index(chip_w, cell)
        j_max = _die_max_index(chip_h, cell)

        qubit_disp: list[float] = []
        segment_disp: list[float] = []
        unresolved = 0

        # --- 큐빗 위반: die 밖 먼저, 그다음 qq. 둘 다 이동 후보를 만들고, proper crossing이
        # 늘면(embedding 훼손) 후보를 버린다 — 페널티가 아니라 후보 배제.
        for qid in oob_q:
            cand = _resolve_qubit_oob(qubits[qid], chip_w, chip_h)
            if cand is not None and _crossing_safe(qubits, edges, {qid: cand}):
                old = qubits[qid]
                qubits[qid] = replace(old, x=cand[0], y=cand[1])
                qubit_disp.append(math.hypot(cand[0] - old.x, cand[1] - old.y))
            else:
                unresolved += 1

        for qa_id, qb_id in qq:
            qa, qb = qubits[qa_id], qubits[qb_id]
            if not _aabb_overlap(_qubit_aabb(qa), _qubit_aabb(qb)):
                continue  # die-oob 수정 등으로 이미 해소됨
            cand = _resolve_qq(qa, qb)
            if cand is not None and _crossing_safe(qubits, edges, {qa_id: cand[0], qb_id: cand[1]}):
                (ax, ay), (bx, by) = cand
                qubit_disp.append(math.hypot(ax - qa.x, ay - qa.y))
                qubit_disp.append(math.hypot(bx - qb.x, by - qb.y))
                qubits[qa_id] = replace(qa, x=ax, y=ay)
                qubits[qb_id] = replace(qb, x=bx, y=by)
            else:
                unresolved += 1

        # --- 세그먼트 위반: die 밖 -> qc -> cc. 자기 박스(coupler_regions) 안에서 현재
        # 위치와 가장 가까운 빈 셀로만 옮긴다 — 세그먼트는 임베딩과 무관하므로 crossing
        # 재검사가 필요 없다.
        def move_segment(key: tuple[int, int], seg_idx: int) -> float | None:
            box = state.coupler_regions.get(key)
            if box is None:
                return None
            coupler = couplers[key]
            seg = next((s for s in coupler.segments if s.idx == seg_idx), None)
            if seg is None:
                return None  # 이미 다른 위반 처리로 옮겨졌거나 없음
            old_cell = _seg_cell(seg, cell)
            segment_occupied.discard(old_cell)
            new_cell = _nearest_free_cell_in_box(
                (seg.x, seg.y), box, cell, coupler, qubit_owner, segment_occupied, i_max, j_max,
            )
            if new_cell is None:
                segment_occupied.add(old_cell)
                return None
            nx, ny = (new_cell[0] + 0.5) * cell, (new_cell[1] + 0.5) * cell
            segment_occupied.add(new_cell)
            new_segments = [replace(s, x=nx, y=ny) if s.idx == seg_idx else s for s in coupler.segments]
            couplers[key] = replace(coupler, segments=new_segments)
            return math.hypot(nx - seg.x, ny - seg.y)

        for key, seg_idx in oob_s:
            d = move_segment(key, seg_idx)
            if d is None:
                unresolved += 1
            else:
                segment_disp.append(d)

        for key, seg_idx, _qid in qc:
            d = move_segment(key, seg_idx)
            if d is None:
                unresolved += 1
            else:
                segment_disp.append(d)

        for _ka, _sa_idx, kb, sb_idx in cc:
            d = move_segment(kb, sb_idx)
            if d is None:
                unresolved += 1
            else:
                segment_disp.append(d)

        logging.info(
            "[LG] %s: 위반(수정전) qubit_oob=%d qq=%d seg_oob=%d qc=%d cc=%d | 못고침=%d | "
            "num_unplaced_couplers=%d | 큐빗변위 avg/max=%.2f/%.2f um | "
            "세그먼트변위 avg/max=%.2f/%.2f um",
            state.processor_name, len(oob_q), len(qq), len(oob_s), len(qc), len(cc), unresolved,
            n_unplaced,
            (sum(qubit_disp) / len(qubit_disp)) if qubit_disp else 0.0,
            max(qubit_disp) if qubit_disp else 0.0,
            (sum(segment_disp) / len(segment_disp)) if segment_disp else 0.0,
            max(segment_disp) if segment_disp else 0.0,
        )

        return replace(state, qubits=qubits, couplers=couplers)


# ---------------------------------------------------------------------------
# 기하 헬퍼 (func/compute.py와 같은 정의 — AABB 겹침은 두 곳에 따로 짜도 어긋날 여지가
# 없는 단순 기하라 그대로 재구현했다. 다만 "자기 큐빗 제외" 같은 판단은 반드시
# Coupler.is_own_qubit()을 통해서만 한다.)
# ---------------------------------------------------------------------------

def _qubit_aabb(q: Qubit) -> tuple[float, float, float, float]:
    return (q.x - q.w / 2.0, q.x + q.w / 2.0, q.y - q.h / 2.0, q.y + q.h / 2.0)


def _segment_aabb(coupler: Coupler, seg: Segment) -> tuple[float, float, float, float]:
    half = coupler.segment_size_um / 2.0
    return (seg.x - half, seg.x + half, seg.y - half, seg.y + half)


# AABB_EPS_UM(core/state.py)만큼 부동소수 허용 오차를 둔다 — func/compute.py의
# _aabb_overlap과 반드시 같은 상수/부호를 써야 한다(두 곳이 어긋나면 LG가 "legal"이라고
# 판단한 상태를 compute_drc가 "위반"이라고 다시 잡아내는 모순이 생긴다).
def _aabb_overlap(a: tuple, b: tuple) -> bool:
    ax0, ax1, ay0, ay1 = a
    bx0, bx1, by0, by1 = b
    return ax0 < bx1 - AABB_EPS_UM and bx0 < ax1 - AABB_EPS_UM \
        and ay0 < by1 - AABB_EPS_UM and by0 < ay1 - AABB_EPS_UM


def _seg_cell(seg: Segment, cell: float) -> tuple[int, int]:
    # GP는 항상 x=(i+0.5)*cell로 세그먼트를 놓으므로(core/globalplacement.py), 역산은
    # round(x/cell - 0.5)로 정확히 i를 복원한다.
    return (int(round(seg.x / cell - 0.5)), int(round(seg.y / cell - 0.5)))


# ---------------------------------------------------------------------------
# 위반 탐지
# ---------------------------------------------------------------------------

def _find_qubit_oob(qubits: dict[int, Qubit], chip_w: float, chip_h: float) -> list[int]:
    out = []
    for qid, q in qubits.items():
        x0, x1, y0, y1 = _qubit_aabb(q)
        if x0 < -1e-6 or x1 > chip_w + 1e-6 or y0 < -1e-6 or y1 > chip_h + 1e-6:
            out.append(qid)
    return sorted(out)


def _find_qq_overlaps(qubits: dict[int, Qubit]) -> list[tuple[int, int]]:
    ids = sorted(qubits)
    out = []
    for a in range(len(ids)):
        box_a = _qubit_aabb(qubits[ids[a]])
        for b in range(a + 1, len(ids)):
            if _aabb_overlap(box_a, _qubit_aabb(qubits[ids[b]])):
                out.append((ids[a], ids[b]))
    return out


def _find_seg_oob(
    couplers: dict[tuple[int, int], Coupler], chip_w: float, chip_h: float,
) -> list[tuple[tuple[int, int], int]]:
    out = []
    for key, c in couplers.items():
        for s in c.segments:
            x0, x1, y0, y1 = _segment_aabb(c, s)
            if x0 < -1e-6 or x1 > chip_w + 1e-6 or y0 < -1e-6 or y1 > chip_h + 1e-6:
                out.append((key, s.idx))
    return sorted(out)


def _find_qc_overlaps(
    qubits: dict[int, Qubit], couplers: dict[tuple[int, int], Coupler],
) -> list[tuple[tuple[int, int], int, int]]:
    out = []
    for key, c in couplers.items():
        for s in c.segments:
            sbox = _segment_aabb(c, s)
            for qid, q in qubits.items():
                if c.is_own_qubit(qid):
                    continue
                if _aabb_overlap(sbox, _qubit_aabb(q)):
                    out.append((key, s.idx, qid))
    return sorted(out)


def _find_cc_overlaps(
    couplers: dict[tuple[int, int], Coupler],
) -> list[tuple[tuple[int, int], int, tuple[int, int], int]]:
    segs = []
    for key, c in couplers.items():
        for s in c.segments:
            segs.append((key, c, s, _segment_aabb(c, s)))

    out = []
    for a in range(len(segs)):
        ka, ca, sa, boxa = segs[a]
        for b in range(a + 1, len(segs)):
            kb, cb, sb, boxb = segs[b]
            if ca.id == cb.id:
                continue
            if _aabb_overlap(boxa, boxb):
                out.append((ka, sa.idx, kb, sb.idx))
    return out


# ---------------------------------------------------------------------------
# 큐빗 이동 후보 + 임베딩 불변 검사
# ---------------------------------------------------------------------------

def _resolve_qubit_oob(q: Qubit, chip_w: float, chip_h: float) -> tuple[float, float] | None:
    x, y = q.x, q.y
    x0, x1, y0, y1 = _qubit_aabb(q)
    if x0 < 0.0:
        x += -x0
    elif x1 > chip_w:
        x += chip_w - x1
    if y0 < 0.0:
        y += -y0
    elif y1 > chip_h:
        y += chip_h - y1
    return (x, y) if (x, y) != (q.x, q.y) else None


# 두 겹치는 큐빗을 벌린다 — penetration이 더 작은 축으로 밀어야 총 변위가 최소가 된다
# (그 축만 밀어도 안 겹치게 되므로 반대 축은 건드릴 필요가 없다). 양쪽을 절반씩 반대
# 방향으로 밀어 "최소 변위"를 두 큐빗이 공평하게 나눠 갖게 한다.
def _resolve_qq(qa: Qubit, qb: Qubit) -> tuple[tuple[float, float], tuple[float, float]] | None:
    dx, dy = qa.x - qb.x, qa.y - qb.y
    overlap_x = (qa.w + qb.w) / 2.0 - abs(dx)
    overlap_y = (qa.h + qb.h) / 2.0 - abs(dy)
    if overlap_x <= 0.0 or overlap_y <= 0.0:
        return None
    eps = 1e-6
    if overlap_x <= overlap_y:
        push = overlap_x / 2.0 + eps
        sign = 1.0 if dx >= 0.0 else -1.0
        return (qa.x + sign * push, qa.y), (qb.x - sign * push, qb.y)
    push = overlap_y / 2.0 + eps
    sign = 1.0 if dy >= 0.0 else -1.0
    return (qa.x, qa.y + sign * push), (qb.x, qb.y - sign * push)


# candidate_positions에 담긴 큐빗들을 그 좌표로 옮겼을 때, 논리 엣지(state.cmap) 사이의
# proper crossing(공선/접촉 제외, X자로 가로지르는 것만)이 이동 전보다 늘어나지 않는지
# 검사한다. 01_mainref의 BFS DRC 마스킹과 같은 원칙 — 늘어나는 후보는 여기서 걸러
# 애초에 적용되지 않는다(페널티 아님).
def _crossing_safe(
    qubits: dict[int, Qubit], edges: list[tuple[int, int]],
    candidate_positions: dict[int, tuple[float, float]],
) -> bool:
    pos_before = {qid: (q.x, q.y) for qid, q in qubits.items()}
    pos_after = dict(pos_before)
    pos_after.update(candidate_positions)
    before = count_crossings(pos_before, edges, proper_only=True)
    after = count_crossings(pos_after, edges, proper_only=True)
    return after <= before


# ---------------------------------------------------------------------------
# 세그먼트 이동 후보 — 자기 박스 안에서 현재 위치와 가장 가까운 빈 셀
# ---------------------------------------------------------------------------

# GP(core/globalplacement.py)의 격자 인프라(_qubit_owner_cells/_cell_blocked/
# _region_index_range/_die_max_index)를 그대로 재사용한다 — "장애물이 뭔가"를 여기서
# 다시 정의하면 GP와 또 어긋난다. GP는 포트에서부터 체인을 훑는 순서(boustrophedon)로
# 셀을 고르지만, LG는 목적이 달라서(변위 최소화, 체인 연결성 아님) 박스 전체를 훑어
# 현재 위치와 가장 가까운 빈 셀을 고른다 — 박스가 12~24칸 정도로 작아(GlobalPlacement
# 검증 참고) 전수 탐색해도 비용이 무시할 수준이다.
def _nearest_free_cell_in_box(
    cur_xy: tuple[float, float], box: tuple[float, float, float, float], cell: float,
    coupler: Coupler, qubit_owner: dict[tuple[int, int], set[int]],
    segment_occupied: set[tuple[int, int]], i_max: int, j_max: int,
) -> tuple[int, int] | None:
    ri = _region_index_range(box[0], box[1], cell, i_max)
    rj = _region_index_range(box[2], box[3], cell, j_max)
    if ri is None or rj is None:
        return None
    i0, i1 = ri
    j0, j1 = rj
    ci = int(round(cur_xy[0] / cell - 0.5))
    cj = int(round(cur_xy[1] / cell - 0.5))

    best_d = None
    best_cell = None
    for i in range(i0, i1 + 1):
        for j in range(j0, j1 + 1):
            if _cell_blocked(coupler, (i, j), qubit_owner, segment_occupied):
                continue
            d = (i - ci) ** 2 + (j - cj) ** 2
            if best_d is None or d < best_d:
                best_d, best_cell = d, (i, j)
    return best_cell
