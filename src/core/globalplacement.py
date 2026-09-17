"""GlobalPlacement: 커플러 세그먼트(Coupler.segments) 위치 확정 — 순차(sequential/greedy) 배치.

FP(core/floorplan.py)가 큐빗 좌표와 포트 배정까지 확정했으므로, GP는 그 위에 커플러 배선의
물리적 footprint(Coupler.segments)만 채운다. 큐빗은 절대 움직이지 않는다 — FP/GP 책임 분리는
floorplan.py 모듈 docstring 참고.

배치 단위는 한 변 segment_size_um(200um)인 정사각형이고, 전역 격자(원점 (0,0) 고정, 셀 크기
= segment_size_um)에 스냅한다. 셀 크기를 세그먼트 크기와 정확히 맞추면 "서로 다른 셀 = 절대
안 겹침"이 기하 계산 없이 성립한다(정확히 타일링되므로 셀 경계에서만 맞닿고 겹치는 넓이는
0) — 그래서 겹침 판정이 AABB 교차 계산 대신 해시셋 멤버십 확인 하나로 끝난다. eagle 기준
세그먼트 ~4,000개 규모에서 이게 그 자체로 공간 인덱스 역할을 한다(격자 버킷과 동일한 구조) —
따로 쿼드트리 등을 얹을 필요가 없다.

세그먼트는 커플러 하나당 체인(idx 0..k-1)으로, FP가 확정해 둔 배치 후보 영역
(ChipState.coupler_regions, core/floorplan.py의 assign_ports 직후 Coupler.region()으로
계산됨)을 지그재그(boustrophedon)로 채워 넣는다 — CPW meander가 실제로 짧은 변 방향으로
접히며 넓게 퍼지는 물리적 모양(region()의 자체 주석과 동일 근거)을 격자 위에서 흉내낸 것이다.

2단계 구조: FP가 커플러 경계 박스를 "선언"하고(region()이 여기선 더 이상 힌트가 아니라
GP의 하드 탐색 범위다), GP는 그 박스 밖으로 한 걸음도 안 나간다 — 박스 안에 k개 세그먼트가
안 들어가면 그 커플러는 그냥 실패 처리(segments=[]로 남기고 계속)한다. 이전 버전은 박스를
못 채우면 바깥으로 넓혀 재시도했는데, 그러면 커플러들이 서로의 영역을 침범하며 die 전체로
흩어질 수 있었다(체인 인접 쌍 8~13%가 먼 "점프"였던 것도 그 결과) — 박스를 하드 제약으로
못박으면 애초에 그런 흩어짐 자체가 구조적으로 불가능해진다.

장애물 판정(Coupler.is_own_qubit 기준, 2026-09-17 추가): 커플러 박스는 자기 큐빗의 포트에서
시작하므로 자기 자신의 q1/q2 풋프린트와 겹치는 건 정상이다(패드 연결) — func/compute.py의
compute_drc도 qc DRC에서 이 경우를 제외한다. GP도 같은 정의를 써서 자기 큐빗 셀은 장애물에서
뺀다 — 이전엔 GP만 이 규칙을 안 따라서(모든 큐빗을 무차별로 차단) 박스 용량을 체계적으로
과소평가했다: 실패했던 커플러의 점유 셀 중 95.8%가 실은 이 own-qubit 겹침이었고, 이 예외를
넣자 실패율이 27.3%에서 크게 줄었다(정확한 수치는 이 모듈을 부른 검증 스크립트/문서 참고).
"""
import logging
import math
import random
from dataclasses import replace

from core.state import ChipState, Coupler, Qubit, Segment


# 세그먼트를 둘 자리를 끝내 못 찾은 커플러가 하나둘 있는 건 정상적인 부분 실패로 다룬다
# (segments=[]로 남기고 계속 — bbox()가 None이 되어 DRC의 num_unplaced_couplers로 잡힘).
# 이 예외는 "칩 전체"가 실패로 간주될 때만(커플러 전부 실패) 올라오고, run()이 그걸 잡아
# skipped에 기록한 뒤 다음 칩으로 넘어간다 — Floorplan.PlacementInfeasibleError와 같은
# 패턴이되, 실패 사유의 성격이 완전히 달라(교차/간격 대 세그먼트 배치) 이 모듈에 따로 둔다.
class SegmentPlacementInfeasibleError(Exception):
    pass


class GlobalPlacement:
    def __init__(self, params):
        self.params = params
        self.skipped: list[tuple[str, str]] = []

    # 칩 목록을 배치한다. SegmentPlacementInfeasibleError가 난 칩(커플러 전부 실패)은
    # skipped에 (이름, 사유)로 기록하고 건너뛴다 — 한 칩 실패로 전체가 죽지 않는다.
    def run(self, states: list[ChipState]) -> list[ChipState]:
        out: list[ChipState] = []
        for state in states:
            try:
                out.append(self._place(state))
            except SegmentPlacementInfeasibleError as e:
                self.skipped.append((state.processor_name, str(e)))
        return out

    def _place(self, state: ChipState) -> ChipState:
        cell = float(self.params.segment_size_um)
        qubit_owner = _qubit_owner_cells(state.qubits, cell, state.chip_width, state.chip_height)
        segment_occupied: set[tuple[int, int]] = set()

        order_strategy = getattr(self.params, "gp_coupler_order", "shortest_first")
        order = _coupler_order(state.couplers, state.qubits, state.port_assignment, order_strategy)

        # segment_occupied는 칩 하나 안에서 커플러들이 순서대로 공유하며 누적하는 상태다 —
        # 이전 커플러가 쓴 셀은 다음 커플러의 탐색에서 항상 점유된 것으로 보여야 순차 배치가
        # "이미 놓인 걸 피해서 놓는다"는 원칙대로 동작한다(재배치/백트래킹 없음). qubit_owner는
        # 반대로 커플러마다 판정이 달라진다(자기 큐빗은 장애물이 아님, _cell_blocked 참고)라
        # 처음부터 아예 별도 자료구조로 뒀다 — 하나의 occupied 셋에 섞으면 "누가 이 셀을
        # 막았는지"를 커플러별로 다시 해석할 방법이 없어진다.
        new_couplers: dict[tuple[int, int], Coupler] = {}
        failed = 0
        for key in order:
            coupler = state.couplers[key]
            segs = _place_chain(coupler, state.qubits, state.port_assignment.get(key),
                                 state.coupler_regions.get(key), cell,
                                 state.chip_width, state.chip_height,
                                 qubit_owner, segment_occupied)
            if segs is None:
                failed += 1
                segs = []
            new_couplers[key] = replace(coupler, segments=segs)

        n_total = len(state.couplers)
        logging.info(
            "[GP] %s: couplers=%d placed=%d failed=%d segments=%d",
            state.processor_name, n_total, n_total - failed, failed,
            sum(len(c.segments) for c in new_couplers.values()),
        )

        if n_total > 0 and failed == n_total:
            raise SegmentPlacementInfeasibleError(
                f"{state.processor_name}: 커플러 {n_total}개 전부 세그먼트를 배치하지 못했습니다."
            )

        return replace(state, couplers=new_couplers)


# ---------------------------------------------------------------------------
# 큐빗 차단 셀
# ---------------------------------------------------------------------------

def _die_max_index(dim: float, cell: float) -> int:
    # dim 안에 완전히 들어가는 셀의 개수는 floor(dim/cell)개(인덱스 0..count-1) — 마지막
    # 유효 인덱스는 그보다 하나 작다(예: dim=9300, cell=200이면 셀은 46개(0..45)만 다
    # 들어가고, i=46은 [9200,9400]이라 9300을 넘어간다 — count 자체를 valid index로 쓰면
    # 이 마지막 한 칸이 die 밖으로 삐져나간다). +1e-9는 dim이 cell의 정확한 배수인데
    # 부동소수 오차로 몫이 살짝 못 미치는 경우(예: 9200/200이 45.99999...로 나오는 경우)를
    # 보정해 count를 한 칸 깎아 먹지 않기 위함이다.
    count = int(math.floor(dim / cell + 1e-9))
    return count - 1


# 큐빗 풋프린트와 조금이라도 겹치는 셀마다 "어느 큐빗이 겹치는지"를 기록한다(전부 차단하는
# 납작한 set이 아니라 cell -> {qubit_id,...} — 여러 큐빗이 한 셀에 동시에 걸치는 경우는
# min_qubit_spacing_um=400 하에서는 현실적으로 안 생기지만, 겹칠 수 있다고 가정해도 안전
# 하도록 set으로 둔다). 어떤 큐빗이 "장애물로 치는 큐빗"인지는 커플러마다 다르므로
# (Coupler.is_own_qubit — 자기 큐빗은 장애물이 아님) 여기서 미리 걸러내지 않는다 —
# 그 판정은 _cell_blocked()에서 커플러를 알고 있을 때만 한다.
def _qubit_owner_cells(
    qubits: dict[int, Qubit], cell: float, chip_w: float, chip_h: float,
) -> dict[tuple[int, int], set[int]]:
    i_max = _die_max_index(chip_w, cell)
    j_max = _die_max_index(chip_h, cell)
    owner: dict[tuple[int, int], set[int]] = {}
    for qid, q in qubits.items():
        x0, x1 = q.x - q.w / 2.0, q.x + q.w / 2.0
        y0, y1 = q.y - q.h / 2.0, q.y + q.h / 2.0
        i_lo = max(0, int(math.floor(x0 / cell + 1e-9)))
        i_hi = min(i_max, int(math.floor((x1 - 1e-9) / cell)))
        j_lo = max(0, int(math.floor(y0 / cell + 1e-9)))
        j_hi = min(j_max, int(math.floor((y1 - 1e-9) / cell)))
        for i in range(i_lo, i_hi + 1):
            for j in range(j_lo, j_hi + 1):
                owner.setdefault((i, j), set()).add(qid)
    return owner


# cell이 이 커플러에게 장애물인지. 다른 커플러의 세그먼트는 무조건 장애물이다(segment_occupied
# 는 이미 그 시점까지 놓인 것만 들어있음). 큐빗은 Coupler.is_own_qubit()로 자기 큐빗을 뺀
# 나머지가 하나라도 걸리면 장애물이다 — func/compute.py의 compute_drc가 qc DRC에서 쓰는
# 것과 동일한 정의(Coupler.is_own_qubit 참고)라 여기서 "자기 큐빗도 장애물"로 잘못 세면
# 그 정의와 어긋난 채로 박스 용량을 실제보다 작게 계산하게 된다(2026-09-17 이전 버전의
# 버그 — 실패 커플러의 점유 셀 95.8%가 실은 own-qubit 겹침이었다).
def _cell_blocked(
    coupler: Coupler, cell_xy: tuple[int, int],
    qubit_owner: dict[tuple[int, int], set[int]], segment_occupied: set[tuple[int, int]],
) -> bool:
    if cell_xy in segment_occupied:
        return True
    blockers = qubit_owner.get(cell_xy)
    if not blockers:
        return False
    return any(not coupler.is_own_qubit(qid) for qid in blockers)


# ---------------------------------------------------------------------------
# 커플러 처리 순서
# ---------------------------------------------------------------------------
#
# 후보 셋: (a) 세그먼트 수 많은 순, (b) 포트 간 거리 짧은 순, (c) 여유 공간 적은 순.
# 포트 거리를 골랐다 — 01_mainref 라우터에서 "짧은 넷 먼저"가 실측으로 확실히 이겼다는
# 선례(2/624 -> 0/624 실패)와 원리가 같다: 긴 커플러는 region()이 넓고 우회할 후보 셀도
# 많아 나중에 놓여도 자리를 찾기 쉽지만, 짧은 커플러는 애초에 두 포트 사이 좁은 영역에
# 갇혀 있어 그 자리를 다른 커플러가 먼저 잠식하면 되돌릴 방법이 없다(이 구현은 백트래킹이
# 없는 순차 배치라 더더욱). (c)(여유 공간)는 실제로 배치를 해봐야 알 수 있는 값이라 "순서를
# 정하기 위한" 기준으로 쓸 수 없어 제외했다 — (a)(세그먼트 수)는 계산은 쉽지만 이 저장소
# 6칩에서 세그먼트 수와 포트 거리가 강하게 상관돼 있어(둘 다 결국 coupler.l에서 나옴)
# 별도 기준으로서의 정보량이 적다고 보고 포트 거리로 통일했다.
#
# shortest_first/longest_first/random 세 값 다 config/params.json의 gp_coupler_order
# 설명에 실측 비교 결과와 최종 선택 근거를 적었다(6칩 x 3전략 비교, 커밋 메시지 참고).
def _port_distance(
    couplers: dict[tuple[int, int], Coupler], qubits: dict[int, Qubit],
    port_assignment: dict[tuple[int, int], tuple[str, str]], key: tuple[int, int],
) -> float:
    c = couplers[key]
    q1, q2 = qubits[c.q1], qubits[c.q2]
    p1_name, p2_name = port_assignment[key]
    p1, p2 = q1.ports[p1_name], q2.ports[p2_name]
    return math.hypot(p2[0] - p1[0], p2[1] - p1[1])


def _coupler_order(
    couplers: dict[tuple[int, int], Coupler], qubits: dict[int, Qubit],
    port_assignment: dict[tuple[int, int], tuple[str, str]], strategy: str,
) -> list[tuple[int, int]]:
    keys = sorted(couplers)
    if strategy == "shortest_first":
        return sorted(keys, key=lambda k: _port_distance(couplers, qubits, port_assignment, k))
    if strategy == "longest_first":
        return sorted(keys, key=lambda k: _port_distance(couplers, qubits, port_assignment, k), reverse=True)
    if strategy == "random":
        rng = random.Random(0)  # 고정 시드 — 매 실행 결과가 바뀌면 재현/디버깅이 불가능해진다
        shuffled = list(keys)
        rng.shuffle(shuffled)
        return shuffled
    raise ValueError(f"unknown gp_coupler_order: {strategy!r}")


# ---------------------------------------------------------------------------
# 세그먼트 체인 배치
# ---------------------------------------------------------------------------

# 박스가 이제 하드 제약(넘어가면 안 됨)이라, 셀이 [lo,hi] 안에 "완전히" 들어갈 때만
# 유효하다(i*cell >= lo and (i+1)*cell <= hi) — 걸치기만 하는 셀은 세그먼트가 박스 밖으로
# 삐져나온다. Floorplan이 coupler_regions를 이 격자에 이미 스냅해 두므로(_snap_box_to_grid,
# core/floorplan.py) 실전에서는 lo/hi가 정확히 cell의 배수라 손실이 없다 — 여기서도 굳이
# "완전 포함"으로 계산하는 건 스냅이 안 된 입력(예: 테스트에서 직접 넘긴 임의 박스)에도
# 안전하게 동작하도록 방어적으로 짠 것.
def _region_index_range(lo: float, hi: float, cell: float, i_max: int) -> tuple[int, int] | None:
    i_lo = max(0, int(math.ceil(lo / cell - 1e-9)))
    i_hi = min(i_max, int(math.floor(hi / cell + 1e-9)) - 1)
    if i_hi < i_lo:
        return None
    return i_lo, i_hi


# region 박스 안 셀을 지그재그(boustrophedon)로 훑는 순서를 만든다. 긴 변을 진행축(primary),
# 짧은 변을 지그재그축(secondary)으로 삼는다 — Coupler.region()의 "미앤더는 포트-포트
# 방향이 아니라 짧은 변 방향으로 퍼진다"는 근거를 격자 위에서 그대로 따른 것이다. p1과
# 가까운 모서리에서 시작해 진행축을 따라가며, 한 줄(secondary)을 다 훑을 때마다 방향을
# 뒤집는다 — 그러면 줄이 바뀌는 지점에서도 이전 셀과 여전히 인접(1칸 차)이 유지된다.
def _boustrophedon_cells(
    x0: float, x1: float, y0: float, y1: float, p1: tuple[float, float],
    cell: float, i_max: int, j_max: int,
) -> list[tuple[int, int]]:
    ri = _region_index_range(x0, x1, cell, i_max)
    rj = _region_index_range(y0, y1, cell, j_max)
    if ri is None or rj is None:
        return []
    i0, i1 = ri
    j0, j1 = rj

    if (i1 - i0) >= (j1 - j0):
        p_lo, p_hi, s_lo, s_hi = i0, i1, j0, j1
        p_near = p1[0] <= (i0 + i1 + 1) * cell / 2.0
        s_near = p1[1] <= (j0 + j1 + 1) * cell / 2.0
        make = lambda p, s: (p, s)
    else:
        p_lo, p_hi, s_lo, s_hi = j0, j1, i0, i1
        p_near = p1[1] <= (j0 + j1 + 1) * cell / 2.0
        s_near = p1[0] <= (i0 + i1 + 1) * cell / 2.0
        make = lambda p, s: (s, p)

    primaries = range(p_lo, p_hi + 1) if p_near else range(p_hi, p_lo - 1, -1)
    cells: list[tuple[int, int]] = []
    forward = s_near
    for p in primaries:
        seconds = range(s_lo, s_hi + 1) if forward else range(s_hi, s_lo - 1, -1)
        for s in seconds:
            cells.append(make(p, s))
        forward = not forward
    return cells


# 커플러 하나의 세그먼트 체인을 배치한다. 박스(FP가 확정한 coupler_regions[key]) 안에서만
# 찾는다 — 밖으로 넓히는 재시도는 없다(2단계 구조: 박스를 넓히는 건 FP의 권한이지 GP의
# 권한이 아니다). 장애물 판정은 _cell_blocked()(자기 큐빗 제외, 제3 큐빗·다른 커플러
# 세그먼트는 포함)로 한다. 성공하면 Segment 리스트를 돌려주고 segment_occupied를 그
# 셀들만큼 갱신한다(호출부가 따로 갱신할 필요 없음 — qubit_owner는 애초에 커플러가 안
# 바꾸므로 갱신 대상이 아니다). 실패하면 아무것도 건드리지 않고 None을 돌려준다 — 이때
# "실패"는 이 박스 하나에 국한된 사실이라, 같은 박스를 먼저 차지한 다른 커플러의 세그먼트를
# 피해 자기 박스의 남은 자리를 쓰는 것까지는 이 함수의 free-cell 필터링만으로 이미 된다 —
# 백트래킹/재배치는 하지 않는다(순차 배치 원칙).
def _place_chain(
    coupler: Coupler, qubits: dict[int, Qubit], assignment: tuple[str, str] | None,
    box: tuple[float, float, float, float] | None,
    cell: float, chip_w: float, chip_h: float,
    qubit_owner: dict[tuple[int, int], set[int]], segment_occupied: set[tuple[int, int]],
) -> list[Segment] | None:
    k = coupler.num_segments
    if k == 0:
        return []
    if assignment is None or box is None:
        return None  # FP가 포트/박스를 못 정한 커플러 — 이론상 skipped 칩에서만 나오므로 여기 안 옴

    q1 = qubits[coupler.q1]
    p1_name, _p2_name = assignment
    p1 = q1.ports[p1_name]

    x0, x1, y0, y1 = box
    i_max = _die_max_index(chip_w, cell)
    j_max = _die_max_index(chip_h, cell)

    cells = _boustrophedon_cells(x0, x1, y0, y1, p1, cell, i_max, j_max)
    free = [c for c in cells if not _cell_blocked(coupler, c, qubit_owner, segment_occupied)]
    if len(free) < k:
        return None

    chosen = free[:k]
    segment_occupied.update(chosen)
    return [Segment(idx=idx, x=(i + 0.5) * cell, y=(j + 0.5) * cell)
            for idx, (i, j) in enumerate(chosen)]
