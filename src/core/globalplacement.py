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
계산됨) 안에서 **q1의 포트 셀 -> q2의 포트 셀로 이어지는, 4-인접한 단순 경로**를 찾아
채운다(_find_chain_path). RT(core/router.py v3)가 이제 Segment.idx 순서를 그대로 이어
꺾은선을 만들기 때문에(셀 안에서 접지 않는다 — QPlacer Figure 8(e), router.py 모듈
docstring 참고) idx가 연속인 두 세그먼트가 격자상 이웃이 아니면 RT가 그 커플러를 통째로
실패 처리한다.

2026-09-20 이전엔 박스를 지그재그(boustrophedon)로 그냥 "훑으며 빈 셀 채우기"였다 — 진행
방향을 따라 셀을 순서대로 보되, 막힌 셀(own-qubit 몸체/제3 큐빗/다른 커플러 세그먼트)은
그냥 건너뛰고 다음 빈 셀을 집었다. 이게 "제일 먼저 만나는 k개의 빈 셀"은 보장해도 "그
k개가 서로 이어져 있다"는 전혀 보장하지 않는다 — 막힌 셀 하나를 건너뛰는 순간 idx가
연속인 두 세그먼트 사이에 격자상 빈틈이 생긴다. 실측(RT v3 도입 후): 커플러 기준
22.5~82.8%가 이런 "점프"를 하나 이상 가졌고, RT 라우팅 성공률이 40.0%까지 떨어졌다 —
GP가 만드는 "체인"이 사실 체인이 아니었다는 뜻이다. 그래서 배치 자체를 "빈 셀 채우기"가
아니라 "포트에서 포트로 이어지는 경로 탐색"으로 바꿨다 — 막힌 셀은 건너뛰는 게 아니라
경로가 아예 피해 가야 한다(_find_chain_path 참고).

2단계 구조: FP가 커플러 경계 박스를 "선언"하고(region()이 여기선 더 이상 힌트가 아니라
GP의 하드 탐색 범위다), GP는 그 박스 밖으로 한 걸음도 안 나간다 — 박스 안에 k개 세그먼트가
안 들어가면 그 커플러는 그냥 실패 처리(segments=[]로 남기고 계속)한다. 이전 버전은 박스를
못 채우면 바깥으로 넓혀 재시도했는데, 그러면 커플러들이 서로의 영역을 침범하며 die 전체로
흩어질 수 있었다(체인 인접 쌍 8~13%가 먼 "점프"였던 것도 그 결과) — 박스를 하드 제약으로
못박으면 애초에 그런 흩어짐 자체가 구조적으로 불가능해진다.

장애물 판정(core.state.coupler_own_port_cell 기준): 커플러 박스는 자기 큐빗의 포트에서
시작하므로 그 포트 근처와 겹치는 건 정상이다(패드 연결) — func/compute.py의 compute_drc도
qc DRC에서 이 경우를 제외한다(둘 다 같은 함수 하나만 쓴다). GP는 그 정의로 자기 큐빗 셀을
장애물에서 뺀다.

이 예외의 범위가 두 번 바뀌었다. (1) 2026-09-17: 그 전엔 GP가 모든 큐빗을 무차별로
차단해 박스 용량을 체계적으로 과소평가했다 — 실패했던 커플러의 점유 셀 중 95.8%가 실은
own-qubit 겹침이었고, "자기 큐빗은 전부 허용"으로 고치자 실패율이 27.3%에서 크게 줄었다.
(2) 2026-09-20: "자기 큐빗은 전부 허용"이 너무 관대했다는 게 드러났다 — 실측 결과
세그먼트가 포트 근처가 아니라 큐빗 중심 코앞(반폭 200um의 최대 91%인 182um)까지 파고들었고,
칩당 1~28개는 몸체 안에 완전히 들어가 있었다(docs/20260920_segment_model_review.md 발견
3). 큐빗 몸체는 실제로는 금속이라 배선이 지날 수 없다 — 그래서 예외를 "배정 포트 근처
셀만"으로 다시 좁혔다(coupler_own_port_cell). 이 좁힘은 (1)이 넓혀준 용량의 일부를 다시
가져가므로 실패율이 다시 오를 것으로 예상된다 — 정확한 수치는 이 변경의 커밋/보고
메시지에 남긴다.
"""
import logging
import math
import random
from dataclasses import replace

from core.floorplan import segments_cross
from core.state import ChipState, Coupler, Qubit, Segment, coupler_own_port_cell, _port_edge_axis


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
        # 칩 이름 -> leg_ok가 실제로 후보를 걸러낸 횟수(2026-09-21 검증 보고용) — 이웃
        # 방향 포트로 바꾼 뒤 케이스 A가 거의 사라져 leg_ok가 얼마나 덜 걸리는지 확인.
        self.leg_ok_blocks: dict[str, int] = {}

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
        qubit_rects = {
            qid: (q.x - q.w / 2.0, q.x + q.w / 2.0, q.y - q.h / 2.0, q.y + q.h / 2.0)
            for qid, q in state.qubits.items()
        }
        segment_occupied: set[tuple[int, int]] = set()

        order_strategy = getattr(self.params, "gp_coupler_order", "shortest_first")
        order = _coupler_order(state.couplers, state.ports, order_strategy)

        # segment_occupied는 칩 하나 안에서 커플러들이 순서대로 공유하며 누적하는 상태다 —
        # 이전 커플러가 쓴 셀은 다음 커플러의 탐색에서 항상 점유된 것으로 보여야 순차 배치가
        # "이미 놓인 걸 피해서 놓는다"는 원칙대로 동작한다(재배치/백트래킹 없음). qubit_owner는
        # 반대로 커플러마다 판정이 달라진다(자기 큐빗은 장애물이 아님, _cell_blocked 참고)라
        # 처음부터 아예 별도 자료구조로 뒀다 — 하나의 occupied 셋에 섞으면 "누가 이 셀을
        # 막았는지"를 커플러별로 다시 해석할 방법이 없어진다.
        new_couplers: dict[tuple[int, int], Coupler] = {}
        failed = 0
        short = 0  # 경로는 찾았지만 k개에 못 미친 커플러(요청 2절: 실패 아니라 미달로 다룸)
        leg_ok_blocks = 0  # leg_ok가 실제로 후보를 걸러낸 횟수(검증 보고용) — _place_chain에 누산기로 전달
        for key in order:
            coupler = state.couplers[key]
            segs, blocks = _place_chain(coupler, state.qubits, state.ports.get(key),
                                         state.coupler_regions.get(key), cell,
                                         state.chip_width, state.chip_height,
                                         qubit_owner, segment_occupied, qubit_rects)
            leg_ok_blocks += blocks
            if segs is None:
                failed += 1
                segs = []
            elif 0 < len(segs) < coupler.num_segments:
                short += 1
            new_couplers[key] = replace(coupler, segments=segs)

        n_total = len(state.couplers)
        logging.info(
            "[GP] %s: couplers=%d placed=%d failed=%d 길이미달=%d segments=%d leg_ok차단=%d",
            state.processor_name, n_total, n_total - failed, failed, short,
            sum(len(c.segments) for c in new_couplers.values()), leg_ok_blocks,
        )
        self.leg_ok_blocks[state.processor_name] = leg_ok_blocks

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
# 는 이미 그 시점까지 놓인 것만 들어있음). 큐빗은 core.state.coupler_own_port_cell()로
# "자기 큐빗이면서 배정 포트 근처"인 경우만 제외한다 — func/compute.py의 compute_drc가
# qc DRC에서 쓰는 것과 정확히 같은 함수라, 여기서 "놓을 수 있다"고 판단한 배치를 DRC가
# 다시 "위반"으로 잡아내는 모순이 생기지 않는다. 2026-09-20 이전엔 자기 큐빗이면 몸체
# 전체를 허용했는데, 실측 결과 세그먼트가 큐빗 중심 코앞까지 파고들었다
# (docs/20260920_segment_model_review.md 발견 3) — coupler_own_port_cell()로 좁힌 이유는
# 그 함수 자체의 docstring 참고. 이 좁힘으로 GP의 박스 용량이 다시 줄어든다(예전 own-qubit
# 확장이 실패율을 27.3%→2.6%로 낮췄던 것의 부분적 반대 방향) — 그만큼 실패율이 다시
# 오를 것으로 예상되고, 실측치는 이 변경의 커밋/보고 메시지에 남긴다.
# 선분(a,b)가 어떤 큐빗의 AABB 내부를 실제로 지나는지 — cell_xy vs qubit AABB(사각형 대
# 사각형) 겹침과는 다른 판정이다. 포트는 큐빗 경계(변) 위의 점이고, 그 포트와 이어지는
# 첫/마지막 인접 셀의 중심은 그 셀 자체가 어떤 큐빗과도 안 겹치도록 이미 보장돼 있는데도
# (own-qubit 예외는 포트가 속한 셀 하나에만 적용되므로), "경계 위의 점 -> 안 겹치는 셀의
# 중심"을 잇는 직선은 여전히 큐빗의 볼록한 내부를 스쳐 지나갈 수 있다(포트가 모서리 근처에
#있고 다음 셀이 인접한 변 너머에 있을 때 — 실측: core/router.py가 첫/마지막 연결점을
# 포트 좌표로 바꾼 뒤에도 6칩에서 347개 라우팅된 커플러 중 148개, leg 171건이 여전히
# 큐빗 AABB를 지났다). 두 셀 중심 사이(중간 세그먼트끼리)는 이 문제가 없다 — 둘 다 이미
# 어떤 큐빗과도 안 겹치는 두 인접 정사각형의 합집합 안에 직선이 갇히므로(각 셀이 큐빗과
# 안 겹치면 합집합도 안 겹친다) 안전하다. 그래서 이 검사는 포트<->첫/마지막 인접 셀 구간
# (아래 _place_chain의 leg_ok)에만 쓴다.
def _segment_crosses_any_qubit(
    a: tuple[float, float], b: tuple[float, float],
    qubit_rects: dict[int, tuple[float, float, float, float]],
) -> bool:
    for x0, x1, y0, y1 in qubit_rects.values():
        if (x0 < a[0] < x1 and y0 < a[1] < y1) or (x0 < b[0] < x1 and y0 < b[1] < y1):
            return True
        corners = ((x0, y0), (x1, y0), (x1, y1), (x0, y1))
        if any(segments_cross(a, b, corners[i], corners[(i + 1) % 4]) for i in range(4)):
            return True
    return False


def _cell_blocked(
    coupler: Coupler, cell_xy: tuple[int, int], cell: float,
    qubits: dict[int, Qubit],
    ports: tuple[tuple[float, float], tuple[float, float]] | None,
    qubit_owner: dict[tuple[int, int], set[int]], segment_occupied: set[tuple[int, int]],
) -> bool:
    if cell_xy in segment_occupied:
        return True
    blockers = qubit_owner.get(cell_xy)
    if not blockers:
        return False
    cell_aabb = (cell_xy[0] * cell, (cell_xy[0] + 1) * cell, cell_xy[1] * cell, (cell_xy[1] + 1) * cell)
    return any(
        not coupler_own_port_cell(coupler, qid, ports, cell_aabb)
        for qid in blockers
    )


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
    ports: dict[tuple[int, int], tuple[tuple[float, float], tuple[float, float]]], key: tuple[int, int],
) -> float:
    p1, p2 = ports[key]
    return math.hypot(p2[0] - p1[0], p2[1] - p1[1])


def _coupler_order(
    couplers: dict[tuple[int, int], Coupler],
    ports: dict[tuple[int, int], tuple[tuple[float, float], tuple[float, float]]], strategy: str,
) -> list[tuple[int, int]]:
    keys = sorted(couplers)
    if strategy == "shortest_first":
        return sorted(keys, key=lambda k: _port_distance(ports, k))
    if strategy == "longest_first":
        return sorted(keys, key=lambda k: _port_distance(ports, k), reverse=True)
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


# 격자 위에서 start(q1 포트 셀)에서 end(q2 포트 셀)까지, 정확히(또는 그보다 짧게) k개의
# 셀을 지나는 4-인접 단순 경로(같은 셀을 두 번 지나지 않음 — 세그먼트는 물리적 정사각형
# 이라 겹칠 수 없다)를 찾는다. RT(core/router.py v3)가 Segment.idx 순서를 그대로 이어
# 꺾은선을 만들므로, 이 함수가 돌려주는 경로의 순서가 곧 idx 순서가 되고, 인접한 두
# 원소는 항상 격자상 이웃이라는 게 RT가 기대하는 유일한 불변식이다.
#
# 왜 DFS+백트래킹인가: "정확히 k개짜리 단순 경로 찾기"는 일반적으로 NP-hard다(해밀턴
# 경로 문제의 변형 — 그래프에 장애물이 있고 길이가 그래프 크기보다 작아도 마찬가지).
# 이 저장소 규모(박스당 자유 셀 수십 개, 칩당 커플러 최대 144개)에서 정확 해를 매번
# 구하는 건 현실적이지 않다. 후보 셋(요청 1절)을 이렇게 평가했다:
#   - 길이를 상태에 넣은 BFS: "방문한 셀 집합"까지 상태에 넣어야 단순 경로를 보장하는데,
#     그러면 상태 공간이 사실상 전수 탐색과 같아진다(방문 집합 자체가 지수개) — BFS의
#     "레벨별로 짧은 것부터"라는 장점이 여기선 안 산다(짧은 경로가 아니라 "정확히 k개"가
#     필요하므로).
#   - 탐욕 + 국소 수정: 매 걸음 목표 방향으로만 가면 막다른 길에 몰렸을 때 되돌아갈
#     방법이 없다 — "국소 수정"이 사실상 백트래킹인데, 그럴 거면 처음부터 백트래킹을
#     제대로 설계하는 게 낫다.
#   - DFS+백트래킹(채택): 아래 두 가지치기를 더하면 실전 규모에서 빠르다(실측은
#     GlobalPlacement._place 호출부의 커밋/보고 메시지 참고). 이 저장소가 이미 비슷한
#     원칙(예: LG의 "후보 배제, 페널티 아님")을 여러 곳에 쓰고 있어 일관적이기도 하다.
#
# 가지치기 둘 다 "참인 해를 걸러내지 않는다"(안전한 필요조건)는 게 핵심이다:
#   1) 맨해튼 거리 가지치기: 어떤 셀에서 end까지 맨해튼 거리보다 적은 걸음으로는 절대
#      못 간다 — 장애물은 필요 걸음 수를 늘릴 뿐 줄이지 않으므로 이 하한은 장애물 유무와
#      무관하게 항상 성립한다.
#   2) 홀짝 가지치기: 격자는 이분 그래프라 한 걸음마다 맨해튼 거리의 홀짝이 뒤집힌다 —
#      "남은 예산 - 남은 거리"가 홀수면 그 지점에서 정확히 그 예산 안에 end로 못 들어온다.
#
# 여유(slack = 남은 예산 - 맨해튼 거리)가 있으면 end에서 "먼" 이웃부터 시도한다(경로를
# 길게 뽑아 k에 다가가려는 목적) — 여유가 0이 되는 순간부터는 매 걸음이 거리를 정확히
# 1씩 줄여야만 하므로(그러지 않으면 예산 초과) end로 "직행"하는 이웃만 시도한다. 이
# "여유 있을 때 방황, 없으면 직행" 전략은 백트래킹을 여유가 있는 초반 구간에 국한시켜,
# 박스가 심하게 막혀 있지 않은 한 대부분 빠르게 끝난다.
#
# max_expansions로 DFS 노드 확장 수를 제한한다 — 그 예산 안에 정확히 k개를 못 찾으면
# 그때까지 찾은 "end에 닿은 가장 긴 경로"를 대신 쓴다(요청 2절: k개보다 짧은 경로만
# 가능하면 그대로 배치하고 미달을 보고한다 — 시간 예산 소진도 같은 원칙으로 다룬다,
# 억지로 더 찾지 않는다). end에 닿는 경로 자체를 하나도 못 찾으면(장애물에 막혀 연결
# 자체가 불가능하거나, k가 start-end 최短 거리보다 작아 애초에 도달 불가능하면) None —
# 이 경우는 "박스 밖으로 나가서라도 k개를 채운다"가 아니라 그대로 실패로 본다(요청 2절:
# "k개를 넘기려고 박스 밖으로 나가지 마라" — 박스 안에서 아예 길이 없으면 실패다).
def _find_chain_path(
    start: tuple[int, int], end: tuple[int, int], is_free, k: int,
    max_expansions: int = 20000,
    leg_ok=None,
) -> list[tuple[int, int]] | None:
    if k <= 0:
        return None
    if not is_free(start) or not is_free(end):
        return None
    if start == end:
        return [start]

    def manhattan(a: tuple[int, int], b: tuple[int, int]) -> int:
        return abs(a[0] - b[0]) + abs(a[1] - b[1])

    # 격자는 이분 그래프라 start->end 단순 경로의 길이(셀 개수)는 항상 같은 홀짝성만
    # 가능하다 — 걸음 수(길이-1)가 manhattan(start,end)와 같은 홀짝이어야 하므로(각 걸음이
    # 거리를 ±1만 바꾸니까), 길이 자체는 항상 manhattan(start,end)+1과 같은 홀짝이다. k가
    # 그 홀짝과 다르면 "정확히 k개"는 장애물과 무관하게 수학적으로 아예 불가능하다 — 실측
    # (grid_25 재현: k=10, manhattan=4일 때 9-4=5로 홀수)으로 처음엔 이걸 놓쳐서, 매 반복의
    # 가지치기가 "고정된 k" 기준으로 계산되는 바람에 틀린 홀짝이 첫 걸음부터 모든 분기를
    # 막아버려 "더 짧은 경로조차" 못 찾고 완전 실패로 돌아갔다(요청 2절이 원하는 "짧은
    # 경로라도 배치"를 아예 시도조차 못 함). 고친 방법: 탐색을 시작하기 전에 목표 길이를
    # k 또는 k-1 중 실제로 가능한 홀짝으로 맞춘다 — k와 k-1은 항상 서로 다른 홀짝이므로
    # 반드시 둘 중 하나는 맞는다. 이렇게 하면 가지치기 공식이 항상 "달성 가능한 목표"
    # 기준으로 서기 때문에(아래 dfs 내부의 홀짝 불변식이 시작점부터 성립하고, 한 걸음마다
    # (remaining - 거리)가 0 또는 -2만큼만 바뀌므로 그 불변식이 끝까지 유지된다), 진짜
    # 장애물 때문에 막히는 경우와 "애초에 숫자가 안 맞아서" 막히는 경우가 더 이상 섞이지
    # 않는다.
    target = k
    if (target - 1 - manhattan(start, end)) % 2 != 0:
        target -= 1
    if target <= 0:
        return None

    best_path: list[tuple[int, int]] | None = None
    path = [start]
    visited = {start}
    expansions = [0]
    k = target

    def dfs() -> bool:  # True = 그만 찾아도 됨(정확히 k개 성공 또는 예산 소진)
        nonlocal best_path
        expansions[0] += 1
        if expansions[0] > max_expansions:
            return True

        cur = path[-1]
        if cur == end:
            if best_path is None or len(path) > len(best_path):
                best_path = list(path)
            if len(path) == k:
                return True

        if len(path) >= k:
            return False

        remaining = k - len(path)
        candidates = []
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            n = (cur[0] + dx, cur[1] + dy)
            if n in visited or not is_free(n):
                continue
            # cur가 start이거나 n이 end면(둘 다인 경우 — start-end 직행 포함) 이 걸음은
            # RT가 실제로 포트 좌표를 잇는 leg가 된다 — leg_ok가 그 leg를 검사한다(모듈
            # 상단 _segment_crosses_any_qubit 참고). 순수 중간 구간(둘 다 아님)은 이미
            # 안전함이 증명돼 있어(어느 큐빗과도 안 겹치는 두 인접 셀의 합집합 안에만
            # 있음) 검사하지 않는다.
            if leg_ok is not None and (cur == start or n == end) and not leg_ok(cur, n):
                continue
            d = manhattan(n, end)
            slack = (remaining - 1) - d
            if slack < 0 or slack % 2 != 0:
                continue  # 이 이웃으로 가면 남은 예산 안에 end 도달 불가(안전한 가지치기)
            candidates.append((n, d))

        cur_slack = remaining - manhattan(cur, end)
        candidates.sort(key=lambda t: t[1], reverse=(cur_slack > 0))

        for n, _d in candidates:
            visited.add(n)
            path.append(n)
            done = dfs()
            path.pop()
            visited.discard(n)
            if done:
                return True
        return False

    dfs()
    return best_path


# 포트가 속한 셀이 이제 항상 박스 안에 있다는 보장이 없다 — core/state.py의 Coupler.region()
# 이 "케이스 B"(포트 쪽 축이 자기 큐빗 반대편으로 뻗어나가는 경우) 박스를 자기 큐빗
# 몸체 앞에서 멈추도록 핀을 걸면서, 포트가 걸친 셀 자체(그 셀은 큐빗 몸체와 겹친다)가
# 박스에서 빠질 수 있게 됐다. 그 경우 체인은 포트가 아니라 포트 바로 바깥의 빈 셀에서
# 시작해야 하고, 그 빈 셀 <-> 포트 사이는 별도의 리드(lead)로 잇는다(RT가 실제로 그릴
# 좌표). 이 함수가 그 "포트 -> 실제 시작 셀" 후보를 찾는다: 먼저 포트가 속한 셀을
# 시도하고(박스가 안 핀됐으면(케이스 A) 여전히 유효 — 이 경우 아무것도 안 바뀐다), 안
# 되면 그 큐빗의 중심에서 포트 쪽으로(=바깥쪽으로) 한 칸씩 옮겨가며 박스 안 + 가용 +
# (포트->그 셀 중심) 직선이 어떤 큐빗도 안 지나는 셀을 찾는다. max_shift 안에 못 찾으면
# None(이 커플러는 실패 — 억지로 더 찾지 않는다, 이 저장소의 일관된 원칙).
def _resolve_endpoint(
    port_pt: tuple[float, float], owner_qubit: Qubit,
    in_box, is_free, qubit_rects: dict[int, tuple[float, float, float, float]],
    cell: float, max_shift: int = 3,
) -> tuple[tuple[int, int], bool] | tuple[None, bool]:
    base = (int(math.floor(port_pt[0] / cell)), int(math.floor(port_pt[1] / cell)))
    if in_box(base) and is_free(base):
        return base, True  # True = 포트가 속한 셀 그대로(케이스 A와 동일한 예전 동작)
    # 큐빗 몸체 반대 방향(바깥쪽)으로 한 칸씩 옮긴다 — 2026-09-21 이웃 방향 포트로 바뀌며
    # 포트가 상/하 변에도 있을 수 있게 됐으므로, core/state.py의 _port_edge_axis로 축을
    # 매번 판별한다(예전엔 포트가 항상 좌/우 변이라 x축 고정이었다).
    axis, qubit_dir = _port_edge_axis(port_pt, owner_qubit)
    sign = -1 if qubit_dir > 0 else 1
    for step in range(1, max_shift + 1):
        cand = (base[0] + sign * step, base[1]) if axis == "x" else (base[0], base[1] + sign * step)
        if not (in_box(cand) and is_free(cand)):
            continue
        cand_pt = ((cand[0] + 0.5) * cell, (cand[1] + 0.5) * cell)
        if _segment_crosses_any_qubit(port_pt, cand_pt, qubit_rects):
            continue
        return cand, False  # False = 포트와 별개인 셀 -- RT가 리드로 잇는다
    return None, False


# 커플러 하나의 세그먼트 체인을 배치한다. 박스(FP가 확정한 coupler_regions[key]) 안에서만
# 찾는다 — 밖으로 넓히는 재시도는 없다(2단계 구조: 박스를 넓히는 건 FP의 권한이지 GP의
# 권한이 아니다). 장애물 판정은 _cell_blocked()(자기 큐빗의 배정 포트 근처만 제외, 그
# 밖의 자기 큐빗 몸체·제3 큐빗·다른 커플러 세그먼트는 전부 포함)로 한다 — _find_chain_path
# 가 그 판정 하나로 "지나갈 수 있는 셀"을 정의하므로, GP가 "놓을 수 있다"고 본 경로는
# 정의상 장애물을 피해 간다. 성공(부분 성공 포함)하면 Segment 리스트를 돌려주고
# segment_occupied를 그 셀들만큼 갱신한다. start/end 자체가 막혀 있거나(또는 _resolve_endpoint가
# 못 찾았거나) 둘이 아예 연결돼 있지 않으면 None — 이때 "실패"는 이 박스 하나에 국한된
# 사실이라, 같은 박스를 먼저 차지한 다른 커플러의 세그먼트를 피해 자기 박스의 남은 자리를
# 쓰는 것까지는 이 함수의 free-cell 판정만으로 이미 된다 — 백트래킹은 경로 탐색 내부에서만
# 하고, 커플러 사이의 재배치는 하지 않는다(순차 배치 원칙 그대로 유지).
def _place_chain(
    coupler: Coupler, qubits: dict[int, Qubit],
    ports: tuple[tuple[float, float], tuple[float, float]] | None,
    box: tuple[float, float, float, float] | None,
    cell: float, chip_w: float, chip_h: float,
    qubit_owner: dict[tuple[int, int], set[int]], segment_occupied: set[tuple[int, int]],
    qubit_rects: dict[int, tuple[float, float, float, float]],
) -> tuple[list[Segment] | None, int]:
    # 반환값의 두 번째 항목(leg_ok가 실제로 후보를 걸러낸 횟수)은 검증 보고용이다 —
    # 2026-09-21 이웃 방향 포트 도입 후 케이스 A가 거의 사라져 이 값이 크게 줄 것으로
    # 예상된다(GlobalPlacement._place가 칩별로 누산해 로그에 남긴다).
    k = coupler.num_segments
    if k == 0:
        return [], 0
    if ports is None or box is None:
        return None, 0  # FP가 포트/박스를 못 정한 커플러 — 이론상 skipped 칩에서만 나오므로 여기 안 옴

    q1, q2 = qubits[coupler.q1], qubits[coupler.q2]
    p1, p2 = ports

    x0, x1, y0, y1 = box
    i_max = _die_max_index(chip_w, cell)
    j_max = _die_max_index(chip_h, cell)
    ri = _region_index_range(x0, x1, cell, i_max)
    rj = _region_index_range(y0, y1, cell, j_max)
    if ri is None or rj is None:
        return None, 0
    i_lo, i_hi = ri
    j_lo, j_hi = rj

    def in_box(c: tuple[int, int]) -> bool:
        return i_lo <= c[0] <= i_hi and j_lo <= c[1] <= j_hi

    def is_free(c: tuple[int, int]) -> bool:
        return in_box(c) and not _cell_blocked(coupler, c, cell, qubits, ports, qubit_owner, segment_occupied)

    # start/end: 포트가 속한 셀이 박스 안 + 가용이면 그대로 쓰고(케이스 A, 예전과 동일),
    # 아니면(케이스 B — Coupler.region()이 자기 큐빗 몸체 앞에서 박스를 핀으로 멈췄음)
    # 그 바깥의 가용한 셀로 옮긴다 — _resolve_endpoint 참고. start_is_port/end_is_port는
    # "그 셀이 실제로 포트가 속한 셀인가"를 아래 leg_ok와 core/router.py(_route_coupler)에
    # 전달한다: 그럴 때만(케이스 A) 그 셀의 중심 대신 포트 좌표를 대신 쓴다(그 셀 중심이
    # 큐빗 몸체 안쪽일 수 있으므로) — 케이스 B는 그 셀 자체가 이미 어떤 큐빗과도 안 겹치므로
    # (그래서 is_free를 통과했다) 중심 그대로 써도 안전하고, RT가 포트<->그 셀 사이에
    # 실제 리드 세그먼트를 그린다(router.py 참고).
    start, start_is_port = _resolve_endpoint(p1, q1, in_box, is_free, qubit_rects, cell)
    end, end_is_port = _resolve_endpoint(p2, q2, in_box, is_free, qubit_rects, cell)
    if start is None or end is None:
        return None, 0

    # RT(core/router.py)가 폴리라인의 첫/마지막 연결점으로 셀 중심 대신 실제 포트 좌표를
    # 쓴다(2026-09-20, 관통 수정) — 그런데 포트(큐빗 경계 위의 점)에서 그 다음 셀 중심까지
    # 잇는 leg는 그 셀 자체가 큐빗과 안 겹쳐도(own-qubit 예외가 없는 자리라 _cell_blocked가
    # 이미 보장) 여전히 큐빗의 볼록한 몸체 모서리를 스칠 수 있다(모듈 상단
    # _segment_crosses_any_qubit 참고). leg_ok로 그 leg 자체를 검사해 후보에서
    # 제외한다 — "포트 셀이 아니면 몸체 안은 못 지나간다"를 셀 겹침뿐 아니라 포트로
    # 이어지는 leg까지 확장한 것이다. a/b가 각각 start/end이고 *그 셀이 실제로 포트
    # 셀일 때만*(start_is_port/end_is_port) 그 끝점을 셀 중심 대신 포트 좌표로 바꿔서
    # 검사한다 — 케이스 B(포트와 다른 셀)는 이미 큐빗과 안 겹치는 셀이라 그 중심 자체를
    # 써도 되고(_resolve_endpoint가 포트->그 셀 리드는 이미 검사해 뒀다), RT가 실제로
    # 그릴 좌표와 정확히 같은 것을 검사해야 하므로.
    blocks = [0]

    def leg_ok(a: tuple[int, int], b: tuple[int, int]) -> bool:
        a_pt = p1 if (a == start and start_is_port) else ((a[0] + 0.5) * cell, (a[1] + 0.5) * cell)
        b_pt = p2 if (b == end and end_is_port) else ((b[0] + 0.5) * cell, (b[1] + 0.5) * cell)
        ok = not _segment_crosses_any_qubit(a_pt, b_pt, qubit_rects)
        if not ok:
            blocks[0] += 1
        return ok

    path = _find_chain_path(start, end, is_free, k, leg_ok=leg_ok)
    if path is None:
        return None, blocks[0]

    segment_occupied.update(path)
    return [Segment(idx=idx, x=(i + 0.5) * cell, y=(j + 0.5) * cell)
            for idx, (i, j) in enumerate(path)], blocks[0]
