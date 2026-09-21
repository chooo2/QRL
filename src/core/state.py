import math
from dataclasses import dataclass, field

# AABB(축정렬 경계상자) 겹침 판정에 쓰는 부동소수 허용 오차 (um) — func/compute.py,
# core/legalization.py, core/detailedplacement.py가 전부 이 값을 쓴다. 물리적 clearance가
# 아니다(그 결정은 그대로다: qc는 zero-margin, 01_mainref의 500um는 도입 안 함,
# docs/20260917_gp_qc_findings.md 참고) — 1e-6um(=1피코미터)는 어떤 실제 설계 기준보다도
# 8~9자리 작아서 "진짜로 안전 여유를 준다"는 의미가 전혀 없고, 순수하게 부동소수 표현
# 오차만 흡수한다. 왜 필요한가: core/detailedplacement.py가 좌표를 배치 중심에 대해
# α배로 스케일하는데, "중심을 각각 스케일한 뒤 빼서 거리(dx)를 구하는" 경로와 "반폭을
# 각각 스케일해 더해서 합(sum_hw)을 구하는" 경로는 수학적으로 같은 값이 나와야 할 때도
# 서로 다른 반올림을 거쳐(부동소수 좌표가 최대 수만 um 규모라 float64 정밀도 한계상)
# ~1e-13um 수준으로 어긋난다 — eps 없이 엄격한 부등호(<)만 쓰면 "정확히 맞닿음"(겹침
# 아님)이던 관계가 무작위로 겹침/비겹침을 오간다(실측 사례: core/detailedplacement.py
# 압축 중 원래 정확히 맞닿아 있던 세그먼트 쌍이 이 잡음 때문에 스퓨리어스하게 겹침으로
# 잡혔다). eps를 빼는 방향(엄격하게)으로 걸어 두면 이 잡음은 흡수하면서 진짜 겹침(1e-6um
# 보다 깊은)은 그대로 잡는다.
AABB_EPS_UM = 1e-6

# Coupler.region()의 가용 셀 반복 확장이 종료를 보장받기 위한 최대 반복 횟수. 한 반복은
# 단변을 정확히 cell(segment_size_um) 한 칸만큼 넓힌다 — 그 새 줄이 제3자 큐빗에 다시
# 막히면 다음 반복에서 또 넓혀야 한다. 고정 폭(cell)만큼만 늘리므로 2026-09-20 버전
# (부족분 비례 증가, 감쇠비가 1에 가까우면 느려짐)과 달리 매 반복 진전량이 최소 0인
# 경우(그 줄이 전부 막혔을 때)만 빼면 항상 최소 1행/열은 새로 확보된다 — 큐빗 폭이
# 유한하므로(400um=cell 2개) 몇 번이면 그 큐빗을 지나쳐 막힘이 풀린다. 40은 그런 경우도
# 여유 있게 흡수하는 값(반복당 비용은 O(박스 안 셀 수 × 큐빗 수)로 작다).
_REGION_GROW_MAX_ITERS = 40


# region()의 반복 확장, 그리고 아래 _count_free_cells_in_box에서 쓰는 사각형 교집합
# 헬퍼(겹치는지 여부만 필요 — 겹치는 사각형 자체는 안 쓰므로 None 여부만 본다).
def _rect_intersect(
    a: tuple[float, float, float, float], b: tuple[float, float, float, float],
) -> tuple[float, float, float, float] | None:
    x0, x1 = max(a[0], b[0]), min(a[1], b[1])
    y0, y1 = max(a[2], b[2]), min(a[3], b[3])
    if x1 <= x0 or y1 <= y0:
        return None
    return (x0, x1, y0, y1)


# Coupler.region()의 축별 "핀"(pin) 헬퍼 — 2026-09-21 포트를 이웃 방향 기반
# (qubit_port_toward)으로 바꾸기 전엔 포트가 항상 큐빗의 좌/우 변(x축)에 있어서 이
# 로직이 x축 전용이었다. 이제 포트가 상/하 변(y축)에도 있을 수 있으므로, 두 축 어느
# 쪽이든 같은 함수로 다룬다(호출부가 x/y 각각 한 번씩 호출).
def _pin_snap_axis(
    lo: float, hi: float, cell: float, pin_lo: bool, pin_hi: bool,
) -> tuple[float, float, bool, bool]:
    # 핀이 걸린 쪽(자기 큐빗 몸체 방향)은 안쪽으로(ceil for lo, floor for hi) 스냅해
    # 그 큐빗의 경계 안 셀을 박스에서 뺀다 — 핀이 없으면 예전처럼 바깥으로(floor/ceil).
    raw_lo, raw_hi = lo, hi
    new_lo = (math.ceil(lo / cell) if pin_lo else math.floor(lo / cell)) * cell
    new_hi = (math.floor(hi / cell) if pin_hi else math.ceil(hi / cell)) * cell

    if new_hi <= new_lo and pin_lo and pin_hi:
        # 두 핀이 서로 충돌(두 큐빗 사이 간격이 흔히 1~2칸 사이라 드물지 않다) — 둘 다
        # 포기하지 않고, 단독 적용 시 더 넓게 남는 쪽 하나만 살린다.
        hi_only_lo = math.ceil(raw_hi / cell) * cell
        lo_only_hi = math.floor(raw_lo / cell) * cell
        width_only_lo = hi_only_lo - new_lo
        width_only_hi = new_hi - lo_only_hi
        if width_only_lo <= 0 and width_only_hi <= 0:
            new_lo, new_hi = math.floor(raw_lo / cell) * cell, math.ceil(raw_hi / cell) * cell
            pin_lo = pin_hi = False
        elif width_only_lo >= width_only_hi:
            new_hi, pin_hi = hi_only_lo, False
        else:
            new_lo, pin_lo = lo_only_hi, False

    if new_hi <= new_lo:
        # 핀 하나만 걸렸는데도 반대쪽 원래 경계와 충돌(두 포트가 1칸보다 가까움) —
        # 안전망으로 완전히 되돌린다.
        new_lo, new_hi = math.floor(raw_lo / cell) * cell, math.ceil(raw_hi / cell) * cell
        pin_lo = pin_hi = False

    return new_lo, new_hi, pin_lo, pin_hi


# region() 2)단계(최소 변 길이 cell 확보)의 핀 인식 버전 — 핀이 걸린 변은 밀 수 없으므로
# 남은 쪽만 넓힌다. 양쪽 다 핀이면(둘 다 못 밈) 핀을 포기하고 반씩 넓힌다.
def _pin_pad_axis(
    lo: float, hi: float, cell: float, pin_lo: bool, pin_hi: bool,
) -> tuple[float, float, bool, bool]:
    if hi - lo < cell:
        pad = cell - (hi - lo)
        if pin_lo and not pin_hi:
            hi += pad
        elif pin_hi and not pin_lo:
            lo -= pad
        else:
            lo, hi = lo - pad / 2.0, hi + pad / 2.0
            pin_lo = pin_hi = False
    return lo, hi, pin_lo, pin_hi


# 포트 p가 큐빗 q의 어느 축 경계에 정확히 닿아 있는지(법선 축)와, 그 축에서 q의 몸체가
# p 기준 어느 부호 방향에 있는지. qubit_port_toward()가 만드는 포트는 항상 두 축 중
# 하나만 ±(w/2 또는 h/2) 그대로이므로(그게 법선), 여기서 그 축을 되짚는다. 포트가
# 좌/우(x축) 또는 상/하(y축) 어느 변에도 있을 수 있으므로 둘 다 판별한다 — qubit_port_toward
# 가 이웃 방향에 따라 어느 변이든 고를 수 있게 만든 것과 짝을 이룬다.
def _port_edge_axis(p: tuple[float, float], q: "Qubit") -> tuple[str, float]:
    if abs(abs(p[0] - q.x) - q.w / 2.0) < 1e-6:
        return "x", (1.0 if p[0] < q.x else -1.0)
    return "y", (1.0 if p[1] < q.y else -1.0)

# cell_aabb(격자 셀 하나)가 큐빗 때문에 막혔는지 — core/globalplacement.py의
# _cell_blocked와 같은 판정(자기 큐빗은 coupler_own_port_cell로 배정 포트가 있는 셀만
# 예외, 제3자는 전부 막힘)을 GP의 segment_occupied 없이(아직 다른 커플러가 뭘 놨는지
# 모르는 단계이므로) 적용한 버전이다. region()이 이제 격자를 직접 세므로(아래
# _count_free_cells_in_box) GP와 정확히 같은 함수(coupler_own_port_cell) 하나로 판정을
# 통일한다 — 예전 _region_blocked_area는 면적을 연속값으로 근사했는데("포트 주변 한 칸
# 넓이"를 면적으로 흉내), 이제는 격자를 직접 알므로 근사가 필요 없다.
def _cell_blocked_by_qubits(
    coupler: "Coupler", qubits: dict[int, "Qubit"],
    ports: tuple[tuple[float, float], tuple[float, float]] | None,
    cell_aabb: tuple[float, float, float, float],
) -> bool:
    for qid, q in qubits.items():
        q_aabb = (q.x - q.w / 2.0, q.x + q.w / 2.0, q.y - q.h / 2.0, q.y + q.h / 2.0)
        if _rect_intersect(cell_aabb, q_aabb) is None:
            continue
        if coupler_own_port_cell(coupler, qid, ports, cell_aabb):
            continue
        return True
    return False


# [lo,hi] 축 하나에 "완전히" 포함되는 셀들의 인덱스 범위(양끝 포함, 없으면 None) —
# core/globalplacement.py의 _region_index_range와 정확히 같은 공식이다(중복은 부득이함:
# 그쪽을 여기서 import하면 계층이 뒤집힌다 — globalplacement.py가 이미 core.state를
# 가져오므로). **바깥쪽으로 안 부풀린다** — 이게 2026-09-17판(_snap_box_to_grid)과의
# 핵심 차이다. [lo,hi]가 격자에 안 맞으면 그 축 양끝의 부분 셀은 그냥 버려진다(예:
# lo=350, cell=200이면 i=1(200~400)은 [350,400)만 [lo,hi] 안이라 완전 포함이 아니므로
# i=2(400~600)부터 시작) — "포함 판정을 정확히 하면 스냅이 필요 없을 수 있다"는 제안을
# 그대로 구현한 것. 이것만 쓰면 2026-09-17 falcon 문제(정렬 손실로 46% 실패)가 그대로
# 재현된다는 것도 실측으로 확인했다(전용 스크립트: 407개 커플러 중 396개가 이걸로는
# 부족). 그래서 이 함수 단독이 아니라 region()의 반복 확장(아래)과 항상 같이 쓴다 —
# 모자라면 그만큼만 정확히 더 키우지, 애초에 여유분을 미리 얹어두지 않는다.
def _cell_index_range(lo: float, hi: float, cell: float) -> tuple[int, int] | None:
    i_lo = math.ceil(lo / cell - 1e-9)
    i_hi = math.floor(hi / cell + 1e-9) - 1
    if i_hi < i_lo:
        return None
    return int(i_lo), int(i_hi)


# box(연속좌표, 격자에 안 맞아도 됨) 안에서 "완전히 포함되면서" 큐빗에 안 막힌 격자 셀
# 개수. GP가 실제로 채울 수 있는 최대 세그먼트 수와 정확히 같은 값이다(다른 커플러의
# 세그먼트 점유는 아직 모르므로 그만큼의 낙관치 — GP 단계에서 경합으로 더 줄 수 있다,
# region()의 책임 밖). 반환값과 함께 실제 쓰인 셀 인덱스 범위(x_range, y_range)도
# 돌려준다 — region()이 "충분하다" 판정 시 그 범위 그대로를 박스 좌표로 쓰기 위함
# (다시 스냅할 필요가 없다 — 이미 셀 경계 그 자체다).
def _count_free_cells_in_box(
    coupler: "Coupler", box: tuple[float, float, float, float],
    qubits: dict[int, "Qubit"],
    ports: tuple[tuple[float, float], tuple[float, float]] | None, cell: float,
) -> tuple[int, tuple[int, int] | None, tuple[int, int] | None]:
    x0, x1, y0, y1 = box
    x_range = _cell_index_range(x0, x1, cell)
    y_range = _cell_index_range(y0, y1, cell)
    if x_range is None or y_range is None:
        return 0, x_range, y_range
    i0, i1 = x_range
    j0, j1 = y_range
    free = 0
    for i in range(i0, i1 + 1):
        for j in range(j0, j1 + 1):
            cell_aabb = (i * cell, (i + 1) * cell, j * cell, (j + 1) * cell)
            if not _cell_blocked_by_qubits(coupler, qubits, ports, cell_aabb):
                free += 1
    return free, x_range, y_range


@dataclass
class Qubit:
    id: int              # node id
    w:  float            # width
    h:  float            # height
    f:  float            # frequency ghz
    pad_inset_um: float   # 포트가 실제 닿는 변을 따라(접선 방향) 모서리에서 안쪽으로 들이는 여백 (um)
                          # — params.qubit_pad_inset_um. 2026-09-21 이전엔 "y-경계에서" 고정
                          # 값이었다(포트가 항상 좌/우 변에 있었으므로) — 아래 qubit_port_toward()
                          # docstring 참고.
    x:  float = 0.0       # x-coordinate(center)
    y:  float = 0.0       # y-coordinate(center)

    # ports 프로퍼티는 2026-09-21 제거했다 — 포트가 이제 "이 큐빗이 어느 이웃을 향하는가"에
    # 의존해(qubit_port_toward 참고) 큐빗 하나만으로는 계산할 수 없다(이웃 좌표가 필요).
    # 좌표는 core/floorplan.py의 build_ports(state 전체를 봄)가 FP 이후 한 번에 계산해
    # ChipState.ports에 담는다 — Qubit 자체는 순수하게 자기 자신의 물리 속성(w/h/f/위치)만
    # 갖는 얇은 데이터로 되돌아간다.

    def __str__(self):
        return (f"id={self.id:>3d} w={self.w:>6.2f} h={self.h:>6.2f} "
                f"f={self.f:>5.2f}GHz x={self.x:>8.2f} y={self.y:>8.2f}")


# 큐빗 중심에서 이웃 쪽으로 쏜 레이가 큐빗 사각형(w x h)과 만나는 경계점 — 01_mainref의
# utils/crosstalk.py::_pad_boundary_point를 이식했다(실측: 이 방식으로 포트를 두면
# 박스-자기큐빗 겹침이 819(202609-21 핀 적용 후 364)에서 0으로, "케이스 A"(포트가 실제
# 이웃 방향과 어긋나는 경우) 185건이 0건으로 사라진다 — 세션 보고서 참고). 고정 4슬롯
# (p_top_left 등, 항상 좌/우 변)과 달리 이 포트는 항상 "그 특정 이웃을 향한" 변 위에
# 있으므로, 박스가 그 방향으로 뻗어나가는 것 자체가 자기 큐빗 몸체에서 멀어지는 방향과
# 항상 일치한다 — case A가 구조적으로 불가능해지는 이유.
#
# pad_inset_um 반영: 고정 4슬롯 시절엔 포트가 항상 좌/우 변(off_x=w/2 정확히)에 있고
# pad_inset은 그 변을 "따라"(접선 방향, y축) 모서리에서 안쪽으로 들이는 용도였다(off_y =
# h/2 - inset) — 실제 qiskit-metal 커넥터 패드가 코너 정확히가 아니라 그 근처에 달리는
# 것을 반영). 레이-교점 방식은 어느 변에 닿을지 미리 알 수 없으므로, "닿은 변에 수직인
# 축(법선, 경계에 정확히 고정)"과 "그 변을 따르는 축(접선, 안쪽으로 inset)"을 매번
# 새로 판별해 같은 의미를 일반화한다: 법선 좌표는 항상 ±(w/2 또는 h/2) 그대로 두고,
# 접선 좌표만 중심 쪽으로 inset만큼 당긴다(코너를 완전히 넘어가지 않도록 접선 raw값의
# 절대값으로 클램프).
def qubit_port_toward(qubit: "Qubit", neighbor_x: float, neighbor_y: float) -> tuple[float, float]:
    dx, dy = neighbor_x - qubit.x, neighbor_y - qubit.y
    dist = math.hypot(dx, dy)
    if dist < 1e-9:
        # 이웃과 완전히 같은 좌표(퇴화) — Coupler.region()이 이 경우 자체를 ValueError로
        # 막으므로 여기까지 오면 안 되지만, 방어적으로 임의의 고정 방향(+x)을 쓴다.
        dx, dy, dist = 1.0, 0.0, 1.0
    dx_n, dy_n = dx / dist, dy / dist
    hw, hh = qubit.w / 2.0, qubit.h / 2.0
    t_x = hw / abs(dx_n) if abs(dx_n) > 1e-9 else math.inf
    t_y = hh / abs(dy_n) if abs(dy_n) > 1e-9 else math.inf
    inset = qubit.pad_inset_um
    if t_x <= t_y:
        # 좌/우 변에서 만남 — x가 법선(±hw 그대로), y가 접선(inset).
        rel_x = math.copysign(hw, dx_n)
        rel_y_raw = dy_n * t_x
        shrink = min(inset, abs(rel_y_raw))
        rel_y = rel_y_raw - math.copysign(shrink, rel_y_raw)
    else:
        # 상/하 변에서 만남 — y가 법선(±hh 그대로), x가 접선(inset).
        rel_y = math.copysign(hh, dy_n)
        rel_x_raw = dx_n * t_y
        shrink = min(inset, abs(rel_x_raw))
        rel_x = rel_x_raw - math.copysign(shrink, rel_x_raw)
    return (qubit.x + rel_x, qubit.y + rel_y)

@dataclass
class Segment:
    # QPlacer 모델(qplacer_bm/design_format.py의 wireblk)을 참고한 커플러 분할 단위.
    # GP가 std-cell처럼 개별 배치하는 대상 — 좌표(x, y)는 GP 이전엔 의미 없는 0.0.
    # coupler_id/size는 없다: 세그먼트는 항상 Coupler.segments 안에만 존재해 부모를 알고,
    # 같은 커플러의 세그먼트는 전부 같은 크기(Coupler.segment_size_um)라 따로 들고 다닐 필요가 없다.
    idx: int          # 커플러 내 순번 (0-based)
    x: float = 0.0    # x-coordinate(center) — GP가 채우기 전엔 0.0
    y: float = 0.0    # y-coordinate(center) — GP가 채우기 전엔 0.0

@dataclass
class Coupler:
    id: int      # edge id
    q1: int      # start qubit id
    q2: int      # end qubit id
    f:  float    # edge(coupler) frequency ghz

    # Coupler.l(반파장 CPW 공진기 길이) 공식에 쓰이는 기판 유전율.
    # params.substrate_epsilon_r 기본값(11.9)은 실리콘의 표준 문헌값이자
    # /home/LabMember/ngchoi/research/qiskit_metal/transmon_draw2.py:27의 calc_half_wave_length
    # 기본값과 동일하다 — 이전 v_phase=1.3e8(출처 불명)보다 유도 경로가 있다.
    epsilon_r: float

    # QPlacer 분할(segment) 모델의 세그먼트 정사각형 한 변 길이 (um) — QPlacer partition_size.
    # 물리적 유도 근거 없음. QPlacer 쪽(qplacer_bm/qplacement_param.py, benchmark_params.json)에도
    # 이 값을 설명하는 주석·문서가 없다 — 즉 QPlacer 코드 자체에도 근거가 없다.
    segment_size_um: float

    # 커플러 배선(meander CPW trace)의 실효 피치 (um). params.meander_spacing_um 기본값(120)은
    # qiskit-metal RouteMeander의 meander=dict(spacing='120um')에서 그대로 가져왔다
    # (transmon_draw2.py:227). 구성: trace_width=10um(223번 줄) + trace_gap 6um×2(224번 줄)
    # = 22um CPW 선폭 + fold간 간격을 합쳐 실효 피치 120um. 이전에 따로 있던
    # wire_width_um(100, QPlacer 유래)과 coupler_meander_pitch_um(40, 출처 불명)은
    # 둘 다 "배선이 차지하는 폭"이라는 같은 개념을 서로 다른 값으로 표현하고 있어서
    # 이 하나로 통합했다.
    meander_spacing_um: float

    # region()이 필요 셀 수(self.num_segments) 위에 얹는 여유 비율. params.box_slack_ratio의
    # description에 전체 경위(스윕 근거, 선택값)를 남겼다 — 요약: 2026-09-20 region()을
    # 연속 면적 기준에서 격자 셀 개수 기준으로 바꾸면서 이전에 격자 스냅이 "우연히"
    # 얹어주던 평균 1.82배 여유(GP의 다른 커플러와의 경합을 흡수하는 버퍼 역할을 겸했다)
    # 가 사라져 경합 실패가 23건→71건으로 뛰었다 — 그 버퍼를 우연한 부산물이 아니라
    # 이 명시적 파라미터로 대체한다.
    box_slack_ratio: float

    # GP가 채우기 전엔 비어 있다 — 이 리스트의 유무 자체가 "세그먼트 배치 이전/이후" 단계 구분이 된다.
    segments: list[Segment] = field(default_factory=list)

    # RT(core/router.py)가 확정하는 실제 배선 경로 — q1의 배정 포트에서 시작해 q2의 배정
    # 포트로 끝나는 절대좌표(um) 꺾은선(polyline). RT 이전엔 비어 있다(segments와 같은 패턴:
    # 리스트의 유무 자체가 "라우팅 이전/이후" 구분). RT가 경로를 못 찾은 커플러도 빈 리스트로
    # 남는다 — 그래서 "RT 이후 waypoints == []"는 실패를 뜻하고(segments가 있는 한 RT는
    # 항상 시도하므로 성공/실패가 모호하지 않다), Router.run()이 실패 사유를 별도로 기록한다.
    # 좌표 리스트(순수 (x,y) 튜플)로 둔 이유: (1) 렌더링(utils/rendering.py)이 바로 꺾은선으로
    # 그릴 수 있고, (2) qiskit-metal 변환 시 RoutePathfinder/RouteAnchors류(명시적 anchor
    # 좌표를 받는 라우터)에 그대로 넘길 수 있다 — RouteMeander(total_length 기반 자동 라우팅)
    # 대신 이쪽을 쓰는 게 이 필드의 존재 이유(RT가 직접 계산한 실제 경로)와 맞는다.
    waypoints: list[tuple[float, float]] = field(default_factory=list)

    @property
    def route_length(self) -> float:
        # waypoints 꺾은선의 총 길이(um) — l(목표 반파장 길이)과 비교해 길이 오차를 구할 때 쓴다.
        # 세그먼트 크기(l)와 달리 이건 RT가 실제로 만든 경로의 관측값이라 property로 둔다
        # (route_length == l이라는 보장은 없다 — RT가 미앤더로 채워도 회랑이 좁으면 못 채울 수
        # 있고, 그 오차 자체가 core/router.py의 검증 지표다).
        if len(self.waypoints) < 2:
            return 0.0
        return sum(
            math.hypot(self.waypoints[i + 1][0] - self.waypoints[i][0],
                       self.waypoints[i + 1][1] - self.waypoints[i][1])
            for i in range(len(self.waypoints) - 1)
        )

    @property
    def l(self) -> float:
        # λ/2 CPW 공진기 길이. 실효 유전율 근사 eps_eff = (eps_r + 1)/2 는 마이크로스트립
        # (microstrip) 라인용 근사라서 CPW(coplanar waveguide)엔 정확히 맞지 않는다 —
        # 실제 CPW eps_eff는 trace/gap 폭 비율에도 의존한다. 그래도 이전의 출처 불명
        # v_phase=1.3e8 상수보다는 유도 경로가 있다. eps_r=11.9, f=6.5GHz 기준으로
        # 이전 v_phase식 대비 길이가 약 10% 차이난다(짧아짐) — grid_25 비교 결과는
        # 이 변경의 커밋/보고 메시지에 남겨둔다.
        eps_eff = (self.epsilon_r + 1.0) / 2.0
        C_UM_PER_S = 2.99792458e14  # 빛의 속도 (um/s)
        return C_UM_PER_S / (2.0 * self.f * 1e9 * math.sqrt(eps_eff))

    @property
    def num_segments(self) -> int:
        # QPlacer의 num_poly = ceil(edge_wirelength * padding_size / partition_size**2)와 동일 구조.
        wire_area_um2 = self.l * self.meander_spacing_um
        return math.ceil(wire_area_um2 / (self.segment_size_um ** 2))

    def make_segments(self) -> list[Segment]:
        return [Segment(idx=i) for i in range(self.num_segments)]

    # 이 큐빗이 커플러 자신의 양 끝(q1 또는 q2)인지. 자기 큐빗이라는 사실 하나만으로 "장애물
    # 아님"을 보장하지는 않는다 — 2026-09-20까지는 그렇게 썼지만(자기 큐빗이면 몸체 전체를
    # 허용), 실측 결과 세그먼트가 포트 근처가 아니라 큐빗 중심 코앞(반폭 200um의 최대 91%인
    # 182um)까지, 칩당 1~28개는 몸체 안에 완전히 파묻히는 것으로 나왔다
    # (docs/20260920_segment_model_review.md 발견 3). 물리적으로 큐빗 몸체는 금속(포켓/패드/
    # 조셉슨 접합)이라 배선이 그 위를 지나갈 수 없다 — 커플러가 실제로 필요로 하는 건 배정된
    # 포트 한 점(패드) 접촉뿐이다. 그래서 "장애물 아님" 판정은 이제 이 메서드(자기 큐빗인가)
    # 단독이 아니라 coupler_own_port_cell()(자기 큐빗이면서 배정 포트 근처인가)로 대체됐다 —
    # 이 메서드 자체는 coupler_own_port_cell() 내부에서, 그리고 다른 "이 큐빗이 내 소유인가"
    # 판정이 필요한 곳(예: hotspot 계산에서 직결 커플러 제외)에서 여전히 쓰인다.
    def is_own_qubit(self, qubit_id: int) -> bool:
        return qubit_id == self.q1 or qubit_id == self.q2

    # boundary()는 제거했다 — sqrt(l * p) 공식이 두 큐빗 사이 거리(span)를 무시하고,
    # 필요 배선 면적(l * meander_spacing_um)보다 훨씬 작은 박스를 만들어 배치 전/후 어느
    # 용도로도 쓸 수 없었다. "배치 전 후보 영역"과 "배치 후 실제 footprint"는 서로 다른
    # 개념이라 region()/bbox()로 분리한다.

    def region(
        self, q1: Qubit, q2: Qubit,
        p1: tuple[float, float] | None, p2: tuple[float, float] | None,
        qubits: dict[int, Qubit], chip_width: float, chip_height: float,
    ) -> tuple[float, float, float, float] | None:
        # 배치 전 제약: GP가 이 커플러의 segments를 채울 때, 그 안에서만 배치하도록 쓸
        # 후보 영역이다. DRC는 bbox()(AABB)를 쓰므로 여기서도 AABB로 유지한다 — 회전
        # 사각형은 겹침 판정이 복잡해져 지금 단계에서 도입하지 않는다.
        #
        # p1/p2는 2026-09-21부터 포트 "이름"이 아니라 좌표를 직접 받는다(core/floorplan.py의
        # build_ports가 qubit_port_toward로 미리 계산해 둔 것) — 포트가 이제 "이 큐빗이
        # 어느 이웃을 향하는가"에 의존해 이름-슬롯 하나로 고정할 수 없어졌기 때문이다.
        # None이면(호출부가 아직 포트를 못 정한 경우) 큐빗 중심으로 대신한다 — 이 경우
        # 반환된 영역이 실제 포트를 덮는다는 보장이 없다.
        #
        # qubits(칩 전체 큐빗)/chip_width/chip_height는 2026-09-20 own-qubit 예외를
        # "포트 셀만"으로 좁힌 뒤 추가됐다 — 아래 3) 참고. None을 반환할 수 있게 된 것도
        # 이때부터다(그 전엔 항상 박스를 돌려줬다): 가용 면적 확장이 die 밖으로 나가야
        # 하거나 반복 한도 안에 못 끝나면, 이 커플러는 이 위치에서 애초에 배선이 안
        # 들어간다는 뜻이라 None으로 실패를 알린다 — 호출부(core/floorplan.py)가 그대로
        # coupler_regions[key]=None으로 남기면, GP의 기존 "박스 없음" 실패 경로
        # (core/globalplacement.py의 _place_chain, `box is None: return None`)를 그대로
        # 타서 이 커플러 하나만 실패로 세고 칩 전체는 계속 진행한다.
        has_p1, has_p2 = p1 is not None, p2 is not None
        if p1 is None:
            p1 = (q1.x, q1.y)
        if p2 is None:
            p2 = (q2.x, q2.y)

        # 정확한 float 등가(p1 == p2) 대신 거리 임계값을 쓴다 — FP가 내놓는 좌표는 반복
        # 최적화 결과라 "같은 점"이어도 부동소수 잔차(예: 0.001um)가 남을 수 있고, 그 정도
        # 거리로 아래 3)에서 required_wire_area / w(또는 h)를 나누면 scale이 폭발한다.
        if math.hypot(p2[0] - p1[0], p2[1] - p1[1]) < 1e-6:
            # region()은 FP 이후(적어도 두 큐빗/포트가 서로 다른 위치일 때)에만 유효하다.
            # 예전엔 이 경우 조용히 segment_size_um을 대신 써서 의미 없는 박스를 반환했다 —
            # DRC에서 num_unplaced_couplers로 미검사를 드러낸 것과 같은 원칙으로, 여기서도
            # 조용히 넘어가지 않고 즉시 실패시킨다.
            raise ValueError(
                f"Coupler(id={self.id}, q1={self.q1}, q2={self.q2}): "
                f"두 포트가 사실상 같은 점입니다 {p1} ~ {p2} — region()은 FP 이후에만 유효합니다."
            )

        # 1) 두 포트를 반드시 덮는 최소 AABB (대각선이어도 x/y 각 축을 직접 min/max로 잡으므로
        #    항상 두 점을 포함한다 — 이전 버전은 abs(dy) > abs(dx)로 축 하나만 골라 나머지
        #    축에서 포트를 놓치는 버그가 있었다).
        x0, x1 = min(p1[0], p2[0]), max(p1[0], p2[0])
        y0, y1 = min(p1[1], p2[1]), max(p1[1], p2[1])
        cell = self.segment_size_um

        # 핀(pin): 포트가 자기 큐빗의 경계(qubit_port_toward가 만드는 "법선" 축)에
        # 있으므로, 그 포트가 세우는 축의 경계값을 바깥으로(floor/ceil) 넓히면 자기
        # 큐빗 몸체로 들어갈 수도, 안 들어갈 수도 있다 — 두 경우가 갈린다:
        #   케이스 B(핀 적용): 큐빗 몸체가 포트에서 "박스가 뻗어가는 반대 방향"에 있으면
        #   원래 박스 구간은 큐빗을 안 건드린다 — floor/ceil 격자 스냅만이 최대 1칸(cell)
        #   침범을 만든다. 스냅 방향을 뒤집어(바깥 대신 안쪽으로) 완전히 막는다 — 그
        #   대가로 포트가 걸친 셀 자체는 박스에서 빠진다(GP가 _resolve_endpoint로 한 칸
        #   바깥의 빈 셀에서 시작해 lead로 포트에 잇는다).
        #   케이스 A(핀 안 함): 큐빗 몸체가 포트에서 "박스가 뻗어가는 같은 방향"에 있으면,
        #   박스는 스냅 여부와 무관하게 이미 큐빗의 그 축 전체 폭을 포함한다 — floor/ceil
        #   방향을 바꿔도 해결되지 않는다. 이 함수의 책임 범위 밖으로 남긴다.
        #
        # 2026-09-21에 이웃 방향 포트(qubit_port_toward)로 바꾸면서 이 핀이 어느 축(x
        # 또는 y)에 걸리는지도 포트마다 다시 판별해야 하게 됐다 — 고정 4슬롯 시절엔 포트가
        # 항상 좌/우 변(x축)에 있어 케이스 A가 819개 중 185개(x축 전용)였는데, 이웃 방향
        # 포트는 그 자체로 케이스 A가 구조적으로 거의 발생하지 않는다(포트가 항상 실제
        # 이웃을 향한 변에 있으므로 "박스가 뻗는 방향 = 큐빗 몸체에서 멀어지는 방향"이
        # 자동으로 일치) — 실측(세션 보고서): 케이스 A 185→0. 그래도 안전망으로 판별
        # 로직 자체는 그대로 둔다(포트가 예외적으로 안 맞는 경우를 대비).
        pin_x0 = pin_x1 = pin_y0 = pin_y1 = False
        for p_own, q_own, p_other, has_port in (
            (p1, q1, p2, has_p1), (p2, q2, p1, has_p2),
        ):
            if not has_port:
                continue
            axis, qubit_dir = _port_edge_axis(p_own, q_own)
            if axis == "x":
                sets_lo = p_own[0] <= p_other[0]
            else:
                sets_lo = p_own[1] <= p_other[1]
            box_dir = 1.0 if sets_lo else -1.0
            if qubit_dir == box_dir:
                continue  # 케이스 A -- 핀으로 못 막음, 그대로 둠
            if axis == "x":
                pin_x0, pin_x1 = pin_x0 or sets_lo, pin_x1 or not sets_lo
            else:
                pin_y0, pin_y1 = pin_y0 or sets_lo, pin_y1 or not sets_lo

        # 1b) 격자 스냅 — 축마다 같은 핀 인식 로직(_pin_snap_axis)을 쓴다. 2026-09-20
        # (연속 면적 대신 격자 셀 개수를 직접 세도록 바꾼 버전) 이후로 3)단계가 "완전히
        # 포함된 셀만" 센다(_count_free_cells_in_box, _cell_index_range) — 안 핀된 축은
        # 예전처럼 바깥으로 최대 반 칸 넓혀 포트가 속한 셀이 박스에서 안 잘려나가게 한다
        # (실측(grid_25 커플러(16,17)): 이 스냅이 없으면 포트 좌표가 박스의 원 경계가 돼
        # 그 포트의 셀이 박스 밖으로 밀려났고, 40개 중 21개가 이 버그로 실패했었다).
        x0, x1, pin_x0, pin_x1 = _pin_snap_axis(x0, x1, cell, pin_x0, pin_x1)
        y0, y1, pin_y0, pin_y1 = _pin_snap_axis(y0, y1, cell, pin_y0, pin_y1)

        # 2) 각 변 하한: segment_size_um — 그보다 좁으면 세그먼트(정사각형)가 물리적으로
        #    안 들어간다. 핀이 걸린 변은 밀 수 없으므로 남은 쪽만 넓힌다.
        x0, x1, pin_x0, pin_x1 = _pin_pad_axis(x0, x1, cell, pin_x0, pin_x1)
        y0, y1, pin_y0, pin_y1 = _pin_pad_axis(y0, y1, cell, pin_y0, pin_y1)

        # 3) 셀 하한: GP가 실제로 놓아야 할 세그먼트 개수(self.num_segments) 이상의 '가용
        #    격자 셀'이 되도록 반복 확장한다. 매 반복 지금 박스([x0,x1]×[y0,y1], 격자에
        #    안 맞아도 됨) 안에 "완전히 포함되는" 셀만 세고(_count_free_cells_in_box,
        #    core/globalplacement.py의 _region_index_range와 같은 공식 — 격자에 안 맞는
        #    가장자리는 버려진다, 안 부풀린다), 그 개수가 충분하면 딱 그 셀들의 경계를
        #    박스로 반환한다. 모자라면 단변을 cell 한 칸만큼만 넓혀 다시 센다 — "얼마나
        #    부족한지"를 셀 단위로 직접 재므로, 필요한 만큼만 넓어지고 미리 여유를
        #    얹어두지 않는다(장변까지 늘리면 포트 바깥으로 삐져나가므로 단변만).
        #
        #    목표 셀 수는 self.num_segments가 아니라 그 위에 box_slack_ratio만큼 얹은
        #    값이다(params.box_slack_ratio의 description에 스윕 근거).
        #    required_segments가 0이면(짧은 커플러) 1+ratio를 곱해도 ceil(0)=0이라
        #    영향이 없다.
        required_segments = math.ceil(self.num_segments * (1.0 + self.box_slack_ratio))
        ports = (p1, p2)
        for _ in range(_REGION_GROW_MAX_ITERS):
            n_free, x_range, y_range = _count_free_cells_in_box(
                self, (x0, x1, y0, y1), qubits, ports, cell,
            )
            if n_free >= required_segments and x_range is not None and y_range is not None:
                i0, i1 = x_range
                j0, j1 = y_range
                return (i0 * cell, (i1 + 1) * cell, j0 * cell, (j1 + 1) * cell)
            w, h = x1 - x0, y1 - y0
            # 짧은 변이어도 양쪽 다 핀이면(자기 큐빗 몸체 앞에서 더 못 밈) 그 축은 못
            # 키우고 다른 축을 키운다 -- 핀은 위 1)/2)에서 이미 자기 큐빗 겹침을 막았으므로,
            # 성장이 그 핀을 다시 깨서(=큐빗 몸체로 도로 들어가서) 되돌리면 안 된다.
            x_growable = not (pin_x0 and pin_x1)
            y_growable = not (pin_y0 and pin_y1)
            if not x_growable and not y_growable:
                return None  # 두 축 다 못 키움(극히 좁은 경우) -- 판단 못 하면 안전하게 실패.
            if (w >= h and y_growable) or not x_growable:
                new_y0 = y0 if pin_y0 else y0 - cell / 2.0
                new_y1 = y1 if pin_y1 else y1 + cell / 2.0
                if new_y0 < 0.0 or new_y1 > chip_height:
                    return None
                y0, y1 = new_y0, new_y1
            else:
                new_x0 = x0 if pin_x0 else x0 - cell / 2.0
                new_x1 = x1 if pin_x1 else x1 + cell / 2.0
                if new_x0 < 0.0 or new_x1 > chip_width:
                    return None
                x0, x1 = new_x0, new_x1

        return None

    def bbox(self) -> tuple[float, float, float, float] | None:
        # 배치 후 관찰: 실제 배치된 세그먼트들의 bounding box. segments가 비어 있으면(GP 이전) None.
        if not self.segments:
            return None
        half = self.segment_size_um / 2.0
        x0 = min(s.x - half for s in self.segments)
        x1 = max(s.x + half for s in self.segments)
        y0 = min(s.y - half for s in self.segments)
        y1 = max(s.y + half for s in self.segments)
        return (x0, x1, y0, y1)

    def __str__(self):
        return (f"id={self.id:>3d} q1={self.q1:>3d} q2={self.q2:>3d} "
                f"f={self.f:>5.2f}GHz segments={len(self.segments)}")


# 큐빗-세그먼트 "장애물 아님" 판정의 유일한 정의 (2026-09-20, is_own_qubit 단독 사용을
# 대체, 같은 날 점-포함 판정으로 재수정 — 아래 "판정 방식" 항목 참고). aabb(세그먼트
# 하나 또는 GP 격자 셀 하나의 (x0,x1,y0,y1))가 qubit_id의 소유이고, 그 큐빗에서 coupler에게
# 배정된 포트 점(px,py)을 담고 있으면 정상(장애물 아님) — 그 밖의 자기 큐빗 겹침(몸체
# 내부)은 이제 장애물/위반으로 센다.
#
# 이전(포트 무관, 자기 큐빗이면 몸체 전체 허용) 대비 좁힌 이유: 큐빗 몸체는 실제로는
# 금속(TransmonPocket 포켓/패드/조셉슨 접합)이라 배선이 그 위를 지날 수 없다 — 커플러가
# 필요로 하는 건 배정 포트 한 점(패드) 접촉뿐이다. 실측(GP 산출물, DP 이전) 결과 이전
# 규칙 아래서 세그먼트가 포트 근처가 아니라 큐빗 중심 코앞까지(반폭 200um의 최대 91%인
# 182um) 파고들었고, 칩당 1~28개는 몸체 안에 완전히 들어가 있었다
# (docs/20260920_segment_model_review.md 발견 3).
#
# 판정 방식(점-포함, cell_size 매개변수 제거): 처음 버전은 "aabb가 포트 중심의
# cell_size 정사각형과 겹치는가"(사각형-사각형 AABB 겹침)였다. 그런데 그 포트-중심
# 정사각형은 GP의 전역 격자에 정렬돼 있지 않다(포트 좌표가 FP 최적화 결과라 격자
# 배수가 아님) — 격자 정렬 사각형(aabb)과 비정렬 사각형이 겹치는지를 물으면, 포트를
# 담은 격자 셀뿐 아니라 그 옆의(때로는 대각선의) 인접 셀까지 "겹친다"고 잡힌다. 실측
# (6칩, 이 수정 직전): 자기 큐빗과 겹치는 것으로 "허용된" 세그먼트 중 30~46%가 실제로는
# 포트를 담은 셀이 아니라 이 기하학적 누출로 통과한 인접 셀이었고, 그중 일부는 배정
# 포트에서 최대 277um(셀 크기 200um의 1.4배) 떨어진, 큐빗 몸체 안에 완전히 파묻힌
# 셀이었다 — "포트 셀 하나만 허용"이라는 애초 의도를 사각형 겹침 판정이 지키지 못한
# 것이다. 점-포함 판정("aabb가 포트 점을 담는가")으로 바꾸면 포트 좌표가 어디에 있든
# 항상 정확히 그 점을 담은 격자 셀 하나만 골라내므로(포트 자체가 큐빗 경계(±w/2)
# 위에 있어 그 셀은 거의 항상 경계에 걸친 셀이 된다) 이 누출이 구조적으로 불가능해진다.
# 그래서 cell_size 매개변수 자체를 없앴다(더 이상 쓰이지 않음).
#
# core/globalplacement.py(GP의 장애물 판정, `_cell_blocked`), func/compute.py
# (`compute_drc`의 qc 판정), core/legalization.py(`_find_qc_overlaps`) 셋 다 이 함수
# 하나만 쓴다 — 이 판정을 여러 곳에서 각자 구현하면 GP가 "놓을 수 있다"고 판단한 배치를
# DRC가 "위반"으로 다시 잡아내는(또는 그 반대) 모순이 생긴다(같은 부류의 실수를
# core/legalization.py 모듈 docstring이 AABB_EPS_UM에 대해 이미 경고한 바 있다).
#
# ports가 None이면(이론상 GP/LG 호출부에서 세그먼트가 있는 커플러는 항상 배정도
# 있으므로 도달 불가능, 방어적 분기) 포트를 알 수 없으니 항상 장애물로 취급한다 — "판단
# 못 하면 안전하게 막는다"는 이 저장소의 일관된 원칙(예: Coupler.region()의 방어적
# ValueError)과 같다.
#
# 반개구간([x0,x1), [y0,y1))으로 담는다 — core/globalplacement.py의 _qubit_owner_cells가
# 격자 셀을 같은 방식(floor 기반, 오른쪽/위쪽 경계 미포함)으로 다루는 것과 규칙을
# 맞춘 것이다. 포트 좌표가 부동소수 최적화 결과라 셀 경계선에 정확히 걸칠 확률은
# 사실상 0이므로 이 선택이 실측 결과에 영향을 주지는 않는다.
#
# 2026-09-21: ports 매개변수가 "포트 이름"(assignment) 대신 좌표 쌍(p1,p2) 자체를 받는다
# — 포트가 이제 이웃 방향으로 계산되고 이름이 없어졌으므로(qubit_port_toward,
# core/floorplan.py의 build_ports) 이름으로 Qubit.ports를 조회할 이유도 없어졌다.
# qubit 매개변수 자체가 필요 없어져 시그니처에서 뺐다.
def coupler_own_port_cell(
    coupler: "Coupler", qubit_id: int,
    ports: tuple[tuple[float, float], tuple[float, float]] | None,
    aabb: tuple[float, float, float, float],
) -> bool:
    if not coupler.is_own_qubit(qubit_id):
        return False
    if ports is None:
        return False
    px, py = ports[0] if qubit_id == coupler.q1 else ports[1]
    ax0, ax1, ay0, ay1 = aabb
    return ax0 <= px < ax1 and ay0 <= py < ay1


@dataclass
class ChipState:
    processor_name: str
    num_qubits: int
    cmap: tuple  # coupling map: ((q1, q2), ...) — (min, max)로 정규화·중복제거된 상태
    chip_width: float
    chip_height: float

    qubits: dict[int, Qubit] = field(default_factory=dict)
    couplers: dict[tuple[int, int], Coupler] = field(default_factory=dict)

    # 커플러 → (q1 쪽 포트 좌표, q2 쪽 포트 좌표). 2026-09-21 이전엔 port_assignment
    # (포트 "이름" 쌍, p_top_left 등 고정 4슬롯)이었다 — 포트가 이웃 방향으로 계산되며
    # (qubit_port_toward, core/floorplan.py의 build_ports) 이름 있는 슬롯 자체가
    # 없어졌으므로 좌표를 직접 담는다. 계산 로직은 FP에서 채운다 — 여기선 필드만 정의.
    ports: dict[tuple[int, int], tuple[tuple[float, float], tuple[float, float]]] = field(default_factory=dict)

    # 커플러 → FP가 확정한 배치 후보 영역(Coupler.region() 결과, AABB (x0,x1,y0,y1)).
    # ports와 같은 패턴(값 자체는 Coupler가 아니라 ChipState의 dict에 둠)을
    # 따른다 — 둘 다 "FP가 큐빗/포트를 보고 계산해 커플러별로 확정하는 사실"이라는 같은
    # 종류의 데이터라, Coupler에 새 필드를 얹기보다 기존 관례를 그대로 잇는 게 일관적이다.
    # GP는 이 안에서만 세그먼트를 배치한다(박스 밖으로 넓히지 않음, core/globalplacement.py
    # 참고) — 즉 여기 저장된 박스가 GP의 하드 탐색 범위 자체가 된다(FP 이전엔 하드 제약
    # 아니었던 region()의 "힌트"를, FP가 확정한 순간부터 하드 제약으로 승격시키는 셈).
    # 값이 None이면(2026-09-20) region()의 가용 면적 확장이 수렴하지 못했다는 뜻 —
    # Coupler.region() docstring 참고. GP는 이를 "박스 없음" 실패로 처리한다.
    coupler_regions: dict[tuple[int, int], tuple[float, float, float, float] | None] = field(default_factory=dict)

    # 아래 3개는 FP(core/floorplan.py)가 확정해 채우는 필드 — GP가 참조한다(예: crossover로
    # 잘려나간 엣지는 GP/routing이 air bridge 등 다른 방식으로 처리해야 함을 알아야 한다).
    # FP 이전(예: init_state)엔 기본값 그대로다.
    embedding: object = None     # networkx.PlanarEmbedding | None — FP가 확정한 평면 조합적 임베딩
    cut_edges: tuple = ()        # crossover 모드에서 평면성 확보를 위해 잘라낸 엣지들 (air bridge 필요)
    fallback_used: bool = False  # FP가 주 경로(휴리스틱) 대신 결정론적 보장 폴백을 썼는지 여부

    @classmethod
    def processor_config(cls, config: dict, params) -> "ChipState":
        num_qubits = config["num_qubits"]
        # chip_width = config.get("chip_width_um") or params.chip_width
        # chip_height = config.get("chip_height_um") or chip_width
        chip_width = config.get("chip_width_um")
        chip_height = config.get("chip_height_um")
        if chip_width is None or chip_height is None:
            raise ValueError(
                f"{config.get('processor', '?')}: chip_width_um/chip_height_um이 input JSON에 "
                f"없습니다 (source={config.get('source_path')})"
            )

        freq_ghz = config["freq_ghz"]
        qubits = {
            i: Qubit(id=i, w=params.qubit_width, h=params.qubit_height, f=freq_ghz[i],
                     pad_inset_um=params.qubit_pad_inset_um)
            for i in range(num_qubits)
        }

        # edge_freq는 (min, max)로 정규화된 키를 쓰는데 cmap은 원본 순서 그대로였다.
        # cmap에 (5,3)처럼 뒤집힌 엣지가 있으면 KeyError, 양방향이 다 있으면 커플러가
        # 중복 생성됐다 — cmap 자체를 정규화·중복제거해서 couplers 키와 항상 일치시킨다.
        cmap = tuple(sorted(set((min(u, v), max(u, v)) for u, v in config["coupling_map"])))
        edge_freq = {(min(u, v), max(u, v)): f for u, v, f in config["edge_freq_ghz"]}
        couplers = {
            (q1, q2): Coupler(
                id=idx, q1=q1, q2=q2, f=edge_freq[(q1, q2)],
                epsilon_r=params.substrate_epsilon_r,
                meander_spacing_um=params.meander_spacing_um,
                segment_size_um=params.segment_size_um,
                box_slack_ratio=params.box_slack_ratio,
            )
            for idx, (q1, q2) in enumerate(cmap)
        }

        return cls(
            processor_name  = config["processor"],
            num_qubits      = num_qubits,
            cmap            = cmap,
            chip_width      = chip_width,
            chip_height     = chip_height,
            qubits          = qubits,
            couplers        = couplers
        )