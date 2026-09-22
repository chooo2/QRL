import math
from dataclasses import dataclass, field

AABB_EPS_UM = 1e-6
_REGION_GROW_MAX_ITERS = 40
_REGION_MAX_ASPECT = 4.0
COUPLER_PHASE_VELOCITY_UM_PER_S = 1.3e14  # QPlacer식의 v0 ~= 1.3e8 m/s
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


# Coupler.region()의 축별 "핀"(pin) 헬퍼
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

def _pin_grow_axis(
    lo: float, hi: float, amount: float, pin_lo: bool, pin_hi: bool,
    lower_bound: float, upper_bound: float,
) -> tuple[float, float] | None:
    if pin_lo and pin_hi:
        return None
    if pin_lo:
        new_lo, new_hi = lo, hi + amount
    elif pin_hi:
        new_lo, new_hi = lo - amount, hi
    else:
        new_lo, new_hi = lo - amount / 2.0, hi + amount / 2.0
        if new_lo < lower_bound:
            shift = lower_bound - new_lo
            new_lo += shift
            new_hi += shift
        if new_hi > upper_bound:
            shift = new_hi - upper_bound
            new_lo -= shift
            new_hi -= shift
    if new_lo < lower_bound or new_hi > upper_bound:
        return None
    return new_lo, new_hi

def _grow_axis_centered(
    lo: float, hi: float, amount: float,
    lower_bound: float, upper_bound: float,
) -> tuple[float, float] | None:
    new_lo, new_hi = lo - amount / 2.0, hi + amount / 2.0
    if new_lo < lower_bound:
        shift = lower_bound - new_lo
        new_lo += shift
        new_hi += shift
    if new_hi > upper_bound:
        shift = new_hi - upper_bound
        new_lo -= shift
        new_hi -= shift
    if new_lo < lower_bound or new_hi > upper_bound:
        return None
    return new_lo, new_hi

def _port_edge_axis(p: tuple[float, float], q: "Qubit") -> tuple[str, float]:
    if abs(abs(p[0] - q.x) - q.w / 2.0) < 1e-6:
        return "x", (1.0 if p[0] < q.x else -1.0)
    return "y", (1.0 if p[1] < q.y else -1.0)

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

def _cell_index_range(lo: float, hi: float, cell: float) -> tuple[int, int] | None:
    i_lo = math.ceil(lo / cell - 1e-9)
    i_hi = math.floor(hi / cell + 1e-9) - 1
    if i_hi < i_lo:
        return None
    return int(i_lo), int(i_hi)


# box(연속좌표, 격자에 안 맞아도 됨) 안에서 "완전히 포함되면서" 큐빗에 안 막힌 격자 셀
# 개수. GP가 실제로 채울 수 있는 최대 세그먼트 수와 정확히 같은 값
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
    x:  float = 0.0       # x-coordinate(center)
    y:  float = 0.0       # y-coordinate(center)

    def __str__(self):
        return (f"id={self.id:>3d} w={self.w:>6.2f} h={self.h:>6.2f} f={self.f:>5.2f}GHz x={self.x:>8.2f} y={self.y:>8.2f}")

def qubit_port_toward(qubit: "Qubit", neighbor_x: float, neighbor_y: float) -> tuple[float, float]:
    dx, dy = neighbor_x - qubit.x, neighbor_y - qubit.y
    hw, hh = qubit.w / 2.0, qubit.h / 2.0
    y_margin = max(0.0, hh - qubit.pad_inset_um)
    x = qubit.x + (hw if dx >= 0.0 else -hw)
    y = qubit.y + (y_margin if dy >= 0.0 else -y_margin)
    return (x, y)

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
    epsilon_r: float
    segment_size_um: float
    meander_spacing_um: float
    box_slack_ratio: float
    region_max_aspect: float
    segments: list[Segment] = field(default_factory=list)
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
        # QPlacer식: f = v0 / (2L), v0 ~= 1.3e8 m/s. self.f는 GHz이므로 1e9를 곱한다.
        return COUPLER_PHASE_VELOCITY_UM_PER_S / (2.0 * self.f * 1e9)

    @property
    def num_segments(self) -> int:
        return math.ceil(self.l / self.segment_size_um)

    @property
    def boundary_grid_shape(self) -> tuple[int, int]:
        n = self.num_segments

        best_nx = 1
        best_ny = n
        best_score = None

        for nx in range(1, n + 1):
            ny = math.ceil(n / nx)

            aspect_diff = abs(nx - ny)
            capacity = nx * ny

            score = (aspect_diff, capacity)

            if best_score is None or score < best_score:
                best_score = score
                best_nx = nx
                best_ny = ny

        return best_nx, best_ny

    @property
    def boundary_size(self) -> tuple[float, float]:
        nx, ny = self.boundary_grid_shape

        return (nx * self.segment_size_um, ny * self.segment_size_um)

    def make_segments(self) -> list[Segment]:
        return [Segment(idx=i) for i in range(self.num_segments)]

    def is_own_qubit(self, qubit_id: int) -> bool:
        return qubit_id == self.q1 or qubit_id == self.q2


    def region(
        self, q1: Qubit, q2: Qubit,
        p1: tuple[float, float] | None, p2: tuple[float, float] | None,
        qubits: dict[int, Qubit], chip_width: float, chip_height: float,
    ) -> tuple[float, float, float, float] | None:
        has_p1, has_p2 = p1 is not None, p2 is not None
        if p1 is None:
            p1 = (q1.x, q1.y)
        if p2 is None:
            p2 = (q2.x, q2.y)


        if math.hypot(p2[0] - p1[0], p2[1] - p1[1]) < 1e-6:

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
        target_nx, target_ny = self.boundary_grid_shape
        ports = (p1, p2)
        for _ in range(_REGION_GROW_MAX_ITERS):
            n_free, x_range, y_range = _count_free_cells_in_box(
                self, (x0, x1, y0, y1), qubits, ports, cell,
            )
            if x_range is not None and y_range is not None:
                nx = x_range[1] - x_range[0] + 1
                ny = y_range[1] - y_range[0] + 1
                aspect = max(nx / max(ny, 1), ny / max(nx, 1))
                shape_ok = (
                    nx >= target_nx
                    and ny >= target_ny
                    and aspect <= self.region_max_aspect
                )
            else:
                nx = ny = 0
                shape_ok = False
            if n_free >= required_segments and x_range is not None and y_range is not None:
                i0, i1 = x_range
                j0, j1 = y_range
                candidate = (i0 * cell, (i1 + 1) * cell, j0 * cell, (j1 + 1) * cell)
                if shape_ok:
                    return candidate

                if nx < target_nx and not (pin_x0 and pin_x1):
                    needed = max(1, target_nx - nx) * cell
                    grown = _pin_grow_axis(x0, x1, needed, pin_x0, pin_x1, 0.0, chip_width)
                    if grown is not None:
                        x0, x1 = grown
                        continue
                if ny < target_ny and not (pin_y0 and pin_y1):
                    needed = max(1, target_ny - ny) * cell
                    grown = _pin_grow_axis(y0, y1, needed, pin_y0, pin_y1, 0.0, chip_height)
                    if grown is not None:
                        y0, y1 = grown
                        continue
                if nx > ny and aspect > self.region_max_aspect:
                    desired_ny = math.ceil(nx / self.region_max_aspect)
                    grown = _pin_grow_axis(y0, y1, max(1, desired_ny - ny) * cell, pin_y0, pin_y1, 0.0, chip_height)
                    if grown is None:
                        grown = _grow_axis_centered(
                            y0, y1, max(1, desired_ny - ny) * cell,
                            0.0, chip_height,
                        )
                    if grown is not None:
                        y0, y1 = grown
                        continue
                if ny > nx and aspect > self.region_max_aspect:
                    desired_nx = math.ceil(ny / self.region_max_aspect)
                    grown = _pin_grow_axis(x0, x1, max(1, desired_nx - nx) * cell, pin_x0, pin_x1, 0.0, chip_width)
                    if grown is None:
                        grown = _grow_axis_centered(
                            x0, x1, max(1, desired_nx - nx) * cell,
                            0.0, chip_width,
                        )
                    if grown is not None:
                        x0, x1 = grown
                        continue
                if aspect > self.region_max_aspect:
                    if nx > ny:
                        grown = _pin_grow_axis(y0, y1, cell, pin_y0, pin_y1, 0.0, chip_height)
                        if grown is None:
                            grown = _grow_axis_centered(y0, y1, cell, 0.0, chip_height)
                        if grown is not None:
                            y0, y1 = grown
                            continue
                    if ny > nx:
                        grown = _pin_grow_axis(x0, x1, cell, pin_x0, pin_x1, 0.0, chip_width)
                        if grown is None:
                            grown = _grow_axis_centered(x0, x1, cell, 0.0, chip_width)
                        if grown is not None:
                            x0, x1 = grown
                            continue
                return candidate

            # 짧은 변이어도 양쪽 다 핀이면(자기 큐빗 몸체 앞에서 더 못 밈) 그 축은 못
            # 키우고 다른 축을 키운다 -- 핀은 위 1)/2)에서 이미 자기 큐빗 겹침을 막았으므로,
            # 성장이 그 핀을 다시 깨서(=큐빗 몸체로 도로 들어가서) 되돌리면 안 된다.
            x_growable = not (pin_x0 and pin_x1)
            y_growable = not (pin_y0 and pin_y1)
            if not x_growable and not y_growable:
                return None

            grow_x = False
            grow_y = False
            if nx < target_nx and x_growable:
                grow_x = True
            if ny < target_ny and y_growable:
                grow_y = True
            if not grow_x and not grow_y:
                if nx > ny and y_growable:
                    grow_y = True
                elif ny > nx and x_growable:
                    grow_x = True
                elif x_growable:
                    grow_x = True
                else:
                    grow_y = True

            if grow_y and not grow_x:
                grown = _pin_grow_axis(y0, y1, cell, pin_y0, pin_y1, 0.0, chip_height)
                if grown is None:
                    if not x_growable:
                        return None
                    grown = _pin_grow_axis(x0, x1, cell, pin_x0, pin_x1, 0.0, chip_width)
                    if grown is None:
                        return None
                    x0, x1 = grown
                else:
                    y0, y1 = grown
            else:
                grown = _pin_grow_axis(x0, x1, cell, pin_x0, pin_x1, 0.0, chip_width)
                if grown is None:
                    if not y_growable:
                        return None
                    grown = _pin_grow_axis(y0, y1, cell, pin_y0, pin_y1, 0.0, chip_height)
                    if grown is None:
                        return None
                    y0, y1 = grown
                else:
                    x0, x1 = grown

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
        return (f"id={self.id:>3d} q1={self.q1:>3d} q2={self.q2:>3d} f={self.f:>5.2f}GHz segments={len(self.segments)}")

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
    ports: dict[tuple[int, int], tuple[tuple[float, float], tuple[float, float]]] = field(default_factory=dict)
    coupler_regions: dict[tuple[int, int], tuple[float, float, float, float] | None] = field(default_factory=dict)
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
        if (
            config.get("processor") == "grid_25"
            and bool(getattr(params, "grid_compact_mode_enabled", False))
        ):
            compact_side = float(getattr(params, "grid_compact_chip_size_um", 7800.0))
            chip_width = compact_side
            chip_height = compact_side

        freq_ghz = config["freq_ghz"]
        qubits = {
            i: Qubit(id=i, w=params.qubit_width, h=params.qubit_height, f=freq_ghz[i],
                     pad_inset_um=params.qubit_pad_inset_um)
            for i in range(num_qubits)
        }
        cmap = tuple(sorted(set((min(u, v), max(u, v)) for u, v in config["coupling_map"])))
        edge_freq = {(min(u, v), max(u, v)): f for u, v, f in config["edge_freq_ghz"]}
        topology_family = config.get("topology_family")
        region_max_aspect = 8.0 if topology_family == "grid" else _REGION_MAX_ASPECT
        couplers = {
            (q1, q2): Coupler(
                id=idx, q1=q1, q2=q2, f=edge_freq[(q1, q2)],
                epsilon_r=params.substrate_epsilon_r,
                meander_spacing_um=params.meander_spacing_um,
                segment_size_um=params.segment_size_um,
                box_slack_ratio=params.box_slack_ratio,
                region_max_aspect=region_max_aspect,
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
