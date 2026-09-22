import logging
import os
import csv

import matplotlib
matplotlib.use('Agg')  # headless — main.py는 디스플레이 없는 환경(서버)에서도 돌아가야 한다
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb, to_rgba
from shapely.geometry import box
from shapely.ops import unary_union

from core.router import find_crossing_legs
from core.state import ChipState

# 요소가 die 밖으로 나갔을 때 강조할 색 — 포트 강조(주황)와 구분되도록 별도로 예약.
_OOB_COLOR = "#d62728"
# RT가 실제로 그은 배선(성공)을 표시하는 색 — FP~DP 단계의 "논리적" 직선(파랑, #4c72b0)과
# 구분해야 "실제 배선이냐 그냥 포트-포트 직선이냐"를 한눈에 알 수 있다. 교차 구간은 die-밖과
# 같은 경고색(_OOB_COLOR)을 재사용한다 — 둘 다 "잘못된 상태"라는 같은 신호이기 때문.
_ROUTE_COLOR = "#2ca02c"


class Rendering:
    """dir에 state(list[ChipState])의 칩마다 {processor_name}.png를 저장하는 pure-matplotlib 렌더러.

    /home/LabMember/ngchoi/research/qiskit_metal/transmon_draw2.py의 QubitMacro.draw_pure /
    render_pure_python 스타일(큐빗=사각형+라벨, 포트=작은 마커)을 참고했다. 그쪽은 5큐빗
    고정 레이아웃이라 그대로 쓸 수 없다 — 이 파이프라인은 최대 127큐빗(eagle)까지 다뤄야
    해서 라벨 생략/축소 로직이 필요하고, segments/bbox()는 GP 이전(FP 직후)엔 없으므로
    있는 것만 그리는 분기가 필요하다.
    """

    def __init__(self, dir, state):
        os.makedirs(dir, exist_ok=True)
        for chip in state:
            self._render_chip(dir, chip)

    # -----------------------------------------------------------------
    # 칩 1개 렌더링
    # -----------------------------------------------------------------
    def _render_chip(self, dir, chip: ChipState):
        fig, ax = plt.subplots(figsize=self._figsize(chip))

        self._draw_die_boundary(ax, chip)

        ports_by_qubit = self._ports_by_qubit(chip)
        total_segments = sum(len(c.segments) for c in chip.couplers.values())
        has_segments = total_segments > 0
        has_boxes = False
        has_bbox = False

        # RT 이전 단계(0_FP~3_DP)는 모든 coupler.waypoints가 비어 있으므로(core/state.py의
        # Coupler.waypoints 기본값) has_routes=False로 자동 폴백 — 이 렌더러 하나로 5단계
        # 전부 그린다는 기존 구조를 그대로 유지한다.
        n_routed = sum(1 for c in chip.couplers.values() if c.waypoints)
        has_routes = n_routed > 0
        n_cross, crossing_legs = find_crossing_legs(chip) if has_routes else (0, {})

        # FP 박스를 먼저 그려 맨 아래 깔아 둔다 — 세그먼트/bbox()가 그 위에 겹쳐 보여야
        # "박스 안에 세그먼트가 머무는지"를 한눈에 비교할 수 있다.
        for key in sorted(chip.couplers):
            box = chip.coupler_regions.get(key)
            if box is not None:
                has_boxes = True
                self._draw_coupler_box(ax, chip, box)

        for key in sorted(chip.couplers):
            coupler = chip.couplers[key]
            if coupler.waypoints:
                self._draw_coupler_route(ax, coupler, crossing_legs.get(key, set()))
            elif coupler.segments:
                # RT가 시도했지만 경로를 못 찾은 커플러(core/router.py의 route_failures) —
                # 보통의 논리적 직선(파랑)과 구분되는 경고색 점선으로 표시한다.
                self._draw_coupler_failed(ax, chip, coupler)
            else:
                self._draw_coupler_line(ax, chip, coupler)
            if coupler.segments:
                self._draw_segments(ax, coupler)
            bbox = coupler.bbox()
            if bbox is not None:
                has_bbox = True
                self._draw_coupler_bbox(ax, chip, bbox)

        label_mode = self._label_mode(chip.num_qubits)
        for qid in sorted(chip.qubits):
            self._draw_qubit(ax, chip, chip.qubits[qid], ports_by_qubit, label_mode)

        title = (f"{chip.processor_name} | qubits={chip.num_qubits} "
                 f"couplers={len(chip.couplers)} boxes={sum(1 for k in chip.couplers if k in chip.coupler_regions)} "
                 f"segments={total_segments}")
        if has_routes:
            title += f" | routed={n_routed} crossings={n_cross}"
        if chip.fallback_used:
            title += " | fallback"
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
        ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.3)
        ax.margins(0.04)
        ax.set_aspect("equal", adjustable="box")

        # self._draw_legend(ax, has_segments, has_boxes, has_bbox)

        fig.tight_layout()
        fig.savefig(os.path.join(dir, f"{chip.processor_name}.png"), dpi=160)
        plt.close(fig)

    # -----------------------------------------------------------------
    # 개별 요소
    # -----------------------------------------------------------------
    def _draw_die_boundary(self, ax, chip: ChipState):
        ax.add_patch(patches.Rectangle(
            (0.0, 0.0), chip.chip_width, chip.chip_height,
            facecolor="none", edgecolor="black", linewidth=1.3, zorder=0,
        ))

    def _draw_qubit(self, ax, chip: ChipState, qubit, ports_by_qubit: dict, label_mode: str):
        x0, x1 = qubit.x - qubit.w / 2.0, qubit.x + qubit.w / 2.0
        y0, y1 = qubit.y - qubit.h / 2.0, qubit.y + qubit.h / 2.0
        oob = not self._in_die(x0, x1, y0, y1, chip)

        ax.add_patch(patches.Rectangle(
            (x0, y0), qubit.w, qubit.h,
            facecolor="#dbe9f6", alpha=0.65,
            edgecolor=_OOB_COLOR if oob else "#1f4e79",
            linewidth=1.8 if oob else 0.8,
            zorder=7 if oob else 4,
        ))

        # 2026-09-21: 포트가 이웃 방향(qubit_port_toward)으로 바뀌며 "미사용 슬롯"이라는
        # 개념 자체가 없어졌다 — 큐빗은 실제 이웃 수만큼만 포트를 갖고, 전부 쓰인다.
        # 그래서 회색(미사용)/주황(사용) 구분 없이 전부 같은 색으로 그린다.
        port_r = (0.045 if label_mode != "none" else 0.03) * min(qubit.w, qubit.h)
        for px, py in ports_by_qubit.get(qubit.id, []):
            ax.add_patch(patches.Circle(
                (px, py), port_r,
                facecolor="#ff7f0e", edgecolor="none", alpha=0.9,
                zorder=6,
            ))

        if label_mode == "full":
            ax.text(qubit.x, qubit.y, f"{qubit.id}\n{qubit.f:.2f}GHz",
                     ha="center", va="center", fontsize=6, zorder=8)
        elif label_mode == "id":
            ax.text(qubit.x, qubit.y, f"{qubit.id}",
                     ha="center", va="center", fontsize=4.5, zorder=8)

    def _draw_coupler_line(self, ax, chip: ChipState, coupler):
        q1, q2 = chip.qubits[coupler.q1], chip.qubits[coupler.q2]
        pts = chip.ports.get((coupler.q1, coupler.q2))
        if pts is not None:
            (x0, y0), (x1, y1) = pts
        else:
            x0, y0 = q1.x, q1.y
            x1, y1 = q2.x, q2.y

        oob = not (self._point_in_die(x0, y0, chip) and self._point_in_die(x1, y1, chip))
        ax.plot(
            [x0, x1], [y0, y1],
            color=_OOB_COLOR if oob else "#4c72b0",
            linewidth=1.6 if oob else 0.9,
            alpha=0.9 if oob else 0.55,
            zorder=6 if oob else 2,
        )

    # RT가 확정한 실제 배선(core/state.py의 Coupler.waypoints) — 꺾은선 그대로 그린다.
    # crossing_legs에 든 leg 인덱스(waypoints[i]->waypoints[i+1])만 경고색으로 굵게 강조한다
    # (core/router.py의 find_crossing_legs가 커플러 키 -> leg 인덱스 집합으로 반환).
    def _draw_coupler_route(self, ax, coupler, crossing_legs: set):
        wp = coupler.waypoints
        for i in range(len(wp) - 1):
            (x0, y0), (x1, y1) = wp[i], wp[i + 1]
            is_crossing = i in crossing_legs
            ax.plot(
                [x0, x1], [y0, y1],
                color=_OOB_COLOR if is_crossing else _ROUTE_COLOR,
                linewidth=2.6 if is_crossing else 1.2,
                alpha=1.0 if is_crossing else 0.8,
                zorder=9 if is_crossing else 3,
                solid_capstyle="round",
            )

    # RT가 시도했으나 경로를 못 찾은 커플러(coupler.segments는 있지만 waypoints가 빈 채로
    # 남음) — 포트-포트 직선을 경고색 점선으로 그려 "여기 배선이 빠졌다"를 표시한다.
    def _draw_coupler_failed(self, ax, chip: ChipState, coupler):
        q1, q2 = chip.qubits[coupler.q1], chip.qubits[coupler.q2]
        pts = chip.ports.get((coupler.q1, coupler.q2))
        if pts is None:
            x0, y0, x1, y1 = q1.x, q1.y, q2.x, q2.y
        else:
            (x0, y0), (x1, y1) = pts
        ax.plot([x0, x1], [y0, y1], color=_OOB_COLOR, linewidth=1.3,
                linestyle=":", alpha=0.75, zorder=5)

    def _draw_segments(self, ax, coupler):
        color = "#808080"
        half = coupler.segment_size_um / 2.0
        # bbox()가 None이 아닌 칩에서만 die 크기를 알 수 있으므로, die 판정은 커플러 영역
        # 그리기(_draw_coupler_region)에서 이미 한 번 수행한다 — 세그먼트 자체는 bbox 안에
        # 있는 게 보장되므로(Coupler.bbox() 정의) 여기선 개별 die-경계 판정을 반복하지 않는다.
        for seg in coupler.segments:
            ax.add_patch(patches.Rectangle(
                (seg.x - half, seg.y - half), coupler.segment_size_um, coupler.segment_size_um,
                facecolor=color, alpha=0.55, edgecolor=color, linewidth=0.4, zorder=3,
            ))

    # FP가 build_ports 직후 확정한 배치 후보 영역(ChipState.coupler_regions) — GP 이전
    # (output/0_FP/)에도 이미 존재하고, GP 이후로는 GP의 하드 탐색 범위 그 자체가 된다.
    # 옅은 채움 + 점선으로, 그 안에 그려질 세그먼트/bbox()보다 눈에 덜 띄게 깔아 둔다.
    # 테두리는 채움과 다른(더 짙은) 회색+alpha를 쓴다 — Rectangle에 공용 alpha= 하나만
    # 주면 테두리도 채움과 같은 옅은 투명도를 먹어 경계가 거의 안 보였다. facecolor/
    # edgecolor를 각각 RGBA(to_rgba)로 만들어 넘기면 그 둘을 독립적으로 조절할 수 있다.
    def _draw_coupler_box(self, ax, chip: ChipState, box):
        x0, x1, y0, y1 = box
        oob = not self._in_die(x0, x1, y0, y1, chip)
        ax.add_patch(patches.Rectangle(
            (x0, y0), x1 - x0, y1 - y0,
            facecolor=to_rgba(_OOB_COLOR if oob else "#999999", 0.15 if oob else 0.07),
            edgecolor=to_rgba(_OOB_COLOR if oob else "#555555", 0.15 if oob else 0.6),
            linewidth=1.2 if oob else 0.6,
            linestyle=":", zorder=6 if oob else 1,
        ))

    # GP가 실제로 배치한 세그먼트들의 bounding box(Coupler.bbox(), 세그먼트 없으면 None —
    # 그래서 GP 이전엔 안 그려짐). 박스(coupler_regions)를 넘지 않아야 정상이므로(GP는
    # 박스 밖에 안 놓는다) 박스보다 진한 실선 대시로 대비를 줘서 "안에 잘 들어갔는지"
    # 비교하기 쉽게 했다.
    def _draw_coupler_bbox(self, ax, chip: ChipState, bbox):
        x0, x1, y0, y1 = bbox
        oob = not self._in_die(x0, x1, y0, y1, chip)
        ax.add_patch(patches.Rectangle(
            (x0, y0), x1 - x0, y1 - y0,
            facecolor="none",
            edgecolor=_OOB_COLOR if oob else "#555555",
            linewidth=1.2 if oob else 0.7,
            linestyle="--", alpha=0.8 if oob else 0.45,
            zorder=6 if oob else 2,
        ))

    def _draw_legend(self, ax, has_segments: bool, has_boxes: bool, has_bbox: bool):
        handles = [
            patches.Patch(facecolor="none", edgecolor="black", linewidth=1.3, label="die boundary"),
            patches.Patch(facecolor="#dbe9f6", edgecolor="#1f4e79", label="qubit"),
            plt.Line2D([], [], marker="o", color="none", markerfacecolor="#ff7f0e", markersize=6, label="port"),
            plt.Line2D([], [], color="#4c72b0", linewidth=1.5, label="coupler"),
        ]
        if has_segments:
            handles.append(patches.Patch(facecolor="gray", alpha=0.55, edgecolor="gray", label="segment"))
        if has_boxes:
            handles.append(patches.Patch(facecolor="#999999", alpha=0.2, edgecolor="#999999",
                                          linestyle=":", label="coupler box (FP)"))
        if has_bbox:
            handles.append(patches.Patch(facecolor="none", edgecolor="#555555", linestyle="--", label="segment bbox (GP)"))
        handles.append(patches.Patch(facecolor="none", edgecolor=_OOB_COLOR, linewidth=1.5, label="out of die"))
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                   fontsize=6, frameon=False, borderaxespad=0.0)

    # -----------------------------------------------------------------
    # 판단 로직 (근거는 주석 참고)
    # -----------------------------------------------------------------
    def _label_mode(self, num_qubits: int) -> str:
        # 이 저장소 6칩의 num_qubits: grid_25=25, falcon=27, aspen_11=40, xtree_53=53,
        # aspen_m=80, eagle=127. id+freq 두 줄 라벨은 큐빗이 서로 붙어 있을 때(400um 큐빗,
        # min_qubit_spacing_um=400) 쉽게 겹친다 — 30개 이하(grid_25/falcon)는 격자/트리형이라
        # 실측상 라벨 공간이 있어 두 줄 다 표시. 30~60개(aspen_11/xtree_53)는 id 한 줄로
        # 줄여 폰트를 낮추고, 60개 초과(aspen_m/eagle)는 라벨을 아예 생략한다 — eagle 127
        # 큐빗을 A4 비율 캔버스에 다 표시하면 8pt 미만 폰트도 서로 겹쳐 글자 형태가 깨진다.
        if num_qubits <= 30:
            return "full"
        if num_qubits <= 60:
            return "id"
        return "none"

    def _figsize(self, chip: ChipState):
        base = 9.0
        w, h = float(chip.chip_width), float(chip.chip_height)
        if w >= h:
            return (base, max(base * h / w, 3.0))
        return (max(base * w / h, 3.0), base)

    def _coupler_color(self, coupler_id: int):
        # 커플러 개수가 chip마다 달라(eagle=144개) matplotlib 기본 10색 순환은 바로 반복돼
        # 인접 커플러가 같은 색이 되기 쉽다. 황금각(golden angle) hue 스텝은 색 개수와
        # 무관하게 인접 id끼리 색이 잘 갈리는 결정론적(재현 가능한) 팔레트를 준다.
        hue = (coupler_id * 0.6180339887498949) % 1.0
        return hsv_to_rgb((hue, 0.65, 0.85))

    # 큐빗 id -> 그 큐빗이 갖는 모든 포트 좌표(이웃마다 하나씩, core/state.py의
    # qubit_port_toward). 2026-09-21 이전엔 "이 큐빗의 4개 고정 슬롯 중 어느 게 실제로
    # 이웃에 배정됐는지"(used_ports)였는데, 이제 포트 자체가 이웃 수만큼만 존재해 전부
    # "사용됨"이므로 그 구분이 없어졌다.
    def _ports_by_qubit(self, chip: ChipState) -> dict:
        by_qubit: dict[int, list] = {}
        for (q1, q2), (p1, p2) in chip.ports.items():
            by_qubit.setdefault(q1, []).append(p1)
            by_qubit.setdefault(q2, []).append(p2)
        return by_qubit

    def _in_die(self, x0, x1, y0, y1, chip: ChipState, eps: float = 1e-6) -> bool:
        return x0 >= -eps and x1 <= chip.chip_width + eps and y0 >= -eps and y1 <= chip.chip_height + eps

    def _point_in_die(self, x, y, chip: ChipState, eps: float = 1e-6) -> bool:
        return -eps <= x <= chip.chip_width + eps and -eps <= y <= chip.chip_height + eps


# TransmonPocket.connection_pads가 물리적으로 지원하는 슬롯은 이 4개(모서리, loc_W/loc_H
# ∈ {-1,+1})뿐이다 — qiskit-metal 소스 직접 확인함(임의 경계 위치 불가, PlacementInfeasible
# 체크 사유이기도 함, core/floorplan.py의 build_ports 참고). (loc_W, loc_H) -> 슬롯 이름.
_QM_SLOTS: dict[tuple[int, int], str] = {
    (-1, +1): "p_top_left",
    (-1, -1): "p_bottom_left",
    (+1, +1): "p_top_right",
    (+1, -1): "p_bottom_right",
}


# state.ports[key]의 연속 좌표(이웃 방향, qubit_port_toward)를 그 큐빗 중심 기준
# 4사분면 중 가장 가까운 것으로 반올림해 (loc_W, loc_H) 슬롯을 고른다 — 정확한 경계
# 위치가 아니라 "어느 모서리에 가장 가까운가"만 필요하므로 부호만 보면 된다.
def _nearest_qm_slot(qubit: "Qubit", port: tuple[float, float]) -> tuple[int, int]:
    return (1 if port[0] >= qubit.x else -1, 1 if port[1] >= qubit.y else -1)


# state(ChipState 하나)로부터 큐빗별 (loc_W,loc_H) 슬롯 배정을 계산한다. 한 큐빗의 이웃이
# 5개 이상이면(qiskit-metal 물리 한계는 4개) 여러 이웃이 같은 슬롯으로 반올림될 수 있다 —
# 그 경우 정렬 순서상 먼저 온 이웃이 그 슬롯을 차지하고, 나머지는 자리를 못 얻는다.
# 반환: qubit_id -> {(loc_W,loc_H): coupler_key}, (coupler_key, qubit_id) -> slot 이름.
def _assign_qm_slots(state: ChipState):
    slot_of_qubit: dict[int, dict[tuple[int, int], tuple[int, int]]] = {qid: {} for qid in state.qubits}
    slot_name_of: dict[tuple[tuple[int, int], int], str] = {}

    for key in sorted(state.ports):
        port1, port2 = state.ports[key]
        for qid, port in ((key[0], port1), (key[1], port2)):
            qubit = state.qubits[qid]
            slot = _nearest_qm_slot(qubit, port)
            taken_by = slot_of_qubit[qid].get(slot)
            if taken_by is not None:
                logging.warning(
                    "[QM] %s: qubit %d의 슬롯 %s를 커플러 %s가 이미 차지 -- 커플러 %s의 "
                    "이 쪽 핀은 근사 렌더링에서 생략(4-슬롯 근사 한계, 실제 배치는 이웃 "
                    "방향 포트라 겹치지 않음)",
                    state.processor_name, qid, _QM_SLOTS[slot], taken_by, key,
                )
                continue
            slot_of_qubit[qid][slot] = key
            slot_name_of[(key, qid)] = _QM_SLOTS[slot]

    return slot_of_qubit, slot_name_of


# ChipState 하나를 qiskit-metal CAD로 렌더링해 {dir}/{processor_name}.png로 저장한다.
#
# TransmonPocket의 connection_pads는 4개 고정 슬롯만 지원하므로 큐빗 패드 위치는 여전히
# 가까운 슬롯으로 근사한다. 하지만 커플러 경로는 RouteMeander로 재합성하지 않고 RT가 확정한
# Coupler.waypoints를 qiskit-metal qgeometry.path에 직접 넣는다. 따라서 4_RT의 라우팅
# 꺾은선과 이 CAD 그림의 커플러 경로는 같은 데이터를 기준으로 한다.
#
# main.py가 Rendering()과 같은 호출 규약(dir, list[ChipState])을 쓰므로 이 함수도 그렇게
# 받는다.
def rendering_qiskit_metal(dir, states):
    os.environ.setdefault('QISKIT_METAL_HEADLESS', '1')
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

    try:
        from qiskit_metal import QComponent, draw
        from qiskit_metal import designs
        from qiskit_metal.qlibrary.qubits.transmon_pocket import TransmonPocket
    except Exception as e:
        logging.warning("[QM] qiskit-metal import 실패(%s) -- qiskit-metal 렌더링을 건너뜁니다.", e)
        return

    os.makedirs(dir, exist_ok=True)

    class RoutedCouplerPath(QComponent):
        default_options = dict(points=[], trace_width='10um', trace_gap='6um')

        def make(self):
            pts = [
                (self.design.parse_value(f'{x}um'), self.design.parse_value(f'{y}um'))
                for x, y in self.options.points
            ]
            if len(pts) < 2:
                return
            line = draw.LineString(pts)
            trace_width = self.design.parse_value(self.options.trace_width)
            trace_gap = self.design.parse_value(self.options.trace_gap)
            self.options._actual_length = f'{line.length} {self.design.get_units()}'
            self.add_qgeometry('path', {'trace': line}, width=trace_width, fillet=0.0)
            self.add_qgeometry(
                'path', {'cut': line},
                width=trace_width + 2.0 * trace_gap,
                fillet=0.0,
                subtract=True,
            )

    length_rows: list[dict[str, object]] = []

    for state in states:
        slot_of_qubit, slot_name_of = _assign_qm_slots(state)
        design = designs.DesignPlanar()
        rendered_keys: list[tuple[int, int]] = []
        missing_keys: list[tuple[int, int]] = []
        qubit_components = {}

        for qid, qubit in state.qubits.items():
            connection_pads = {
                _QM_SLOTS[slot]: dict(loc_W=slot[0], loc_H=slot[1], pad_width='30um', pad_gap='10um')
                for slot in slot_of_qubit[qid]
            }
            qubit_components[qid] = TransmonPocket(design, f"Q{qid}", options=dict(
                pos_x=f'{qubit.x}um', pos_y=f'{qubit.y}um',
                pocket_width=f'{qubit.w}um', pocket_height=f'{qubit.h}um',
                pad_width='200um', pad_height='80um', pad_gap='30um',
                connection_pads=connection_pads,
            ))

        for key, coupler in state.couplers.items():
            if not coupler.waypoints:
                missing_keys.append(key)
                continue
            points = _qmetal_pin_aligned_waypoints(
                coupler.waypoints,
                qubit_components.get(key[0]), slot_name_of.get((key, key[0])),
                qubit_components.get(key[1]), slot_name_of.get((key, key[1])),
            )
            rendered_keys.append(key)
            RoutedCouplerPath(
                design,
                f"C{key[0]}_{key[1]}",
                options=dict(points=points, trace_width='10um', trace_gap='6um'),
            )

        for key, coupler in state.couplers.items():
            if not coupler.waypoints:
                length_rows.append({
                    "processor": state.processor_name,
                    "coupler": f"{key[0]}-{key[1]}",
                    "status": "missing_waypoints",
                    "target_um": f"{coupler.l:.6f}",
                    "rt_length_um": "",
                    "qmetal_length_um": "",
                    "err_pct": "",
                    "qmetal_visual_err_pct": "",
                })
                continue
            qmetal_length_um = _qmetal_path_length_um(design, f"C{key[0]}_{key[1]}")
            rt_length_um = coupler.route_length
            length_rows.append({
                "processor": state.processor_name,
                "coupler": f"{key[0]}-{key[1]}",
                "status": "rendered",
                "target_um": f"{coupler.l:.6f}",
                "rt_length_um": f"{rt_length_um:.6f}",
                "qmetal_length_um": "" if qmetal_length_um is None else f"{qmetal_length_um:.6f}",
                "err_pct": f"{((rt_length_um - coupler.l) / coupler.l * 100.0):.6f}",
                "qmetal_visual_err_pct": "" if qmetal_length_um is None else f"{((qmetal_length_um - coupler.l) / coupler.l * 100.0):.6f}",
            })

        _log_qmetal_length_summary(state, length_rows)
        if missing_keys:
            logging.warning(
                "[QM] %s: RT waypoints 없는 커플러 %d/%d개 -- PNG에 빨간 점선으로 표시합니다.",
                state.processor_name, len(missing_keys), len(state.couplers),
            )

        try:
            fig, ax = plt.subplots(figsize=(12, 12))
            for _, table in design.qgeometry.tables.items():
                if not table.empty:
                    table.plot(ax=ax, alpha=0.6, edgecolor='blue')
            _draw_qmetal_cpw_overlay(ax, design)
            _draw_missing_qmetal_couplers(ax, state, missing_keys)
            ax.set_aspect('equal')
            ax.set_title(
                f"{state.processor_name} | Qiskit-Metal CAD "
                f"(routed={len(rendered_keys)}, missing={len(missing_keys)}, 4-slot approx.)"
            )
            plt.grid(True, linestyle='--', alpha=0.4)
            plt.savefig(os.path.join(dir, f"{state.processor_name}.png"), dpi=200, bbox_inches='tight')
            plt.close(fig)
        except Exception as e:
            logging.warning("[QM] %s: qiskit-metal qgeometry plot 실패(%s)", state.processor_name, e)
            plt.close('all')

    if length_rows:
        report_path = os.path.join(dir, "qmetal_length_report.csv")
        with open(report_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "processor", "coupler", "status",
                    "target_um", "rt_length_um", "qmetal_length_um",
                    "err_pct", "qmetal_visual_err_pct",
                ],
            )
            writer.writeheader()
            writer.writerows(length_rows)
        logging.info("[QM] qiskit-metal length report 저장: %s", report_path)


def rendering_qiskit_metal_native(dir, states, params=None):
    os.environ.setdefault('QISKIT_METAL_HEADLESS', '1')
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

    try:
        from qiskit_metal import QComponent, draw
        from qiskit_metal import designs
        from qiskit_metal.qlibrary.qubits.transmon_pocket import TransmonPocket
        from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
        from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
    except Exception as e:
        logging.warning("[QM-native] qiskit-metal import 실패(%s) -- native 렌더링을 건너뜁니다.", e)
        return

    os.makedirs(dir, exist_ok=True)

    class FallbackRoutedCouplerPath(QComponent):
        default_options = dict(points=[], trace_width='10um', trace_gap='6um')

        def make(self):
            pts = [
                (self.design.parse_value(f'{x}um'), self.design.parse_value(f'{y}um'))
                for x, y in self.options.points
            ]
            if len(pts) < 2:
                return
            line = draw.LineString(pts)
            trace_width = self.design.parse_value(self.options.trace_width)
            trace_gap = self.design.parse_value(self.options.trace_gap)
            self.options._actual_length = f'{line.length} {self.design.get_units()}'
            self.add_qgeometry('path', {'trace': line}, width=trace_width, fillet=0.0)
            self.add_qgeometry(
                'path', {'cut': line},
                width=trace_width + 2.0 * trace_gap,
                fillet=0.0,
                subtract=True,
            )

    report_rows: list[dict[str, object]] = []
    meander_spacings_um = _qmetal_native_meander_spacings(params)
    lead_straight_um = _qmetal_native_lead_straight_um(params)
    length_tolerance_pct = _qmetal_native_length_tolerance_pct(params)
    cut_half_width_mm = (
        float(getattr(params, "cpw_trace_width_um", 10.0)) / 2.0
        + float(getattr(params, "cpw_trace_gap_um", 6.0))
    ) / 1000.0
    min_route_spacing_mm = (
        float(getattr(params, "qmetal_native_min_route_spacing_um", 0.0)) / 1000.0
    )

    for state in states:
        slot_of_qubit, slot_name_of = _assign_qm_slots(state)
        design = designs.DesignPlanar()
        accepted_native_routes = []
        qubit_components = {}
        qubit_bodies_mm = {
            qid: (
                (qubit.x - qubit.w / 2.0) / 1000.0,
                (qubit.y - qubit.h / 2.0) / 1000.0,
                (qubit.x + qubit.w / 2.0) / 1000.0,
                (qubit.y + qubit.h / 2.0) / 1000.0,
            )
            for qid, qubit in state.qubits.items()
        }

        for qid, qubit in state.qubits.items():
            connection_pads = {
                _QM_SLOTS[slot]: dict(loc_W=slot[0], loc_H=slot[1], pad_width='30um', pad_gap='10um')
                for slot in slot_of_qubit[qid]
            }
            qubit_components[qid] = TransmonPocket(design, f"Q{qid}", options=dict(
                pos_x=f'{qubit.x}um', pos_y=f'{qubit.y}um',
                pocket_width=f'{qubit.w}um', pocket_height=f'{qubit.h}um',
                pad_width='200um', pad_height='80um', pad_gap='30um',
                connection_pads=connection_pads,
            ))

        rendered = 0
        skipped = 0
        for key, coupler in state.couplers.items():
            pin1 = slot_name_of.get((key, key[0]))
            pin2 = slot_name_of.get((key, key[1]))
            row = {
                "processor": state.processor_name,
                "coupler": f"{key[0]}-{key[1]}",
                "target_um": f"{coupler.l:.6f}",
                "actual_um": "",
                "err_pct": "",
                "meander_spacing_um": "",
                "native_drc": "",
                "status": "",
            }
            if pin1 is None or pin2 is None:
                skipped += 1
                row["status"] = "missing_pin_slot"
                report_rows.append(row)
                continue
            try:
                route, spacing_um = _add_native_route_meander(
                    RouteMeander, design, key, pin1, pin2, coupler.l,
                    meander_spacings_um, lead_straight_um, length_tolerance_pct,
                    accepted_native_routes, qubit_bodies_mm,
                    cut_half_width_mm, min_route_spacing_mm,
                )
                actual_um = _qmetal_component_actual_length_um(design, route)
                line_geom, cut_geom = _qmetal_route_geometries(
                    design, route.name, cut_half_width_mm
                )
                if line_geom is not None and cut_geom is not None:
                    accepted_native_routes.append((key, line_geom, cut_geom))
                rendered += 1
                row["status"] = "rendered_native"
                row["native_drc"] = "clean"
                row["meander_spacing_um"] = f"{spacing_um:.6f}"
                if actual_um is not None:
                    row["actual_um"] = f"{actual_um:.6f}"
                    row["err_pct"] = f"{((actual_um - coupler.l) / coupler.l * 100.0):.6f}"
            except Exception as e:
                if coupler.waypoints:
                    fallback_name = f"F{key[0]}_{key[1]}"
                    points = _qmetal_pin_aligned_waypoints(
                        coupler.waypoints,
                        qubit_components.get(key[0]), pin1,
                        qubit_components.get(key[1]), pin2,
                    )
                    FallbackRoutedCouplerPath(
                        design, fallback_name,
                        options=dict(points=points, trace_width='10um', trace_gap='6um'),
                    )
                    line_geom, cut_geom = _qmetal_route_geometries(
                        design, fallback_name, cut_half_width_mm
                    )
                    if line_geom is not None and cut_geom is not None:
                        accepted_native_routes.append((key, line_geom, cut_geom))
                    actual_um = _qmetal_path_length_um(design, fallback_name)
                    rendered += 1
                    row["status"] = "fallback_pnr_path"
                    row["native_drc"] = f"RouteMeander rejected:{type(e).__name__}"
                    if actual_um is not None:
                        row["actual_um"] = f"{actual_um:.6f}"
                        row["err_pct"] = f"{((actual_um - coupler.l) / coupler.l * 100.0):.6f}"
                    logging.warning(
                        "[QM-native] %s: coupler %s RouteMeander 충돌, PnR qgeometry path로 대체(%s)",
                        state.processor_name, key, e,
                    )
                else:
                    skipped += 1
                    row["status"] = f"route_error:{type(e).__name__}"
                    logging.warning("[QM-native] %s: coupler %s RouteMeander 실패(%s)", state.processor_name, key, e)
            report_rows.append(row)

        try:
            fig, ax = plt.subplots(figsize=(12, 12))
            renderer = QMplRenderer(None, design, logging.getLogger())
            renderer.render(ax)
            _fit_axis_to_qmetal_geometry(ax, design)
            ax.set_aspect('equal')
            ax.set_title(
                f"{state.processor_name} | Qiskit-Metal native RouteMeander "
                f"(rendered={rendered}, skipped={skipped})"
            )
            ax.grid(True, linestyle='-', linewidth=0.5, alpha=0.18)
            fig.savefig(os.path.join(dir, f"{state.processor_name}.png"), dpi=200, bbox_inches='tight')
            plt.close(fig)
        except Exception as e:
            logging.warning("[QM-native] %s: QMplRenderer 실패(%s)", state.processor_name, e)
            plt.close('all')

        chip_rows = [
            r for r in report_rows
            if r["processor"] == state.processor_name
            and r["status"] in {"rendered_native", "fallback_pnr_path"}
            and r["err_pct"]
        ]
        if chip_rows:
            errs = [abs(float(r["err_pct"])) for r in chip_rows]
            logging.info(
                "[QM-native] %s: RouteMeander rendered=%d/%d max_abs_err=%.3f%% within5=%.1f%% spacings=%s",
                state.processor_name, rendered, len(state.couplers), max(errs),
                sum(1 for e in errs if e <= 5.0) / len(errs) * 100.0,
                _qmetal_native_spacing_histogram(chip_rows),
            )
        else:
            logging.info("[QM-native] %s: RouteMeander rendered=0/%d", state.processor_name, len(state.couplers))

    if report_rows:
        report_path = os.path.join(dir, "qmetal_native_route_report.csv")
        with open(report_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "processor", "coupler", "status", "target_um", "actual_um",
                    "err_pct", "meander_spacing_um", "native_drc",
                ],
            )
            writer.writeheader()
            writer.writerows(report_rows)
        logging.info("[QM-native] qiskit-metal native route report 저장: %s", report_path)


def _qmetal_native_meander_spacings(params) -> list[float]:
    values = getattr(params, "qmetal_native_meander_spacings_um", None)
    if values is None:
        values = [200.0]
    if isinstance(values, (int, float)):
        values = [float(values)]
    out = []
    for value in values:
        try:
            spacing = float(value)
        except (TypeError, ValueError):
            continue
        if spacing > 0.0 and spacing not in out:
            out.append(spacing)
    return sorted(out or [200.0], reverse=True)


def _qmetal_native_lead_straight_um(params) -> float:
    try:
        return max(0.0, float(getattr(params, "qmetal_native_lead_straight_um", 80.0)))
    except (TypeError, ValueError):
        return 80.0


def _qmetal_native_length_tolerance_pct(params) -> float:
    try:
        return max(0.0, float(getattr(params, "qmetal_native_length_tolerance_pct", 5.0)))
    except (TypeError, ValueError):
        return 5.0


def _add_native_route_meander(
    route_meander_cls, design, key: tuple[int, int], pin1: str, pin2: str,
    target_length_um: float, spacings_um: list[float], lead_straight_um: float,
    length_tolerance_pct: float, accepted_routes, qubit_bodies_mm: dict,
    cut_half_width_mm: float, min_route_spacing_mm: float,
):
    last_error: Exception | None = None
    best = None
    for spacing_um in spacings_um:
        component_name = f"N{key[0]}_{key[1]}_S{int(round(spacing_um))}"
        try:
            route = route_meander_cls(design, component_name, options=dict(
                pin_inputs=dict(
                    start_pin=dict(component=f"Q{key[0]}", pin=pin1),
                    end_pin=dict(component=f"Q{key[1]}", pin=pin2),
                ),
                trace_width='10um',
                trace_gap='6um',
                fillet='30um',
                total_length=f'{target_length_um}um',
                meander=dict(spacing=f'{spacing_um}um', asymmetry='0um'),
                lead=dict(
                    start_straight=f'{lead_straight_um}um',
                    end_straight=f'{lead_straight_um}um',
                ),
            ))
            actual_um = _qmetal_component_actual_length_um(design, route)
            if actual_um is None or target_length_um <= 0.0:
                actual_ok = True
                abs_err_pct = 0.0
            else:
                abs_err_pct = abs((actual_um - target_length_um) / target_length_um * 100.0)
                actual_ok = abs_err_pct <= length_tolerance_pct
            line_geom, cut_geom = _qmetal_route_geometries(
                design, component_name, cut_half_width_mm
            )
            conflict = _native_route_conflict(
                key, line_geom, cut_geom, accepted_routes, qubit_bodies_mm,
                min_route_spacing_mm,
            )
            if conflict is None:
                if best is None or abs_err_pct < best[0]:
                    best = (abs_err_pct, spacing_um)
                if actual_ok:
                    return route, spacing_um
            _delete_qmetal_component(design, component_name)
            if conflict is not None:
                last_error = RuntimeError(f"spacing {spacing_um}um native DRC conflict: {conflict}")
                continue
            last_error = RuntimeError(
                f"spacing {spacing_um}um length error {abs_err_pct:.3f}% "
                f"> {length_tolerance_pct:.3f}%"
            )
        except Exception as e:
            _delete_qmetal_component(design, component_name)
            last_error = e
    if best is not None:
        _abs_err_pct, spacing_um = best
        component_name = f"N{key[0]}_{key[1]}_S{int(round(spacing_um))}"
        route = route_meander_cls(design, component_name, options=dict(
            pin_inputs=dict(
                start_pin=dict(component=f"Q{key[0]}", pin=pin1),
                end_pin=dict(component=f"Q{key[1]}", pin=pin2),
            ),
            trace_width='10um',
            trace_gap='6um',
            fillet='30um',
            total_length=f'{target_length_um}um',
            meander=dict(spacing=f'{spacing_um}um', asymmetry='0um'),
            lead=dict(
                start_straight=f'{lead_straight_um}um',
                end_straight=f'{lead_straight_um}um',
            ),
        ))
        return route, spacing_um
    if last_error is not None:
        raise last_error
    raise RuntimeError("no RouteMeander spacing candidates")


def _qmetal_route_geometries(design, component_name: str, cut_half_width_mm: float):
    try:
        component = design.components[component_name]
    except Exception:
        return None, None
    table = design.qgeometry.tables.get("path")
    if table is None or table.empty or "component" not in table:
        return None, None
    rows = table[table["component"] == component.id]
    if rows.empty or "geometry" not in rows:
        return None, None
    trace_rows = rows[rows.index == "trace"] if "trace" in rows.index else rows
    line_geom = trace_rows.iloc[0]["geometry"]
    if line_geom is None or line_geom.is_empty:
        return None, None
    cut_geom = line_geom.buffer(cut_half_width_mm, cap_style=2, join_style=2)
    return line_geom, cut_geom


def _native_route_conflict(
    key: tuple[int, int], line_geom, cut_geom, accepted_routes,
    qubit_bodies_mm: dict, min_route_spacing_mm: float,
) -> str | None:
    if line_geom is None or cut_geom is None:
        return "missing_path_geometry"
    for other_key, other_line, other_cut in accepted_routes:
        shared = set(key) & set(other_key)
        test_cut = cut_geom
        test_other_cut = other_cut
        if shared:
            shared_geoms = []
            for qid in shared:
                body = qubit_bodies_mm.get(qid)
                if body is None:
                    continue
                x0, y0, x1, y1 = body
                shared_geoms.append(box(x0, y0, x1, y1))
            if shared_geoms:
                shared_union = unary_union(shared_geoms)
                test_cut = test_cut.difference(shared_union)
                test_other_cut = test_other_cut.difference(shared_union)
        if test_cut.intersection(test_other_cut).area > 1e-12:
            return f"overlap:{other_key[0]}-{other_key[1]}"
        if min_route_spacing_mm > 0.0 and test_cut.distance(test_other_cut) + 1e-12 < min_route_spacing_mm:
            return f"spacing:{other_key[0]}-{other_key[1]}"
        if not shared and line_geom.crosses(other_line):
            return f"centerline_cross:{other_key[0]}-{other_key[1]}"
    return None


def _delete_qmetal_component(design, component_name: str) -> None:
    try:
        design.delete_component(component_name)
    except Exception:
        pass


def _delete_qmetal_component(design, component_name: str) -> None:
    try:
        design.delete_component(component_name)
    except Exception:
        pass


def _qmetal_native_spacing_histogram(rows: list[dict[str, object]]) -> dict[str, int]:
    hist: dict[str, int] = {}
    for row in rows:
        spacing = row.get("meander_spacing_um")
        if not spacing:
            continue
        label = f"{float(spacing):.0f}um"
        hist[label] = hist.get(label, 0) + 1
    return hist


def _qmetal_component_actual_length_um(design, component) -> float | None:
    actual = getattr(component.options, "_actual_length", None)
    if not actual:
        return None
    try:
        return float(design.parse_value(actual)) * 1000.0
    except Exception:
        return None


def _fit_axis_to_qmetal_geometry(ax, design, margin_ratio: float = 0.04) -> None:
    bounds = []
    for table in design.qgeometry.tables.values():
        if table.empty or "geometry" not in table:
            continue
        for geom in table["geometry"]:
            if geom is not None and not geom.is_empty:
                bounds.append(geom.bounds)
    if not bounds:
        return
    min_x = min(b[0] for b in bounds)
    min_y = min(b[1] for b in bounds)
    max_x = max(b[2] for b in bounds)
    max_y = max(b[3] for b in bounds)
    span_x = max(max_x - min_x, 1e-3)
    span_y = max(max_y - min_y, 1e-3)
    pad = max(span_x, span_y) * margin_ratio
    ax.set_xlim(min_x - pad, max_x + pad)
    ax.set_ylim(min_y - pad, max_y + pad)


def _qmetal_path_length_um(design, component_name: str) -> float | None:
    try:
        component = design.components[component_name]
    except Exception:
        return None
    table = design.qgeometry.tables.get("path")
    if table is None or table.empty or "component" not in table:
        return None
    rows = table[table["component"] == component.id]
    if rows.empty or "geometry" not in rows:
        return None
    trace_rows = rows[rows.index == "trace"] if "trace" in rows.index else rows
    geometry = trace_rows.iloc[0]["geometry"]
    return float(geometry.length) * 1000.0


def _qmetal_pin_aligned_waypoints(
    waypoints: list[tuple[float, float]],
    start_component,
    start_pin: str | None,
    end_component,
    end_pin: str | None,
    lead_um: float = 20.0,
) -> list[tuple[float, float]]:
    if len(waypoints) < 2:
        return waypoints

    start = _qmetal_pin_lead_points(start_component, start_pin, lead_um)
    end = _qmetal_pin_lead_points(end_component, end_pin, lead_um)
    if start is None and end is None:
        return waypoints

    interior = list(waypoints[1:-1])
    out: list[tuple[float, float]] = []
    if start is None:
        out.append(waypoints[0])
    else:
        pin_pt, lead_pt = start
        out.extend([pin_pt, lead_pt])
    out.extend(interior)
    if end is None:
        out.append(waypoints[-1])
    else:
        pin_pt, lead_pt = end
        out.extend([lead_pt, pin_pt])
    return _dedupe_points(out)


def _qmetal_pin_lead_points(component, pin_name: str | None, lead_um: float) -> tuple[tuple[float, float], tuple[float, float]] | None:
    if component is None or pin_name is None:
        return None
    pin = component.pins.get(pin_name)
    if pin is None:
        return None
    middle = pin["middle"]
    normal = pin["normal"]
    pin_pt = (float(middle[0]) * 1000.0, float(middle[1]) * 1000.0)
    lead_pt = (
        pin_pt[0] + float(normal[0]) * lead_um,
        pin_pt[1] + float(normal[1]) * lead_um,
    )
    return pin_pt, lead_pt


def _dedupe_points(points: list[tuple[float, float]], eps: float = 1e-6) -> list[tuple[float, float]]:
    if not points:
        return []
    out = [points[0]]
    for point in points[1:]:
        if (abs(point[0] - out[-1][0]) > eps) or (abs(point[1] - out[-1][1]) > eps):
            out.append(point)
    return out


def _log_qmetal_length_summary(state: ChipState, rows: list[dict[str, object]]) -> None:
    chip_rows = [r for r in rows if r["processor"] == state.processor_name and r["status"] == "rendered"]
    if not chip_rows:
        logging.info("[QM] %s: qiskit-metal rendered couplers=0", state.processor_name)
        return
    errs = [abs(float(r["err_pct"])) for r in chip_rows]
    within_5 = sum(1 for e in errs if e <= 5.0)
    logging.info(
        "[QM] %s: rendered=%d/%d RT_length max_abs_err=%.3f%% within5=%.1f%%",
        state.processor_name,
        len(chip_rows),
        len(state.couplers),
        max(errs),
        within_5 / len(chip_rows) * 100.0,
    )


def _draw_qmetal_cpw_overlay(ax, design) -> None:
    table = design.qgeometry.tables.get("path")
    if table is None or table.empty or "geometry" not in table:
        return
    for _, row in table.iterrows():
        geom = row["geometry"]
        if geom.is_empty:
            continue
        subtract = bool(row.get("subtract", False))
        color = "#d8ecff" if subtract else "#1f77b4"
        alpha = 0.22 if subtract else 0.82
        zorder = 1 if subtract else 5
        _plot_linestring(ax, geom, color=color, alpha=alpha, linewidth=2.4 if subtract else 1.8, zorder=zorder)


def _plot_linestring(ax, geom, color: str, alpha: float, linewidth: float, zorder: int) -> None:
    if geom.geom_type == "LineString":
        xs, ys = geom.xy
        ax.plot(xs, ys, color=color, alpha=alpha, linewidth=linewidth, solid_capstyle="round", zorder=zorder)
        return
    if geom.geom_type == "MultiLineString":
        for line in geom.geoms:
            _plot_linestring(ax, line, color, alpha, linewidth, zorder)


def _draw_missing_qmetal_couplers(ax, state: ChipState, missing_keys: list[tuple[int, int]]) -> None:
    for key in missing_keys:
        pts = state.ports.get(key)
        if pts is None:
            q1, q2 = state.qubits[key[0]], state.qubits[key[1]]
            pts = ((q1.x, q1.y), (q2.x, q2.y))
        (x0, y0), (x1, y1) = pts
        ax.plot(
            [x0, x1], [y0, y1],
            color=_OOB_COLOR,
            linewidth=1.4,
            linestyle=":",
            alpha=0.7,
            zorder=6,
        )
