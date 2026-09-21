import logging
import os

import matplotlib
matplotlib.use('Agg')  # headless — main.py는 디스플레이 없는 환경(서버)에서도 돌아가야 한다
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb, to_rgba

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
        # ax.set_xlabel("x (µm)")
        # ax.set_ylabel("y (µm)")
        # ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.3)
        ax.margins(0.04)
        ax.set_aspect("equal", adjustable="box")
        ax.set_axis_off()

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


# ChipState 하나를 qiskit-metal CAD로 근사 렌더링해 {dir}/{processor_name}.png로 저장한다.
#
# 근사인 이유: 이 파이프라인의 실제 포트(core/state.py의 qubit_port_toward)는 큐빗 경계
# 위의 연속한 점(이웃 방향으로의 ray-intersection)이다. qiskit-metal의 TransmonPocket은
# connection_pads로 4개의 고정 모서리 슬롯(loc_W,loc_H ∈ {-1,+1})만 지원하고 임의 경계
# 위치를 표현할 수 없다(소스 직접 확인, TransmonPocket6도 6슬롯이 전부 고정 위치라 근본
# 한계는 같음) — 커스텀 QComponent를 새로 만들면 임의 위치를 표현할 수 있지만 그건 별도
# 작업으로 남겨둔다. 그래서 여기서는 각 포트를 그 큐빗 중심 기준 가장 가까운 4-슬롯 중
# 하나로 반올림해서만 그린다 — **이 그림의 핀 위치·배선 길이는 근사이지 실제 배치가 아니다.**
# 큐빗 하나에 이웃이 5개 이상이면 두 이웃이 같은 슬롯으로 반올림돼 충돌할 수 있다 —
# 그 경우 예외를 던지지 않고 경고 로그만 남기고 그 커플러의 해당 쪽 핀(과 그 커플러 전체
# RouteMeander)을 생략한다. 실제 배치(포트/DRC/라우팅)에는 전혀 영향 없음 — 이 함수는
# 시각적 참고용 CAD 그림만 만든다.
#
# main.py가 Rendering()과 같은 호출 규약(dir, list[ChipState])을 쓰므로 이 함수도 그렇게
# 받는다.
def rendering_qiskit_metal(dir, states):
    os.environ.setdefault('QISKIT_METAL_HEADLESS', '1')
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

    from qiskit_metal import designs
    from qiskit_metal.qlibrary.qubits.transmon_pocket import TransmonPocket
    from qiskit_metal.qlibrary.tlines.meandered import RouteMeander

    os.makedirs(dir, exist_ok=True)

    for state in states:
        slot_of_qubit, slot_name_of = _assign_qm_slots(state)
        design = designs.DesignPlanar()

        for qid, qubit in state.qubits.items():
            connection_pads = {
                _QM_SLOTS[slot]: dict(loc_W=slot[0], loc_H=slot[1], pad_width='30um', pad_gap='10um')
                for slot in slot_of_qubit[qid]
            }
            TransmonPocket(design, f"Q{qid}", options=dict(
                pos_x=f'{qubit.x}um', pos_y=f'{qubit.y}um',
                pocket_width=f'{qubit.w}um', pocket_height=f'{qubit.h}um',
                pad_width='200um', pad_height='80um', pad_gap='30um',
                connection_pads=connection_pads,
            ))

        for key, coupler in state.couplers.items():
            pin1 = slot_name_of.get((key, key[0]))
            pin2 = slot_name_of.get((key, key[1]))
            if pin1 is None or pin2 is None:
                logging.warning(
                    "[QM] %s: coupler %s 근사 렌더링 생략(슬롯 충돌로 한쪽 핀이 없음)",
                    state.processor_name, key,
                )
                continue
            try:
                RouteMeander(design, f"C{key[0]}_{key[1]}", options=dict(
                    pin_inputs=dict(
                        start_pin=dict(component=f"Q{key[0]}", pin=pin1),
                        end_pin=dict(component=f"Q{key[1]}", pin=pin2),
                    ),
                    trace_width='10um',
                    trace_gap='6um',
                    fillet='30um',
                    total_length=f'{coupler.l}um',
                    meander=dict(spacing='120um', asymmetry='0um'),
                    lead=dict(start_straight='150um', end_straight='150um'),
                ))
            except Exception as e:
                logging.warning(
                    "[QM] %s: coupler %s RouteMeander 실패(%s) -- 이 커플러만 생략",
                    state.processor_name, key, e,
                )

        fig, ax = plt.subplots(figsize=(12, 12))
        for _, table in design.qgeometry.tables.items():
            if not table.empty:
                table.plot(ax=ax, alpha=0.6, edgecolor='blue')
        ax.set_aspect('equal')
        ax.set_title(f"{state.processor_name} | Qiskit-Metal CAD (4-slot approx.)")
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.savefig(os.path.join(dir, f"{state.processor_name}.png"), dpi=200, bbox_inches='tight')
        plt.close(fig)
