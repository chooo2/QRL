import os

import matplotlib
matplotlib.use('Agg')  # headless — main.py는 디스플레이 없는 환경(서버)에서도 돌아가야 한다
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb

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

        used_ports = self._used_ports(chip)
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
            self._draw_qubit(ax, chip, chip.qubits[qid], used_ports, label_mode)

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

    def _draw_qubit(self, ax, chip: ChipState, qubit, used_ports: set, label_mode: str):
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

        if label_mode != "none":
            port_r = 0.045 * min(qubit.w, qubit.h)
            for name, (px, py) in qubit.ports.items():
                is_used = (qubit.id, name) in used_ports
                ax.add_patch(patches.Circle(
                    (px, py), port_r,
                    facecolor="#ff7f0e" if is_used else "#888888",
                    edgecolor="none", alpha=0.9 if is_used else 0.55,
                    zorder=6 if is_used else 5,
                ))
        else:
            # 라벨 없는 초밀집 모드(eagle 등)에서도 포트는 남긴다 — 마커만 조금 더 작게.
            port_r = 0.03 * min(qubit.w, qubit.h)
            for name, (px, py) in qubit.ports.items():
                is_used = (qubit.id, name) in used_ports
                ax.add_patch(patches.Circle(
                    (px, py), port_r,
                    facecolor="#ff7f0e" if is_used else "#888888",
                    edgecolor="none", alpha=0.9 if is_used else 0.45,
                    zorder=6 if is_used else 5,
                ))

        if label_mode == "full":
            ax.text(qubit.x, qubit.y, f"{qubit.id}\n{qubit.f:.2f}GHz",
                     ha="center", va="center", fontsize=6, zorder=8)
        elif label_mode == "id":
            ax.text(qubit.x, qubit.y, f"{qubit.id}",
                     ha="center", va="center", fontsize=4.5, zorder=8)

    def _draw_coupler_line(self, ax, chip: ChipState, coupler):
        q1, q2 = chip.qubits[coupler.q1], chip.qubits[coupler.q2]
        assignment = chip.port_assignment.get((coupler.q1, coupler.q2))
        if assignment is not None:
            p1_name, p2_name = assignment
            x0, y0 = q1.ports[p1_name]
            x1, y1 = q2.ports[p2_name]
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
        assignment = chip.port_assignment.get((coupler.q1, coupler.q2))
        if assignment is None:
            x0, y0, x1, y1 = q1.x, q1.y, q2.x, q2.y
        else:
            p1_name, p2_name = assignment
            (x0, y0), (x1, y1) = q1.ports[p1_name], q2.ports[p2_name]
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

    # FP가 assign_ports 직후 확정한 배치 후보 영역(ChipState.coupler_regions) — GP 이전
    # (output/0_FP/)에도 이미 존재하고, GP 이후로는 GP의 하드 탐색 범위 그 자체가 된다.
    # 옅은 채움 + 점선으로, 그 안에 그려질 세그먼트/bbox()보다 눈에 덜 띄게 깔아 둔다.
    def _draw_coupler_box(self, ax, chip: ChipState, box):
        x0, x1, y0, y1 = box
        oob = not self._in_die(x0, x1, y0, y1, chip)
        ax.add_patch(patches.Rectangle(
            (x0, y0), x1 - x0, y1 - y0,
            facecolor=_OOB_COLOR if oob else "#999999",
            alpha=0.15 if oob else 0.07,
            edgecolor=_OOB_COLOR if oob else "#999999",
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
            plt.Line2D([], [], marker="o", color="none", markerfacecolor="#888888", markersize=6, label="port (unused)"),
            plt.Line2D([], [], marker="o", color="none", markerfacecolor="#ff7f0e", markersize=6, label="port (assigned)"),
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

    def _used_ports(self, chip: ChipState) -> set:
        used = set()
        for (q1, q2), (p1, p2) in chip.port_assignment.items():
            used.add((q1, p1))
            used.add((q2, p2))
        return used

    def _in_die(self, x0, x1, y0, y1, chip: ChipState, eps: float = 1e-6) -> bool:
        return x0 >= -eps and x1 <= chip.chip_width + eps and y0 >= -eps and y1 <= chip.chip_height + eps

    def _point_in_die(self, x, y, chip: ChipState, eps: float = 1e-6) -> bool:
        return -eps <= x <= chip.chip_width + eps and -eps <= y <= chip.chip_height + eps


def rendering_qiskit_metal(dir, state):
    pass

    

""" Qiskit-metal rendering reference

import os
import math
import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple, List

# Qt/GUI 헤드리스 오프스크린 환경 설정
os.environ['QISKIT_METAL_HEADLESS'] = '1'
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

from qiskit_metal import designs
from qiskit_metal.qlibrary.qubits.transmon_pocket import TransmonPocket
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander


# ------------------------------------------------------------------------------
# 1. 물리 계산 및 경로(Waypoints) 산출
# ------------------------------------------------------------------------------
def calc_half_wave_length(
    f_ghz: float,           # 목표 공진 주파수 (GHz)
    epsilon_r: float = 11.9 # 기판 유전율 (기본값: 사파이어/실리콘 기판)
) -> float:
#    공진기 주파수(GHz) 기준 반파장 물리 길이(um) 계산
    c = 2.99792458e14  # 빛의 속도 (um/s)
    eps_eff = (epsilon_r + 1.0) / 2.0  # 마이크로스트립 실효 유전율 근사
    return round(c / (2.0 * (f_ghz * 1e9) * math.sqrt(eps_eff)), 2)  # λ/2 = c / (2f·√ε_eff)


def generate_exact_meander_path(
    p_start: Tuple[float, float],   # 경로 시작 좌표 (x, y)
    p_end: Tuple[float, float],     # 경로 끝 좌표 (x, y)
    target_length: float,           # 맞춰야 할 목표 전체 길이 (um)
    orientation: str = 'horizontal',# 지그재그 진행 축: 'horizontal'(좌우 진행/상하 꺾임) or 'vertical'(상하 진행/좌우 꺾임)
    meander_height: float = 250.0   # 한 번 꺾일 때 튀어나오는 폭(진폭) (um)
) -> List[Tuple[float, float]]:
    # 목표 길이에 맞춘 Meander 좌표계 산출
    x1, y1 = p_start
    x2, y2 = p_end

    direct_dist = math.hypot(x2 - x1, y2 - y1)  # 시작-끝 직선 거리
    extra_length = max(0.0, target_length - direct_dist)  # 지그재그로 채워야 할 남는 길이
    num_turns = max(2, int(extra_length / (2.0 * meander_height)))  # 꺾이는 횟수 (최소 2회)

    waypoints = [p_start]
    lead_in = 150.0  # 시작/끝 지점의 직선 리드 구간 길이

    if orientation == 'horizontal':
        # 좌우로 이동하며 위/아래로 번갈아 꺾는 지그재그 경로
        start_x, end_x = x1 + lead_in, x2 - lead_in
        step_x = (end_x - start_x) / (num_turns + 1)
        waypoints.append((start_x, y1))

        curr_x, sign = start_x, 1.0
        for _ in range(num_turns):
            curr_x += step_x
            waypoints.append((curr_x, y1 + sign * meander_height))  # sign이 매번 반전되며 위/아래로 꺾임
            sign *= -1.0

        waypoints.append((end_x, y1))
    else:
        # 상하로 이동하며 좌/우로 번갈아 꺾는 지그재그 경로 (horizontal과 축만 반대)
        start_y, end_y = y1 + lead_in, y2 - lead_in
        step_y = (end_y - start_y) / (num_turns + 1)
        waypoints.append((x1, start_y))

        curr_y, sign = start_y, 1.0
        for _ in range(num_turns):
            curr_y += step_y
            waypoints.append((x1 + sign * meander_height, curr_y))
            sign *= -1.0

        waypoints.append((x1, end_y))

    waypoints.append(p_end)
    return waypoints


# ------------------------------------------------------------------------------
# 2. 데이터 구조
# ------------------------------------------------------------------------------
@dataclass
class QubitMacro:
    name: str            # 큐비트 식별 이름 (예: "Q0")
    x: float             # 중심 x좌표 (um)
    y: float             # 중심 y좌표 (um)
    freq_ghz: float      # 큐비트 자체 주파수 (GHz, 라벨 표시용)
    width: float = 400.0  # Pocket 가로 폭 (um)
    height: float = 400.0 # Pocket 세로 높이 (um)

    @property
    def ports(self) -> Dict[str, Tuple[float, float]]:
        # 큐비트 핀 포트 위치 정보
        hw, off_y = self.width / 2.0, 180.0  # hw: 좌우 절반 폭, off_y: 상/하 포트 y축 오프셋
        return {
            # 큐비트 좌/우 양쪽 모서리에 상단·하단 커플링 핀 4개를 배치
            "p_top_left":     (self.x - hw, self.y + off_y),
            "p_bottom_left":  (self.x - hw, self.y - off_y),
            "p_top_right":    (self.x + hw, self.y + off_y),
            "p_bottom_right": (self.x + hw, self.y - off_y),
        }

    def draw_pure(self, ax):
        # Pure Python 미리보기용 큐비트 시각화
        hw, hh = self.width / 2.0, self.height / 2.0  # Pocket 절반 폭/높이
        pad_w, pad_h, pad_gap = 200.0, 80.0, 30.0     # 커패시터 Pad 폭/높이/상하 간격

        # Pocket(바깥 사각형 공동)과 상/하단 커패시터 Pad
        ax.add_patch(patches.Rectangle((self.x - hw, self.y - hh), self.width, self.height, linewidth=1.2, edgecolor='#1f77b4', facecolor='#a6cee3', alpha=0.5, zorder=1))
        top_y, bot_y = self.y + (pad_h + pad_gap)/2.0, self.y - (pad_h + pad_gap)/2.0
        ax.add_patch(patches.Rectangle((self.x - pad_w/2, top_y - pad_h/2), pad_w, pad_h, edgecolor='#1f77b4', facecolor='#1f77b4', alpha=0.8, zorder=3))
        ax.add_patch(patches.Rectangle((self.x - pad_w/2, bot_y - pad_h/2), pad_w, pad_h, edgecolor='#1f77b4', facecolor='#1f77b4', alpha=0.8, zorder=3))
        ax.plot([self.x, self.x], [top_y - pad_h/2, bot_y + pad_h/2], color='#33a02c', linewidth=2.0, zorder=4)  # 두 Pad를 잇는 조셉슨 접합(junction) 선

        # 좌우 4개 커플링 핀을 작은 사각형으로 표시 (바깥쪽을 향하도록 방향 조정)
        for px, py in self.ports.values():
            dx = 30 if px > self.x else -30
            ax.add_patch(patches.Rectangle((px - (dx/2), py - 15), dx, 30, edgecolor='#1f77b4', facecolor='#1f77b4', zorder=3))

        ax.text(self.x, self.y + 15, self.name, fontsize=10, fontweight='bold', ha='center', va='center', zorder=5)
        ax.text(self.x, self.y - 15, f"{self.freq_ghz}GHz", fontsize=8, color='#023e8a', ha='center', va='center', zorder=5)


@dataclass
class CouplerEdge:
    name: str             # 커플러(연결선) 식별 이름 (예: "C01")
    q_start: str          # 시작 큐비트 이름
    p_start: str          # 시작 큐비트의 포트 이름 (예: "p_top_right")
    q_end: str            # 끝 큐비트 이름
    p_end: str            # 끝 큐비트의 포트 이름
    freq_ghz: float       # 이 커플러(공진기)의 목표 주파수 (GHz)
    orientation: str      # 미앤더 진행 방향 ('horizontal' / 'vertical')
    target_length_um: float = 0.0            # 계산된 목표 전체 길이 (um)
    waypoints: List[Tuple[float, float]] = None  # 실제 경로를 구성하는 좌표 목록


# ------------------------------------------------------------------------------
# 3. 렌더링 엔진 (Pure Python & Qiskit Metal)
# ------------------------------------------------------------------------------
def render_pure_python(
    qubits: Dict[str, QubitMacro],  # 큐비트 이름 -> QubitMacro 매핑
    edges: List[CouplerEdge],       # 그릴 커플러(연결선) 목록
    output_file: str                # 저장할 이미지 파일 경로
):
    # Matplotlib 기반 커스텀 미리보기
    fig, ax = plt.subplots(figsize=(12, 6))

    for q in qubits.values():
        q.draw_pure(ax)

    for e in edges:
        pts = e.waypoints
        # 꺾이는 지점(waypoint)마다 각진 모서리를 최대 30um 반경의 둥근 곡선(fillet)으로 대체
        path_verts, path_codes = [pts[0]], [Path.MOVETO]
        for i in range(1, len(pts) - 1):
            p_prev, p_curr, p_next = np.array(pts[i-1]), np.array(pts[i]), np.array(pts[i+1])
            v1, v2 = p_prev - p_curr, p_next - p_curr  # 현재 점 기준 이전/다음 방향 벡터
            l1, l2 = np.linalg.norm(v1), np.linalg.norm(v2)
            r = min(30.0, l1 / 2.0, l2 / 2.0)  # 인접 선분 길이를 넘지 않도록 라운드 반경 제한
            path_verts.extend([p_curr + (v1/l1)*r, p_curr, p_curr + (v2/l2)*r])
            path_codes.extend([Path.LINETO, Path.CURVE3, Path.CURVE3])
        path_verts.append(pts[-1])
        path_codes.append(Path.LINETO)

        cpw_path = Path(path_verts, path_codes)
        # CPW(코플레인 도파관)를 두 겹으로 그림: 연한 넓은 선(gap) 위에 진한 얇은 선(trace)
        ax.add_patch(patches.PathPatch(cpw_path, fill=False, edgecolor='#a6cee3', linewidth=14.0, alpha=0.6, zorder=2))
        ax.add_patch(patches.PathPatch(cpw_path, fill=False, edgecolor='blue', linewidth=6.0, alpha=0.8, zorder=3))

        mid_x, mid_y = pts[len(pts)//2]  # 경로 중앙에 이름·길이 라벨 표시
        ax.text(mid_x, mid_y + 40, f"{e.name}\n({e.target_length_um/1000:.2f}mm)", 
                fontsize=8, color='darkred', fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="red", lw=0.8, alpha=0.8), zorder=6)

    ax.set_aspect('equal')
    ax.set_title("Pure Python PnR View", fontsize=12)
    ax.set_xlabel("X Position (um)")
    ax.set_ylabel("Y Position (um)")
    ax.set_xlim(-400, 5000)
    ax.set_ylim(-600, 2200)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def render_qiskit_metal(
    qubits: Dict[str, QubitMacro],  # 큐비트 이름 -> QubitMacro 매핑
    edges: List[CouplerEdge],       # 연결할 커플러(연결선) 목록
    output_file: str                # 저장할 이미지 파일 경로
):
    # Qiskit Metal CAD 렌더링
    design = designs.DesignPlanar()

    # 큐비트마다 실제 CAD 컴포넌트(TransmonPocket) 배치, 4방향 커플링 핀(connection_pads) 함께 정의
    # loc_W(-1/+1: 좌/우), loc_H(-1/+1: 하/상)로 핀이 붙는 모서리 위치를 지정
    for q in qubits.values():
        TransmonPocket(design, q.name, options=dict(
            pos_x=f'{q.x}um', pos_y=f'{q.y}um',                     # 큐비트 중심 좌표
            pocket_width=f'{q.width}um', pocket_height=f'{q.height}um', # 바깥 Pocket(공동) 크기
            pad_width='200um', pad_height='80um', pad_gap='30um',   # 큐비트 자체 커패시터 Pad 크기/간격
            connection_pads=dict(
                # loc_W: 좌(-1)/우(+1), loc_H: 하(-1)/상(+1) 모서리 위치
                # pad_width/pad_gap: 각 커플링 핀 패드의 폭/간격
                p_top_left     = dict(loc_W=-1, loc_H=+1, pad_width='30um', pad_gap='10um'),
                p_bottom_left  = dict(loc_W=-1, loc_H=-1, pad_width='30um', pad_gap='10um'),
                p_top_right    = dict(loc_W=+1, loc_H=+1, pad_width='30um', pad_gap='10um'),
                p_bottom_right = dict(loc_W=+1, loc_H=-1, pad_width='30um', pad_gap='10um'),
            )
        ))

    # 커플러마다 두 큐비트 핀을 목표 길이(total_length)에 맞춘 미앤더 CPW 라인으로 연결
    for e in edges:
        RouteMeander(design, e.name, options=dict(
            pin_inputs=dict(
                start_pin=dict(component=e.q_start, pin=e.p_start),  # 연결 시작 큐비트/핀
                end_pin=dict(component=e.q_end, pin=e.p_end)         # 연결 끝 큐비트/핀
            ),
            trace_width='10um',                    # 도선(trace) 폭
            trace_gap='6um',                        # 도선과 그라운드 사이 간격(gap)
            fillet='30um',                          # 꺾이는 모서리 라운드 반경
            total_length=f'{e.target_length_um}um', # 맞춰야 할 전체 길이 (반파장 길이)
            meander=dict(spacing='120um', asymmetry='0um'),         # 지그재그 줄 간 간격 / 좌우 비대칭 오프셋
            lead=dict(start_straight='150um', end_straight='150um') # 시작/끝 직선 리드 구간 길이
        ))

    # Qiskit Metal이 생성한 도형 테이블들을 순회하며 하나의 그림에 렌더링
    fig, ax = plt.subplots(figsize=(12, 6))
    for _, table in design.qgeometry.tables.items():
        if not table.empty:
            table.plot(ax=ax, alpha=0.6, edgecolor='blue')

    ax.set_aspect('equal')
    ax.set_title("Qiskit Metal CAD View", fontsize=12)
    ax.set_xlabel("X Position (mm)")
    ax.set_ylabel("Y Position (mm)")
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()



















# ------------------------------------------------------------------------------
# 4. 메인 실행부
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    # 큐비트 5개 배치 (Q0-Q1-Q2-Q3 일렬 + Q2 위에 Q4, IBM 5큐비트 유사 토폴로지)
    qubits = {
        "Q0": QubitMacro("Q0", x=0.0,    y=0.0,    freq_ghz=5.0),
        "Q1": QubitMacro("Q1", x=1500.0, y=0.0,    freq_ghz=5.2),
        "Q2": QubitMacro("Q2", x=3000.0, y=0.0,    freq_ghz=4.9),
        "Q3": QubitMacro("Q3", x=4500.0, y=0.0,    freq_ghz=5.1),
        "Q4": QubitMacro("Q4", x=3000.0, y=1500.0, freq_ghz=5.3)
    }

    # 커플러 정의: (이름, 시작 큐비트, 시작 포트, 끝 큐비트, 끝 포트, 공진 주파수GHz, 배치 방향)
    raw_edges = [
        ("C01", "Q0", "p_top_right",    "Q1", "p_top_left",     6.8, 'horizontal'),
        ("C12", "Q1", "p_bottom_right", "Q2", "p_bottom_left",  6.5, 'horizontal'),
        ("C23", "Q2", "p_top_right",    "Q3", "p_top_left",     7.0, 'horizontal'),
        ("C24", "Q2", "p_top_left",     "Q4", "p_bottom_left",   6.6, 'vertical'),
    ]

    # 커플러별로 목표 반파장 길이를 계산하고, 그 길이에 맞는 미앤더 경로를 생성
    edges = []
    for name, q_s, p_s, q_e, p_e, f_c, orient in raw_edges:
        L_target = calc_half_wave_length(f_c)
        p_start_pos = qubits[q_s].ports[p_s]
        p_end_pos = qubits[q_e].ports[p_e]

        wps = generate_exact_meander_path(p_start_pos, p_end_pos, L_target, orientation=orient, meander_height=250.0)
        edges.append(CouplerEdge(name, q_s, p_s, q_e, p_e, f_c, orient, L_target, wps))

    # 두 가지 방식(순수 matplotlib 미리보기 / Qiskit Metal CAD)으로 각각 렌더링
    render_pure_python(qubits, edges, "ibm_5q_pure_python_fixed.png")
    render_qiskit_metal(qubits, edges, "ibm_5q_qiskit_metal_fixed.png")
    print("렌더링 완료: ibm_5q_pure_python_fixed.png / ibm_5q_qiskit_metal_fixed.png")

"""