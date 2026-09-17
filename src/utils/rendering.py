import os

import matplotlib
matplotlib.use('Agg')  # headless — main.py는 디스플레이 없는 환경(서버)에서도 돌아가야 한다
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb

from core.state import ChipState

# 요소가 die 밖으로 나갔을 때 강조할 색 — 포트 강조(주황)와 구분되도록 별도로 예약.
_OOB_COLOR = "#d62728"


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

        # FP 박스를 먼저 그려 맨 아래 깔아 둔다 — 세그먼트/bbox()가 그 위에 겹쳐 보여야
        # "박스 안에 세그먼트가 머무는지"를 한눈에 비교할 수 있다.
        for key in sorted(chip.couplers):
            box = chip.coupler_regions.get(key)
            if box is not None:
                has_boxes = True
                self._draw_coupler_box(ax, chip, box)

        for key in sorted(chip.couplers):
            coupler = chip.couplers[key]
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
        if chip.fallback_used:
            title += " | fallback"
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.3)
        ax.margins(0.04)
        ax.set_aspect("equal", adjustable="box")

        self._draw_legend(ax, has_segments, has_boxes, has_bbox)

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

    def _draw_segments(self, ax, coupler):
        color = self._coupler_color(coupler.id)
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
