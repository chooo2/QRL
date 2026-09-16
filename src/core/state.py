import math
from dataclasses import dataclass, field

@dataclass
class Qubit:
    id: int         # node id
    w: float        # width
    h: float        # height
    f: float        # frequency ghz
    x: float = 0.0  # x-coordinate(center)
    y: float = 0.0  # y-coordinate(center)

    @property
    def ports(self):
        off_x = self.w / 2.0
        off_y = 180.0
        return {
            "top_left":     (self.x - off_x, self.y + off_y),
            "top_right":    (self.x + off_x, self.y + off_y),
            "bottom_left":  (self.x - off_x, self.y - off_y),
            "bottom_right": (self.x + off_x, self.y - off_y),
        }

@dataclass
class Coupler:
    id: int     # edge id
    q1: int     # start qubit id
    q2: int     # end qubit id
    f: float    # edge(coupler) frequency ghz — half-wave 길이 공식의 입력
    p: float    # meander pitch (um)

    @property
    def l(self):
        v_phase = 1.3e8
        return (v_phase / (2.0 * self.f * 1e9)) * 1e6

    def boundary(self, q1: Qubit, q2: Qubit) -> tuple[float, float, float, float]:
        long_side = math.sqrt(self.l * self.p)
        short_side = self.p
        dx, dy = q2.x - q1.x, q2.y - q1.y
        half_w, half_h = (short_side / 2.0, long_side / 2.0) if abs(dy) > abs(dx) \
            else (long_side / 2.0, short_side / 2.0)
        cx, cy = (q1.x + q2.x) / 2.0, (q1.y + q2.y) / 2.0
        return (cx - half_w, cx + half_w, cy - half_h, cy + half_h)

@dataclass
class ChipState:
    processor_name: str
    num_qubits: int
    cmap: tuple  # coupling map: ((q1, q2), ...)

    qubits: dict[int, Qubit] = field(default_factory=dict)
    couplers: dict[tuple[int, int], Coupler] = field(default_factory=dict)

    chip_width: float = 0.0
    chip_height: float = 0.0

    # 결함 등으로 비활성화된 큐빗을 걸러내기 위한 마스크. None이면 전체가 active.
    active_qubits: set[int] | None = None
