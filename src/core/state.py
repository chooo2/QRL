import math
from dataclasses import dataclass, field

@dataclass
class Qubit:
    id: int         # node id
    w:  float        # width
    h:  float        # height
    f:  float        # frequency ghz
    x:  float = 0.0  # x-coordinate(center)
    y:  float = 0.0  # y-coordinate(center)

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

    def __str__(self):
        return f"id={self.id:>3d} w={self.w:>6.2f} h={self.h:>6.2f} "f"f={self.f:>5.2f}GHz x={self.x:>8.2f} y={self.y:>8.2f}"
@dataclass
class Coupler:
    id: int      # edge id
    q1: int      # start qubit id
    q2: int      # end qubit id
    f:  float    # edge(coupler) frequency ghz
    p:  float    # meander pitch (um)

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

    def __str__(self):
        return f"id={self.id:>3d} q1={self.q1:>3d} q2={self.q2:>3d} " f"f={self.f:>5.2f}GHz p={self.p:>6.2f}um"

@dataclass
class ChipState:
    processor_name: str
    num_qubits: int
    cmap: tuple  # coupling map: ((q1, q2), ...)
    chip_width: float
    chip_height: float

    qubits: dict[int, Qubit] = field(default_factory=dict)
    couplers: dict[tuple[int, int], Coupler] = field(default_factory=dict)

    @classmethod
    def processor_config(cls, config: dict, params) -> "ChipState":
        num_qubits = config["num_qubits"]
        cmap = tuple(config["coupling_map"])
        # chip_width = config.get("chip_width_um") or params.chip_width
        # chip_height = config.get("chip_height_um") or chip_width
        chip_width = config.get("chip_width_um")
        chip_height = config.get("chip_height_um")
    
        freq_ghz = config["freq_ghz"]
        qubits = {
            i: Qubit(id=i, w=params.qubit_width, h=params.qubit_height, f=freq_ghz[i])
            for i in range(num_qubits)
        }

        edge_freq = {(min(u, v), max(u, v)): f for u, v, f in config["edge_freq_ghz"]}
        couplers = {
            (q1, q2): Coupler(id=idx, q1=q1, q2=q2, f=edge_freq[(q1, q2)], p=params.coupler_meander_pitch_um)
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
