import logging
from dataclasses import dataclass, asdict

metric_label: dict[str, tuple[str, str]] = {
    "freq_hotspot_proportion":  ("Frequency Hotspot Proportion", ""),
    "num_hotspot_qubits":       ("Number of Hotspot Qubits", ""),
    "fidelity":                 ("Fidelity (TODO: 계산식 미정)", ""),
    "area_utilization":         ("Area Efficiency", ""),
    "edge_cross_point":         ("Edge Cross Point", ""),
    "coupler_cross_point":      ("Coupler Cross Point", ""),

    "drc_qq_overlap":           ("Qubit-Qubit Overlap", ""),
    "drc_qc_overlap":           ("Qubit-Coupler Overlap", ""),
    "drc_cc_overlap":           ("Coupler-Coupler Overlap", ""),
    "drc_edge_cross":           ("Edge Cross", "")
}

@dataclass
class Metric:
    # Performance
    freq_hotspot_proportion: float  # 주파수 hotspot 큐빗 비율
    num_hotspot_qubits: int         # 주파수 hotspot에 관여한 큐빗 수
    area_utilization: float         # 칩 면적 대비 배치 효율
    edge_cross_point: int           # 커플링맵 논리 엣지(직선 근사) 교차 쌍 개수
    coupler_cross_point: int        # 커플러 경계 상자(AABB) 겹침 쌍 개수

    # Design Rule Check - True = 위반 없음(통과)
    drc_qq_overlap: bool            # 큐빗-큐빗 겹침 없음 여부
    drc_qc_overlap: bool            # 큐빗-커플러 겹침 없음 여부
    drc_cc_overlap: bool            # 커플러-커플러 겹침 없음 여부 (coupler_cross_point == 0)
    drc_edge_cross: bool            # 논리 엣지 교차 없음 여부 (edge_cross_point == 0)

    fidelity: float | None = None   # Fidelity (계산식 미정 — compute_fidelity 참고)

    def dict(self):
        return asdict(self)

    def log_result(self):
        logging.getLogger(__name__).info('=' * 70)
        logging.getLogger(__name__).info("Processor: %s", str)
        for f, (disp, unit) in metric_label.items():
            v = getattr(self, f)
            val_str = f"{v:.4g}{unit}" if isinstance(v, float) else str(v)
            logging.getLogger(__name__).info("%-45s %s", disp, val_str)
        logging.getLogger(__name__).info('=' * 70)