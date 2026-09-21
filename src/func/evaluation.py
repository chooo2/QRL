import csv
import os

from core.state import ChipState
from func.compute import DEFAULT_DETUNING_GHZ, compute_metrics
from func.metric import Metric


# 단계별(5절 검증표 기준) Performance | Design Rule 지표 목록. 값 자체는 전부
# func/compute.py(→ func/metric.py의 Metric)가 이미 계산해둔 필드를 그대로 읽을 뿐,
# 여기서 새 지표를 계산하지 않는다 — 같은 지표를 두 번 구현하지 않기 위함.
STAGE_METRICS: dict[str, list[str]] = {
    "FloorPlanning": [
        "area_utilization",
        "edge_cross_point",
    ],

    "GlobalPlacement": [
        "area_utilization",
        "freq_hotspot_proportion",
        "num_hotspot_qubits",
        "drc_qq_overlap",
        "drc_qc_overlap",
        "drc_cc_overlap",
        "drc_edge_cross",
    ],

    "Legalization": [
        "freq_hotspot_proportion",
        "num_hotspot_qubits",
        "drc_qq_overlap",
        "drc_qc_overlap",
        "drc_cc_overlap",
        "drc_edge_cross",
    ],

    "DetailedPlacement": [
        "area_utilization",
        "freq_hotspot_proportion",
        "num_hotspot_qubits",
        "fidelity",
        "drc_qq_overlap",
        "drc_qc_overlap",
        "drc_cc_overlap",
        "drc_edge_cross",
    ],

    # Coupler Crossing = 실제 라우팅된 배선끼리의 교차 (compute_route_cross_point,
    # RT 이전 단계에서는 의미가 없어 이 단계에만 둔다).
    "Routing": [
        "route_cross_point",
    ],
}

METRIC_NAMES: dict[str, str] = {
    "area_utilization": "Area Utilization (%)",
    "freq_hotspot_proportion": "Frequency Hotspot Proportion (%)",
    "num_hotspot_qubits": "Number of Hotspot Qubits",
    # Fidelity: 이 저장소엔 계산 코드가 없다. QPlacer(Eq.13)는 회로 시뮬레이션이 있어야
    # 정의되는 값이라 이 파이프라인(배치/라우팅까지만 다룸) 범위 밖으로 판단했다 —
    # Metric에 fidelity 필드 자체가 없어 getattr(..., None)이 항상 None을 반환하고,
    # format_value가 그걸 그대로 "-"로 남긴다(0이나 추정치로 채우지 않는다).
    "fidelity": "Fidelity",
    "edge_cross_point": "Edge Crossing",
    "route_cross_point": "Coupler Crossing",
    "drc_qq_overlap": "Qubit-Qubit Overlap",
    "drc_qc_overlap": "Qubit-Coupler Overlap",
    "drc_cc_overlap": "Coupler-Coupler Overlap",
    "drc_edge_cross": "Edge Crossing",
}


def _format_value(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "PASS" if value else "FAIL"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


# 배치된 칩 상태(ChipState)를 평가 지표(Metric)로 계산하는 최상위 클래스.
# 실제 계산은 func/compute.py의 함수들이 담당하고, 이 클래스는 그것들을 감싸는
# 진입점 역할만 한다 — main.py는 이 클래스 하나만 임포트해서 쓰면 된다.
class Evaluation:

    # 물리/평가 파라미터를 저장 (예: detuning_ghz 등 향후 보정치)
    def __init__(self, params=None, detuning_ghz: float = DEFAULT_DETUNING_GHZ):
        self.params = params
        self.detuning_ghz = detuning_ghz

    # 칩 상태 하나를 모든 지표(Metric)로 계산해 반환
    def compute(self, state: ChipState, active: set[int] | None = None) -> Metric:
        return compute_metrics(state, active=active, detuning_ghz=self.detuning_ghz)

    # 상태 하나를 단일 스칼라 점수로 축약 (낮을수록 좋음) — 주파수 hotspot 비율
    def score(self, state: ChipState) -> float:
        return self.compute(state).freq_hotspot_proportion

    # 후보 상태들 중 점수(score)가 가장 낮은 최적 상태와 그 지표를 선택
    def pick_best(self, candidates: list[ChipState]) -> tuple[ChipState | None, Metric | None]:
        best_state, best_metric, best_score = None, None, None
        for candidate in candidates:
            metric = self.compute(candidate)
            if best_score is None or metric.freq_hotspot_proportion < best_score:
                best_state, best_metric, best_score = candidate, metric, metric.freq_hotspot_proportion
        return best_state, best_metric


    # 한 단계(stage)의 표를 (header, rows)로 만든다 — 콘솔 출력과 CSV 저장이 이 결과를
    # 그대로 공유한다(표 구성 로직을 두 번 만들지 않기 위함). rows의 각 값은 이미
    # _format_value로 문자열화되어 있다.
    @staticmethod
    def _stage_rows(evaluation_results, benchmark_names, stage):
        header = ["Metric", *benchmark_names, "Mean"]
        rows = []

        for metric_name in STAGE_METRICS[stage]:
            row = [METRIC_NAMES[metric_name]]
            values = []

            for benchmark in benchmark_names:
                metric = evaluation_results.get(stage, {}).get(benchmark)
                value = getattr(metric, metric_name, None) if metric is not None else None

                row.append(_format_value(value))

                if value is not None and not isinstance(value, bool):
                    values.append(value)

            row.append(_format_value(sum(values) / len(values)) if values else "-")
            rows.append(row)

        return header, rows

    def print_evaluation_table(self, evaluation_results, benchmark_names):
        metric_width = 40
        benchmark_width = 12
        mean_width = 12
        widths = [metric_width, *([benchmark_width] * len(benchmark_names)), mean_width]

        for stage in STAGE_METRICS:
            header, rows = self._stage_rows(evaluation_results, benchmark_names, stage)

            print(f"\n{'=' * 20}")
            print(stage)
            print(f"{'=' * 20}")

            print(" | ".join(f"{column:>{width}}" for column, width in zip(header, widths)))
            print("-+-".join("-" * width for width in widths))

            for row in rows:
                print(" | ".join(f"{value:>{width}}" for value, width in zip(row, widths)))

    # 단계별 CSV를 output/eval/<Stage>.csv로 저장한다 — 콘솔 표와 같은 내용을 표 옮겨
    # 적기(transcription) 용도로 남긴다. 행=지표, 열=칩(+Mean)은 콘솔 표와 동일하다.
    def save_csv(self, evaluation_results, benchmark_names, out_dir: str = "output/eval"):
        os.makedirs(out_dir, exist_ok=True)

        for stage in STAGE_METRICS:
            header, rows = self._stage_rows(evaluation_results, benchmark_names, stage)

            path = os.path.join(out_dir, f"{stage}.csv")
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(header)
                writer.writerows(rows)