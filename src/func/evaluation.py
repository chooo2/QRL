from core.state import ChipState
from func.compute import DEFAULT_DETUNING_GHZ, compute_metrics
from func.metric import Metric


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
