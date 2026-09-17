"""DetailedPlacement: LG 산출물의 배치 bounding box 면적을 등방 축소로 압축한다.

목적은 면적 하나뿐이다. crosstalk 완화나 밀도 평탄화 같은 다른 목적을 같이 넣지
않는다 — 01_mainref Phase 0에서 제약과 목적을 한 스칼라에 섞었다가 hotspot이
4->18로 악화된 채 전량 롤백된 전례가 있다(이 저장소 docs에도 같은 교훈을 남긴다).
DP는 "면적을 최소화하되 DRC/임베딩을 절대 깨지 않는다"는 하나의 목적-제약 쌍만 갖는다.

방법: 큐빗+세그먼트+커플러 박스 전 좌표를, 배치 bbox 중심(cx, cy)에 대해 α배 축소하는
닮음변환(similarity transform) 하나로 압축한다. α=1.0(LG 산출물, 이미 legal)에서
이분탐색으로 α를 내려가며 DRC와 "세그먼트가 자기 박스 안"을 검사해 α_min(가장 작은
합격 α)을 찾는다.

왜 이분탐색이 성립하는가: 같은 중심에서 같은 비율로 모든 점을 미는 변환이므로, 두 점
사이 거리는 정확히 α배가 된다 — α가 작아질수록 모든 쌍의 거리가 단조 감소하므로,
"어떤 α에서 위반이면 그보다 작은 α에서도 위반"이 항상 성립한다(단조성 보장, 이분탐색의
전제조건). proper crossing은 더 강하게 보장된다: 닮음변환은 방향(orientation) 판정
자체를 뒤집지 않으므로 crossing 여부가 **모든 α>0에서 수학적으로 동일**하다(FP의
_scale_chip이 이미 쓰는 것과 같은 논리) — 그래서 crossing은 이분탐색 대상에 넣지 않고,
최종 α에서 한 번만 재확인한다(내부 버그 방지용 assertion — floorplan.py의 기존 스타일).

큐빗 크기(w/h)는 절대 축소하지 않는다 — TransmonPocket 실제 물리 치수라 알고리즘이
임의로 줄일 권한이 없다. 반면 세그먼트 "크기"(Coupler.segment_size_um)는 위치와 함께
축소한다 — 이유(격자 정렬 문제에 대한 판단, 아래 상세):

GP는 세그먼트를 segment_size_um 격자에 정확히 스냅해 배치했고, 실측 결과 이 저장소
6칩 전부에서 서로 다른 세그먼트 쌍 수백~수천 개가 정확히 한 칸(=segment_size_um) 거리로
"맞닿아" 있다(체인 연결성을 위해 GP가 의도적으로 인접 배치한 결과). 세그먼트 크기를
고정한 채 위치만 축소하면(선택지 B만 단독 적용), α<1이 되는 순간 이 맞닿은 쌍들이
전부 겹침으로 바뀐다 — 즉 α_min이 사실상 정확히 1.0으로 고정되어 모든 칩에서 압축이
전혀 안 되는 퇴화된 결과가 나온다(실측으로 확인함, 아래 검증 절 참고). segment_size_um
자체가 QPlacer에서 그대로 가져온 값으로 "물리적 유도 근거 없음"이라고 이 저장소에
이미 명시돼 있다(core/state.py의 Coupler.segment_size_um 주석) — 즉 세그먼트 "칸 크기"는
GP의 배치 알고리즘을 위한 이산화 단위일 뿐 고정된 물리 치수가 아니다. 그래서 위치와
함께 크기도 α배로 줄인다(선택지 C) — 같은 비율로 줄이므로 세그먼트-세그먼트 관계(맞닿음
포함)와 "세그먼트가 자기 박스 안"은 **어떤 α에서도 원래 상태와 완전히 동일하게 보존된다**
(둘 다 같은 닮음변환을 받으므로 — 증명은 아래 검증 절 참고, cc와 box-containment 위반은
전 범위 α에서 0으로 실측 확인됨). 반면 큐빗(고정 크기)과 세그먼트(축소 크기)가 섞이는
qc, 그리고 큐빗끼리인 qq는 실제로 α에 따라 달라지는 진짜 제약으로 남는다 — 이게 압축의
진짜 병목이 된다는 뜻이고, 실측 결과도 그렇게 나온다.

주의(한계로 남겨둠): segment_size_um을 줄이는 건 GP의 배치 이산화 단위를 줄이는 것이지,
meander_spacing_um(CPW 실제 폴드 피치, 어느 정도 물리적 근거가 있는 값)을 재검증하는
게 아니다 — 즉 DP는 "이 배치 장부(placement bookkeeping)를 더 촘촘히 다시 그릴 수
있는가"를 확인할 뿐, "실제 배선이 이만큼 촘촘하게 접힐 수 있는가"까지 물리적으로
검증하지는 않는다. 이 간극은 이 파이프라인의 segment 모델 자체가 처음부터 갖고 있던
근사(core/state.py Coupler 클래스 주석 참고)이고, DP가 새로 만든 문제가 아니다 — 다만
DP를 거치면 압축 폭만큼 그 근사가 더 벌어진다는 점은 후속 물리 검증이 필요할 때 참고할
것.

DP 이후 Coupler.segment_size_um은 더 이상 params.segment_size_um(전역값)과 같지 않다
(칩별로 α_min배 줄어든 값) — Coupler.num_segments 프로퍼티는 그 새 segment_size_um으로
재계산되므로 len(coupler.segments)(실제 배치된 개수, 불변)와 더 이상 일치하지 않을 수
있다. 이건 GP가 참조하는 "몇 개 놓아야 하는가" 질문에 대한 답이 DP 이후로는 더 이상
유효하지 않다는 뜻일 뿐(DP가 세그먼트를 추가/제거하지 않으므로 실제 배치 개수는 항상
GP가 정한 그대로다), 별도로 바로잡지 않는다.
"""
import logging
import time
from dataclasses import replace

import numpy as np

from core.state import AABB_EPS_UM, ChipState, Coupler, Qubit, Segment
from core.floorplan import count_crossings, sort_edges

# 이분탐색 하한(더 못 내려간다는 보장이 아니라 탐색 범위의 안전한 바닥) — qq만으로도
# 큐빗들이 중심 한 점으로 뭉개지기 전에 훨씬 일찍 걸리므로 실전에서 이 바닥에 닿을 일은
# 없다. 40회 반복이면 (1-FLOOR)/2^40 수준의 해상도라 실용적으로 충분하다.
_ALPHA_FLOOR = 0.01
_BISECT_ITERS = 40


class DetailedPlacementInfeasibleError(Exception):
    pass


class DetailedPlacement:
    def __init__(self, params):
        self.params = params
        self.skipped: list[tuple[str, str]] = []

    def run(self, states: list[ChipState]) -> list[ChipState]:
        out: list[ChipState] = []
        for state in states:
            try:
                out.append(self._place(state))
            except DetailedPlacementInfeasibleError as e:
                self.skipped.append((state.processor_name, str(e)))
        return out

    def _place(self, state: ChipState) -> ChipState:
        t0 = time.perf_counter()
        model = _LayoutModel.build(state)

        area_before = model.bbox_area(1.0)
        edges = sort_edges(state.cmap)
        pos_before = {i: (q.x, q.y) for i, q in state.qubits.items()}
        crossings_before = count_crossings(pos_before, edges, proper_only=True)

        if not model.feasible(1.0)[0]:
            # LG가 이미 legal을 보장했어야 한다 — 여기 걸리면 그 보장이 깨진 것.
            raise DetailedPlacementInfeasibleError(
                f"{state.processor_name}: DP 입력(LG 산출물)이 이미 DRC 위반을 갖고 있습니다 "
                "— LG의 legal 보장이 깨졌거나 DP 모델링에 불일치가 있습니다(내부 버그)."
            )

        alpha_min, first_violation = _bisect_alpha_min(model)

        pos_after_arr = model.qubit_positions(alpha_min)
        pos_after = {qid: (float(pos_after_arr[0][i]), float(pos_after_arr[1][i]))
                     for i, qid in enumerate(model.qubit_ids)}
        crossings_after = count_crossings(pos_after, edges, proper_only=True)
        if crossings_after != crossings_before:
            raise RuntimeError(
                f"{state.processor_name}: 닮음변환 후 proper crossing이 {crossings_before} -> "
                f"{crossings_after}로 변했습니다 — 수학적으로 불가능해야 하는 일이라 내부 버그입니다."
            )

        new_state = model.apply(state, alpha_min)
        area_after = model.bbox_area(alpha_min)
        dt = time.perf_counter() - t0

        n_unplaced = sum(1 for c in state.couplers.values() if c.num_segments > 0 and not c.segments)
        logging.info(
            "[DP] %s: alpha_min=%.4f bbox_area %.3e -> %.3e um^2 (%.1f%% 감소) | "
            "첫 위반 제약(그 직전 α에서)=%s | num_unplaced_couplers=%d | %.3fs",
            state.processor_name, alpha_min, area_before, area_after,
            100.0 * (1.0 - area_after / area_before) if area_before > 0 else 0.0,
            first_violation, n_unplaced, dt,
        )

        return new_state


# ---------------------------------------------------------------------------
# 레이아웃 모델 — 큐빗/세그먼트/박스를 numpy 배열로 들고, 임의의 α에서 변환된 좌표와
# 위반 개수를 빠르게 계산한다(이분탐색 40회 x 6칩을 감당할 만큼 벡터화).
# ---------------------------------------------------------------------------

class _LayoutModel:
    def __init__(
        self,
        qubit_ids, qx, qy, qhw, qhh,
        seg_keys, seg_idx, seg_coupler_idx, sx, sy, seg_half0,
        box_keys, box_x0, box_x1, box_y0, box_y1,
        own_mask,  # own_mask[i, j] = segment i의 커플러가 qubit j를 소유하는지
        qubit_id_to_col: dict[int, int],
        coupler_key_to_box_row: dict[tuple[int, int], int],
        chip_w: float, chip_h: float,
        cx: float, cy: float,
    ):
        self.qubit_ids = qubit_ids
        self.qx, self.qy, self.qhw, self.qhh = qx, qy, qhw, qhh
        self.seg_keys = seg_keys
        self.seg_idx = seg_idx
        self.seg_coupler_idx = seg_coupler_idx
        self.sx, self.sy, self.seg_half0 = sx, sy, seg_half0
        self.box_keys = box_keys
        self.box_x0, self.box_x1, self.box_y0, self.box_y1 = box_x0, box_x1, box_y0, box_y1
        self.own_mask = own_mask
        self.qubit_id_to_col = qubit_id_to_col
        self.coupler_key_to_box_row = coupler_key_to_box_row
        self.chip_w, self.chip_h = chip_w, chip_h
        self.cx, self.cy = cx, cy

    @classmethod
    def build(cls, state: ChipState) -> "_LayoutModel":
        qubit_ids = sorted(state.qubits)
        qx = np.array([state.qubits[i].x for i in qubit_ids], dtype=float)
        qy = np.array([state.qubits[i].y for i in qubit_ids], dtype=float)
        qhw = np.array([state.qubits[i].w / 2.0 for i in qubit_ids], dtype=float)
        qhh = np.array([state.qubits[i].h / 2.0 for i in qubit_ids], dtype=float)
        qubit_id_to_col = {qid: i for i, qid in enumerate(qubit_ids)}

        seg_keys: list[tuple[int, int]] = []
        seg_idx_list: list[int] = []
        seg_coupler_idx: list[int] = []
        sx_list: list[float] = []
        sy_list: list[float] = []
        seg_half_list: list[float] = []
        coupler_ordinal: dict[tuple[int, int], int] = {}
        for ordinal, key in enumerate(sorted(state.couplers)):
            coupler_ordinal[key] = ordinal
        for key, c in state.couplers.items():
            half = c.segment_size_um / 2.0
            for s in c.segments:
                seg_keys.append(key)
                seg_idx_list.append(s.idx)
                seg_coupler_idx.append(coupler_ordinal[key])
                sx_list.append(s.x)
                sy_list.append(s.y)
                seg_half_list.append(half)
        sx = np.array(sx_list, dtype=float)
        sy = np.array(sy_list, dtype=float)
        seg_half0 = np.array(seg_half_list, dtype=float)
        seg_coupler_idx_arr = np.array(seg_coupler_idx, dtype=int)

        box_keys = sorted(state.coupler_regions)
        coupler_key_to_box_row = {k: i for i, k in enumerate(box_keys)}
        box_x0 = np.array([state.coupler_regions[k][0] for k in box_keys], dtype=float)
        box_x1 = np.array([state.coupler_regions[k][1] for k in box_keys], dtype=float)
        box_y0 = np.array([state.coupler_regions[k][2] for k in box_keys], dtype=float)
        box_y1 = np.array([state.coupler_regions[k][3] for k in box_keys], dtype=float)

        # own_mask: 세그먼트 i의 소속 커플러가 큐빗 j를 소유하는지(자기 큐빗 예외, qc에서
        # 제외 대상) — α와 무관하므로 한 번만 계산해 둔다.
        n_seg, n_q = len(seg_keys), len(qubit_ids)
        own_mask = np.zeros((n_seg, n_q), dtype=bool)
        for i, key in enumerate(seg_keys):
            coupler = state.couplers[key]
            for j, qid in enumerate(qubit_ids):
                if coupler.is_own_qubit(qid):
                    own_mask[i, j] = True

        all_x0 = list(qx - qhw) + list(sx - seg_half0)
        all_x1 = list(qx + qhw) + list(sx + seg_half0)
        all_y0 = list(qy - qhh) + list(sy - seg_half0)
        all_y1 = list(qy + qhh) + list(sy + seg_half0)
        cx = (min(all_x0) + max(all_x1)) / 2.0 if all_x0 else 0.0
        cy = (min(all_y0) + max(all_y1)) / 2.0 if all_y0 else 0.0

        return cls(
            qubit_ids, qx, qy, qhw, qhh,
            seg_keys, seg_idx_list, seg_coupler_idx_arr, sx, sy, seg_half0,
            box_keys, box_x0, box_x1, box_y0, box_y1,
            own_mask, qubit_id_to_col, coupler_key_to_box_row,
            state.chip_width, state.chip_height, cx, cy,
        )

    # α배 축소된 큐빗 좌표 (2, n_qubit) — 크기(qhw/qhh)는 절대 안 바뀐다.
    def qubit_positions(self, alpha: float):
        nx = self.cx + alpha * (self.qx - self.cx)
        ny = self.cy + alpha * (self.qy - self.cy)
        return nx, ny

    # α배 축소된 세그먼트 좌표 + α배 축소된 세그먼트 반폭 (크기도 함께 줄인다 — 모듈
    # docstring의 "격자 정렬 문제" 판단 참고).
    def segment_positions(self, alpha: float):
        nx = self.cx + alpha * (self.sx - self.cx)
        ny = self.cy + alpha * (self.sy - self.cy)
        nhalf = alpha * self.seg_half0
        return nx, ny, nhalf

    def box_bounds(self, alpha: float):
        x0 = self.cx + alpha * (self.box_x0 - self.cx)
        x1 = self.cx + alpha * (self.box_x1 - self.cx)
        y0 = self.cy + alpha * (self.box_y0 - self.cy)
        y1 = self.cy + alpha * (self.box_y1 - self.cy)
        return x0, x1, y0, y1

    # 큐빗+세그먼트 전체를 감싸는 bounding box 면적 — DP가 최소화하는 목적함수.
    def bbox_area(self, alpha: float) -> float:
        qx, qy = self.qubit_positions(alpha)
        sx, sy, shalf = self.segment_positions(alpha)
        x0 = min(np.min(qx - self.qhw), np.min(sx - shalf) if len(sx) else np.inf)
        x1 = max(np.max(qx + self.qhw), np.max(sx + shalf) if len(sx) else -np.inf)
        y0 = min(np.min(qy - self.qhh), np.min(sy - shalf) if len(sy) else np.inf)
        y1 = max(np.max(qy + self.qhh), np.max(sy + shalf) if len(sy) else -np.inf)
        return float((x1 - x0) * (y1 - y0))

    # α에서 위반 여부 + 어느 제약이 먼저 걸렸는지(qq/qc/cc/die/box 순으로 검사, 이 순서
    # 자체가 우선순위는 아니고 그냥 보고용 — 개수만 중요하다). (feasible, 위반명) 반환.
    def feasible(self, alpha: float) -> tuple[bool, str | None]:
        qx, qy = self.qubit_positions(alpha)
        sx, sy, shalf = self.segment_positions(alpha)
        bx0, bx1, by0, by1 = self.box_bounds(alpha)

        # 전부 AABB_EPS_UM(core/state.py) 하나로 통일한다 — func/compute.py·core/
        # legalization.py의 _aabb_overlap과 정확히 같은 기준이어야, 여기서 "feasible"이라고
        # 판단한 α가 최종적으로 legalization.py의 _find_qq_overlaps 등(권위 있는 단일
        # 정의)으로 재검증했을 때도 반드시 통과한다. 한때 qq/qc는 eps=0(엄격), cc만 eps
        # 적용으로 나눠서 짰었는데, 그러면 cc는 무해한 부동소수 잡음을 피하면서 qq/qc는
        # legalization.py 쪽(그때는 eps 없이 완전 엄격)과 정확히 같은 값을 요구하게 되고,
        # 그 결과 여기서 계산한 dx가 legalization.py가 재계산한 dx와 극미하게(~1e-13)
        # 달라 "여기선 안전, 저기선 위반"이라는 모순이 실제로 발생했다(grid_25 큐빗
        # 23/24 쌍 실측). AABB_EPS_UM을 두 파일 모두에 같은 값으로 박아 넣은 지금은 그
        # 모순이 구조적으로 불가능하다 — qq/qc(진짜 임계값)든 cc(맞닿음 보존, 모듈
        # docstring 참고)든 같은 잣대로 재는 게 맞다.
        if _pairwise_overlap_count(qx, qy, self.qhw, self.qhh, self.qhw, self.qhh, eps=AABB_EPS_UM) > 0:
            return False, "qq"

        if len(sx) and _cross_overlap_count(sx, sy, shalf, shalf, qx, qy, self.qhw, self.qhh,
                                             self.own_mask, eps=AABB_EPS_UM) > 0:
            return False, "qc"

        if len(sx) and _pairwise_overlap_count(sx, sy, shalf, shalf, shalf, shalf,
                                                self.seg_coupler_idx, eps=AABB_EPS_UM) > 0:
            return False, "cc"

        if np.any(qx - self.qhw < -AABB_EPS_UM) or np.any(qx + self.qhw > self.chip_w + AABB_EPS_UM) \
                or np.any(qy - self.qhh < -AABB_EPS_UM) or np.any(qy + self.qhh > self.chip_h + AABB_EPS_UM):
            return False, "die(qubit)"
        if len(sx) and (np.any(sx - shalf < -AABB_EPS_UM) or np.any(sx + shalf > self.chip_w + AABB_EPS_UM)
                        or np.any(sy - shalf < -AABB_EPS_UM) or np.any(sy + shalf > self.chip_h + AABB_EPS_UM)):
            return False, "die(segment)"

        if len(sx):
            box_row = np.array([self.coupler_key_to_box_row[k] for k in self.seg_keys])
            in_box = (
                (sx - shalf >= bx0[box_row] - AABB_EPS_UM) & (sx + shalf <= bx1[box_row] + AABB_EPS_UM)
                & (sy - shalf >= by0[box_row] - AABB_EPS_UM) & (sy + shalf <= by1[box_row] + AABB_EPS_UM)
            )
            if not np.all(in_box):
                return False, "segment_in_box"

        return True, None

    def apply(self, state: ChipState, alpha: float) -> ChipState:
        qx, qy = self.qubit_positions(alpha)
        new_qubits = dict(state.qubits)
        for i, qid in enumerate(self.qubit_ids):
            q = new_qubits[qid]
            new_qubits[qid] = replace(q, x=float(qx[i]), y=float(qy[i]))

        sx, sy, _shalf = self.segment_positions(alpha)
        per_coupler_new_segs: dict[tuple[int, int], list[Segment]] = {}
        for i, key in enumerate(self.seg_keys):
            new_seg = Segment(idx=self.seg_idx[i], x=float(sx[i]), y=float(sy[i]))
            per_coupler_new_segs.setdefault(key, []).append(new_seg)

        new_couplers = {}
        for key, c in state.couplers.items():
            segs = per_coupler_new_segs.get(key)
            new_couplers[key] = replace(
                c,
                segments=segs if segs is not None else [],
                segment_size_um=alpha * c.segment_size_um,
            )

        bx0, bx1, by0, by1 = self.box_bounds(alpha)
        new_regions = {
            k: (float(bx0[i]), float(bx1[i]), float(by0[i]), float(by1[i]))
            for i, k in enumerate(self.box_keys)
        }

        return replace(state, qubits=new_qubits, couplers=new_couplers, coupler_regions=new_regions)


# ---------------------------------------------------------------------------
# 벡터화된 겹침 카운트 (numpy) — 이분탐색 40회 x 최대 ~1300 세그먼트(eagle)를 감당하려면
# 순수 파이썬 이중 루프로는 느리다(느낀 O(n^2) 반복이 40회 곱해짐). floorplan.py의
# 브로드캐스트 스타일을 그대로 따른다.
# ---------------------------------------------------------------------------

# 같은 집합 안에서 서로 다른 원소 쌍의 AABB 겹침 개수. group_ids를 주면 같은 그룹(같은
# 커플러)끼리는 세지 않는다(cc는 "다른 커플러 간" 정의라서). eps는 항상 AABB_EPS_UM을
# 넘긴다 — func/compute.py·core/legalization.py의 _aabb_overlap과 같은 기준이어야 하는
# 이유는 _LayoutModel.feasible()의 호출부 주석 참고.
def _pairwise_overlap_count(xs, ys, hws, hhs, hws2, hhs2, group_ids=None, eps: float = 0.0) -> int:
    n = len(xs)
    if n < 2:
        return 0
    dx = np.abs(xs[:, None] - xs[None, :])
    dy = np.abs(ys[:, None] - ys[None, :])
    sum_hw = hws[:, None] + hws2[None, :]
    sum_hh = hhs[:, None] + hhs2[None, :]
    overlap = (dx < sum_hw - eps) & (dy < sum_hh - eps)
    iu = np.triu_indices(n, k=1)
    mask = overlap[iu]
    if group_ids is not None:
        same_group = group_ids[iu[0]] == group_ids[iu[1]]
        mask = mask & ~same_group
    return int(np.count_nonzero(mask))


# 서로 다른 두 집합(세그먼트 vs 큐빗) 사이 AABB 겹침 개수. exclude_mask[i,j]=True면 그
# 쌍은(자기 큐빗) 원천적으로 위반 대상에서 뺀다. eps는 항상 AABB_EPS_UM.
def _cross_overlap_count(xs1, ys1, hws1, hhs1, xs2, ys2, hws2, hhs2, exclude_mask, eps: float = 0.0) -> int:
    if len(xs1) == 0 or len(xs2) == 0:
        return 0
    dx = np.abs(xs1[:, None] - xs2[None, :])
    dy = np.abs(ys1[:, None] - ys2[None, :])
    sum_hw = hws1[:, None] + hws2[None, :]
    sum_hh = hhs1[:, None] + hhs2[None, :]
    overlap = (dx < sum_hw - eps) & (dy < sum_hh - eps) & ~exclude_mask
    return int(np.count_nonzero(overlap))


# ---------------------------------------------------------------------------
# 이분탐색
# ---------------------------------------------------------------------------

def _bisect_alpha_min(model: _LayoutModel) -> tuple[float, str | None]:
    lo, hi = _ALPHA_FLOOR, 1.0
    ok_floor, violation_at_floor = model.feasible(lo)
    if ok_floor:
        # 바닥까지도 압축 가능 — 이 저장소 규모에서 실제로 나올 일은 없지만(큐빗이 중심
        # 한 점으로 뭉개지기 훨씬 전에 qq가 걸림), 수학적으로는 가능한 입력을 방어한다.
        return lo, None

    first_violation = violation_at_floor
    for _ in range(_BISECT_ITERS):
        mid = (lo + hi) / 2.0
        ok, violation = model.feasible(mid)
        if ok:
            hi = mid
        else:
            lo = mid
            first_violation = violation
    return hi, first_violation
