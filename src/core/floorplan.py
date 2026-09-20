"""Floorplan: 커플링 그래프 평면성 판정 + 교차-0 큐빗 좌표 생성 + 최소 간격 확보 + 포트 배정.

/home/LabMember/ngchoi/research/01_mainref/src/core/globalplacement.py(16칩에서 교차 0 검증됨)를
이 저장소의 ChipState 스키마에 맞게 가져왔다. 그 파일의 구조 설명을 그대로 옮긴다:

주 경로(Floorplan._primary_layout)는 spectral/tree-aware/ring-aware 초기배치 + 힘 기반 정제로
사후 검증 후 채택하는 필터일 뿐 교차 0을 수학적으로 보장하지 않는다. 이 모듈은 그 경로가
실패했을 때만 쓰이는 결정론적 보장 폴백을 함께 제공한다:

  1. 평면성 판정 (nx.check_planarity) — 비평면이면 단일 레이어에서 교차 0은 애초에 불가능.
  2. 비평면 처리 (handle_non_planar) — strict: 예외, crossover: 최대 평면 부분그래프(MPS) +
     잘린 엣지를 크로스오버(air bridge) 대상으로 반환.
  3. 보장된 배치 (guaranteed_embedding) — Chrobak-Payne 알고리즘(nx.combinatorial_
     embedding_to_pos)으로 평면 임베딩을 직선 배치로 변환 — 3-연결을 요구하지 않는 일반
     평면 그래프에 대해 교차 0을 보장.
  4. 국소 정제 (impred_refine) — ImPrEd 계열: 각 정점의 이동을 "교차를 만들지 않는" 최대
     크기로 클리핑한 뒤 힘을 적용, 임베딩을 불변으로 유지하며 엣지 길이를 균일화.

이 저장소로 옮기며 의도적으로 뺀 것: 원본의 커플러 박스(boundary) 관련 코드 전부
(_coupler_refinement, _eliminate_overlaps, build_coupler_boundary, _coupler_overlap,
widen_thin_coupler_spans, target_span_by_edge 기반 _scale_chip 분기). 커플러를 세그먼트로
쪼개 배치하는 건 GP의 일이다(core/state.py의 Coupler.segments/region()/bbox() 참고) — FP는
큐빗 좌표와 포트 배정까지만 확정한다. self를 쓰지 않는 함수는 전부 모듈 레벨에 뒀다(원본도
이미 그렇게 정리돼 있었음).

FP/GP 책임 분리 (2026-09-16, aspen_11 FP 실패 조사 이후 명문화):
  FP의 간격 책임은 "큐빗 몸체(footprint)끼리 물리적으로 겹치지 않는다"는 것 하나뿐이다
  (min_qubit_spacing_um, 아래 참고). 라우팅(커플러 배선)에 필요한 여유 공간은 FP 책임이
  아니다 — FP는 애초에 커플러 세그먼트/폭(meander_spacing_um 등)을 모르므로 그 여유를
  만들어 줄 방법이 없다. 01_mainref는 FP(원문 GlobalPlacement._primary_layout) 안에서
  커플러 박스를 알고 있는 _coupler_refinement/_eliminate_overlaps가 큐빗 간격도 함께
  벌려줬는데(트레이싱 결과: aspen_11에서 이 두 함수가 min_pairwise_distance를 317um ->
  680um -> 최종 424um으로 끌어올림 — 방금 위에서 뺐다고 적은 바로 그 코드다), 이 저장소는
  그 코드를 의도적으로 뺐으므로 그만큼의 여유가 사라졌다. 이 저장소 설계에서는 그 여유를
  GP가 확보해야 한다 — GP가 Coupler.region()으로 세그먼트 배치 후보 영역을 계산할 때 필요한
  라우팅 공간이 부족하면, 그건 GP 단계의 실패(또는 GP가 큐빗을 추가로 밀어내는 등의 조치)로
  다뤄야지 FP의 min_qubit_spacing_um을 부풀려서 미리 여유를 만들어두면 안 된다 — 그러면
  "무슨 목적의 여유인지"가 다시 섞인다.
"""
import logging
import math
from collections import deque, defaultdict
from dataclasses import replace

import networkx as nx
import numpy as np

from core.state import ChipState, Coupler, Qubit


# strict 모드에서 비평면 그래프를 만났을 때 발생 — 단일 레이어 교차 0이 불가능함을 뜻함.
class NonPlanarError(Exception):
    pass


# 교차 0 직선 배치는 얻었으나, 현재 die 크기에서 최소 큐빗 간격을 만족시키며 담을 수 없거나
# (포트 4개 한도 안에서) 포트 배정이 불가능함 — 이 칩은 이 크기/토폴로지로는 단일 레이어
# fab-ready floorplan이 불가능하다. Floorplan.run()이 이 칩을 skipped에 기록하고 계속한다.
class PlacementInfeasibleError(Exception):
    pass


# ---------------------------------------------------------------------------
# 그래프 유틸
# ---------------------------------------------------------------------------

# 커플링 맵으로부터 큐빗별 인접 집합(인접 리스트)을 구성
def build_adj(edges: list[tuple[int, int]], num_qubits: int) -> dict[int, set[int]]:
    adj: dict[int, set] = {i: set() for i in range(num_qubits)}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    return adj


# 커플링 맵의 엣지를 중복 제거 후 정렬된 리스트로 변환
def sort_edges(cmap: tuple) -> list[tuple[int, int]]:
    edges = {(min(u, v), max(u, v)) for u, v in cmap if u != v}
    return sorted(edges)


# ---------------------------------------------------------------------------
# 선분 교차 판정
# ---------------------------------------------------------------------------

def segments_cross(
    p1: tuple[float, float], p2: tuple[float, float],
    p3: tuple[float, float], p4: tuple[float, float],
) -> bool:
    o1 = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])
    o2 = (p2[0] - p1[0]) * (p4[1] - p1[1]) - (p2[1] - p1[1]) * (p4[0] - p1[0])
    o3 = (p4[0] - p3[0]) * (p1[1] - p3[1]) - (p4[1] - p3[1]) * (p1[0] - p3[0])
    o4 = (p4[0] - p3[0]) * (p2[1] - p3[1]) - (p4[1] - p3[1]) * (p2[0] - p3[0])
    return o1 * o2 < 0.0 and o3 * o4 < 0.0


# 공용 교차 판정(벡터화) — segments_cross가 놓치는 공선/T자 접촉도 잡는다. eps는 각 엣지 쌍
# 자신의 길이 기준(로컬)이어야 함 — 전역 bbox 기준이면 무관한 먼 노드가 움직여도 판정이 뒤집힘.
# include_touch=False면 proper(X자로 서로를 가로지르는) 교차만 — 공선/끝점 접촉 제외.
def _crossing_mask(
    p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray, rel_tol: float,
    include_touch: bool = True,
) -> np.ndarray:
    def orient(o: np.ndarray, p: np.ndarray, q: np.ndarray) -> np.ndarray:
        return (p[:, 0] - o[:, 0]) * (q[:, 1] - o[:, 1]) - (p[:, 1] - o[:, 1]) * (q[:, 0] - o[:, 0])

    def seg_len(p: np.ndarray, q: np.ndarray) -> np.ndarray:
        return np.maximum(np.hypot(q[:, 0] - p[:, 0], q[:, 1] - p[:, 1]), 1e-12)

    len12 = seg_len(p1, p2)
    len34 = seg_len(p3, p4)
    eps1 = rel_tol * len12
    eps2 = rel_tol * len34

    o1 = orient(p1, p2, p3) / len12
    o2 = orient(p1, p2, p4) / len12
    o3 = orient(p3, p4, p1) / len34
    o4 = orient(p3, p4, p2) / len34

    proper = (
        (((o1 > eps1) & (o2 < -eps1)) | ((o1 < -eps1) & (o2 > eps1)))
        & (((o3 > eps2) & (o4 < -eps2)) | ((o3 < -eps2) & (o4 > eps2)))
    )
    if not include_touch:
        return proper

    def on_segment(p: np.ndarray, q: np.ndarray, r: np.ndarray, eps: np.ndarray) -> np.ndarray:
        return (
            (np.minimum(p[:, 0], q[:, 0]) - eps <= r[:, 0]) & (r[:, 0] <= np.maximum(p[:, 0], q[:, 0]) + eps)
            & (np.minimum(p[:, 1], q[:, 1]) - eps <= r[:, 1]) & (r[:, 1] <= np.maximum(p[:, 1], q[:, 1]) + eps)
        )

    touch = (
        ((np.abs(o1) <= eps1) & on_segment(p1, p2, p3, eps1))
        | ((np.abs(o2) <= eps1) & on_segment(p1, p2, p4, eps1))
        | ((np.abs(o3) <= eps2) & on_segment(p3, p4, p1, eps2))
        | ((np.abs(o4) <= eps2) & on_segment(p3, p4, p2, eps2))
    )
    return proper | touch


# 커플링맵 전체에서 서로 교차하는 엣지 쌍의 개수 계산 (dict 좌표, 공차 포함 — 보장 경로/
# 사후검증용). proper_only=True면 X자로 서로를 가로지르는 교차만("proper crossing").
def count_crossings(
    pos: dict[int, tuple[float, float]], edges, rel_tol: float = 1e-9, proper_only: bool = False,
) -> int:
    uniq = sorted({(min(u, v), max(u, v)) for u, v in edges if u != v and u in pos and v in pos})
    E = len(uniq)
    if E < 2:
        return 0
    edge_arr = np.array(uniq, dtype=int)
    a, b = edge_arr[:, 0], edge_arr[:, 1]
    pi, pj = np.triu_indices(E, k=1)
    ai, bi, aj, bj = a[pi], b[pi], a[pj], b[pj]
    keep = ~((ai == aj) | (ai == bj) | (bi == aj) | (bi == bj))
    ai, bi, aj, bj = ai[keep], bi[keep], aj[keep], bj[keep]
    if len(ai) == 0:
        return 0

    def gather(idx: np.ndarray) -> np.ndarray:
        return np.array([pos[int(i)] for i in idx], dtype=float)

    p1, p2, p3, p4 = gather(ai), gather(bi), gather(aj), gather(bj)
    return int(np.count_nonzero(_crossing_mask(p1, p2, p3, p4, rel_tol, include_touch=not proper_only)))


# edges_a x edges_b 쌍의 교차 개수(정점 공유 쌍 제외) — v 이동 시 인접 엣지의 새 교차만 싸게 확인.
def count_crossings_between(
    pos: dict[int, tuple[float, float]], edges_a, edges_b, rel_tol: float = 1e-9,
) -> int:
    if not edges_a or not edges_b:
        return 0
    A = np.array(edges_a, dtype=int)
    B = np.array(edges_b, dtype=int)
    ai = np.repeat(A[:, 0], len(B)); bi = np.repeat(A[:, 1], len(B))
    aj = np.tile(B[:, 0], len(A));   bj = np.tile(B[:, 1], len(A))
    keep = ~((ai == aj) | (ai == bj) | (bi == aj) | (bi == bj))
    ai, bi, aj, bj = ai[keep], bi[keep], aj[keep], bj[keep]
    if len(ai) == 0:
        return 0

    def gather(idx: np.ndarray) -> np.ndarray:
        return np.array([pos[int(i)] for i in idx], dtype=float)

    p1, p2, p3, p4 = gather(ai), gather(bi), gather(aj), gather(bj)
    return int(np.count_nonzero(_crossing_mask(p1, p2, p3, p4, rel_tol, include_touch=True)))


# 엣지 쌍 간 선분 교차 개수(배열 좌표, 공차 없음 — _primary_layout 후보 스코어링에서
# candidate_cnt번 반복 호출되므로 dict 기반 count_crossings보다 싼 배열 버전을 따로 둔다).
def _cnt_edge_crossing(pos: np.ndarray, edges: list[tuple[int, int]]) -> int:
    E = len(edges)
    if E < 2:
        return 0
    edge_arr = np.array(edges, dtype=int)
    a, b = edge_arr[:, 0], edge_arr[:, 1]
    pi, pj = np.triu_indices(E, k=1)
    ai, bi, aj, bj = a[pi], b[pi], a[pj], b[pj]
    keep = ~((ai == aj) | (ai == bj) | (bi == aj) | (bi == bj))
    ai, bi, aj, bj = ai[keep], bi[keep], aj[keep], bj[keep]

    p1, p2, p3, p4 = pos[ai], pos[bi], pos[aj], pos[bj]
    o1 = (p2[:, 0]-p1[:, 0])*(p3[:, 1]-p1[:, 1]) - (p2[:, 1]-p1[:, 1])*(p3[:, 0]-p1[:, 0])
    o2 = (p2[:, 0]-p1[:, 0])*(p4[:, 1]-p1[:, 1]) - (p2[:, 1]-p1[:, 1])*(p4[:, 0]-p1[:, 0])
    o3 = (p4[:, 0]-p3[:, 0])*(p1[:, 1]-p3[:, 1]) - (p4[:, 1]-p3[:, 1])*(p1[:, 0]-p3[:, 0])
    o4 = (p4[:, 0]-p3[:, 0])*(p2[:, 1]-p3[:, 1]) - (p4[:, 1]-p3[:, 1])*(p2[:, 0]-p3[:, 0])
    crossing = (o1 * o2 < 0.0) & (o3 * o4 < 0.0)
    return int(np.count_nonzero(crossing))


# ---------------------------------------------------------------------------
# 최소 간격
# ---------------------------------------------------------------------------

# 최소 쌍거리 (중심간, Chebyshev/L-infinity 거리) — FP의 유일한 간격 책임인 "큐빗 몸체끼리
# 물리적으로 안 겹침"을 판정하는 정확한 메트릭이다. 정사각형(qubit_width == qubit_height)
# 큐빗 두 개가 안 겹칠 필요충분조건은 max(|dx|, |dy|) >= 한 변 길이다 — 한 축에서만 그만큼
# 떨어져 있어도 겹치지 않는다("OR" 조건). 유클리드 거리(hypot)는 이 조건과 다르다: 예를 들어
# dx=dy=qubit_width/sqrt(2)인 정확히 대각선 배치는 유클리드 거리로는 딱 qubit_width라
# "안 겹침"으로 통과시키지만 실제로는 두 정사각형이 겹친다 — 반대로 축 정렬 배치는 필요
# 이상으로 엄격하게 본다. 2026-09-16 FP 책임을 "물리적 비중첩"으로 명확히 좁히면서
# (core/floorplan.py 모듈 docstring 참고) 유클리드에서 Chebyshev로 바꿨다.
def min_pairwise_distance(pos: dict) -> float:
    P = np.array([pos[i] for i in sorted(pos)], dtype=float)
    if len(P) < 2:
        return float("inf")
    dx = np.abs(P[:, None, 0] - P[None, :, 0])
    dy = np.abs(P[:, None, 1] - P[None, :, 1])
    d = np.maximum(dx, dy)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


# min_sep 미달 쌍이 있으면 die(margin 포함) 안에 들어가는 한도에서 배치를 균일 확대한다.
# 균일 상사변환이라 교차/비교차 구조를 정확히 보존한다. 확대해도 여전히 미달이면 그대로
# 반환 — 호출부가 min_pairwise_distance로 재확인해 PlacementInfeasibleError를 던진다.
def uniform_expand_to_spacing(
    pos: dict, min_sep: float, chip_w: float, chip_h: float, margin: float,
) -> dict:
    keys = sorted(pos)
    P = np.array([pos[k] for k in keys], dtype=float)
    if len(P) < 2:
        return pos
    d = min_pairwise_distance(pos)
    if d >= min_sep or d < 1e-9:
        return pos
    # 배치 자신의 bbox 중심을 기준으로 확대한 뒤 die 중심에 놓는다 (코너에 몰린 배치도
    # 정상 확대되도록 — die 중심 기준으로 재면 ext가 커져 fit이 과소평가됨).
    lo, hi = P.min(axis=0), P.max(axis=0)
    layout_center = (lo + hi) / 2.0
    rel = P - layout_center
    ext = np.max(np.abs(rel), axis=0)
    ext[ext < 1e-9] = 1e-9
    fit = min((chip_w / 2.0 - margin) / ext[0], (chip_h / 2.0 - margin) / ext[1])
    factor = min(min_sep / d, max(fit, 1e-9))
    die_center = np.array([chip_w / 2.0, chip_h / 2.0])
    Q = die_center + rel * factor
    return {k: (float(Q[i, 0]), float(Q[i, 1])) for i, k in enumerate(keys)}


# ---------------------------------------------------------------------------
# 평면성 판정 + 비평면 처리
# ---------------------------------------------------------------------------

# 그래프 구조만으로 낸 근사 좌표를 MPS 증분 삽입의 "엣지 길이" 정렬 기준으로 사용 —
# 이 시점엔 아직 물리 좌표가 없으므로(평면성 판정은 FP가 배치를 계산하기 전에 함).
def _edge_length_proxy(G: nx.Graph) -> dict:
    try:
        return nx.spring_layout(G, seed=0)
    except Exception:
        return {n: (0.0, 0.0) for n in G.nodes()}


# 짧은 엣지부터 증분 삽입해 평면성을 유지하는 최대 부분그래프 (단순 길이순 정렬 휴리스틱)
def maximal_planar_subgraph(G: nx.Graph) -> tuple[nx.Graph, list[tuple[int, int]]]:
    pos = _edge_length_proxy(G)

    def edge_len(e: tuple[int, int]) -> float:
        u, v = e
        return math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])

    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    cut: list[tuple[int, int]] = []
    for u, v in sorted(G.edges(), key=edge_len):
        H.add_edge(u, v)
        if not nx.check_planarity(H)[0]:
            H.remove_edge(u, v)
            cut.append((min(u, v), max(u, v)))
    return H, cut


# 비평면 그래프 처리. mode='strict'는 예외를 던져 처리 불가로 분류하고,
# mode='crossover'는 MPS를 추출해 (평면 부분그래프, 그 임베딩, 잘린 엣지 목록)을 반환한다.
def handle_non_planar(
    G: nx.Graph, mode: str,
) -> tuple[nx.Graph, "nx.PlanarEmbedding", list[tuple[int, int]]]:
    if mode == "strict":
        raise NonPlanarError(
            f"Coupling graph is non-planar ({G.number_of_nodes()} qubits, "
            f"{G.number_of_edges()} edges) — a single-layer zero-crossing layout "
            "is mathematically impossible. Use floorplan_planarity_mode='crossover' to "
            "extract a maximal planar subgraph and route the remaining edges as "
            "crossovers (air bridges) instead."
        )
    if mode != "crossover":
        raise ValueError(f"unknown non-planar handling mode: {mode!r}")

    H, cut = maximal_planar_subgraph(G)
    ok, emb = nx.check_planarity(H)
    if not ok:
        raise RuntimeError("maximal_planar_subgraph produced a non-planar graph (internal bug)")
    return H, emb, cut


# ---------------------------------------------------------------------------
# 보장된 배치 (연속 좌표 — FP 좌표가 연속값임을 전제)
# ---------------------------------------------------------------------------

# 컴포넌트가 여럿이면(예: 고립 큐빗) 대표 정점끼리 가상 엣지로 이어 연결 그래프로 만든다.
# 가상 엣지는 반환 좌표에만 반영, 교차 판정에서는 제외해야 함(호출부 책임).
def _connect_components_virtually(G: nx.Graph) -> nx.Graph:
    components = [sorted(c) for c in nx.connected_components(G)]
    if len(components) <= 1:
        return G
    components.sort(key=lambda c: (-len(c), c[0]))
    reps = [max(c, key=lambda q: (G.degree(q), -q)) for c in components]
    G_aug = G.copy()
    for i in range(len(reps) - 1):
        G_aug.add_edge(reps[i], reps[i + 1])
    return G_aug


# 평면 임베딩 -> 교차-0 직선 배치 (Chrobak-Payne, nx.combinatorial_embedding_to_pos, O(n)).
# Tutte는 3-연결 그래프에서만 교차 0 보장 — heavy-hex의 degree-2 큐빗 때문에 실측 교차 발생 → 폐기.
# Chrobak-Payne은 임의 단순 평면 그래프에서 교차 0 (정수 격자 좌표지만 연속 좌표로 그대로 사용).
def guaranteed_embedding(G: nx.Graph) -> dict[int, tuple[float, float]]:
    nodes = list(G.nodes())
    n = len(nodes)
    if n == 0:
        return {}
    if n == 1:
        return {nodes[0]: (0.0, 0.0)}
    if G.number_of_edges() == 0:
        return {
            v: (math.cos(2.0 * math.pi * i / n), math.sin(2.0 * math.pi * i / n))
            for i, v in enumerate(nodes)
        }

    G_aug = _connect_components_virtually(G)
    ok, emb = nx.check_planarity(G_aug)
    if not ok:
        raise RuntimeError(
            "component-connecting virtual edges made the graph non-planar — should be "
            "impossible (joining two planar-embedded regions with one edge is always "
            "planarity-preserving); indicates a bug upstream"
        )
    return nx.combinatorial_embedding_to_pos(emb)


# ---------------------------------------------------------------------------
# ImPrEd 스타일 국소 정제 (단순 버전 — 쿼드트리 등 성능 최적화는 실측 전에는 하지 않음)
# ---------------------------------------------------------------------------

# 정점 v를 delta만큼 옮길 때 새 교차를 안 만드는 최대 이동 비율 t∈[0,1] (이분탐색 근사)
def _safe_move_fraction(
    v: int, delta: tuple[float, float],
    pos: dict[int, tuple[float, float]],
    edges_incident_v: list[tuple[int, int]],
    edges_other: list[tuple[int, int]],
    max_bisect: int = 10,
) -> float:
    if not edges_incident_v or not edges_other:
        return 1.0
    if abs(delta[0]) < 1e-12 and abs(delta[1]) < 1e-12:
        return 0.0

    # v를 비율 t만큼 옮겼을 때 생기는 교차 개수. pos를 통째 복사하지 않고 v만
    # 임시 이동 후 복원한다 (127q+ 폴백에서 dict(pos) 전체복사가 병목이었음).
    _v_orig = pos[v]

    def crossings_at(t: float) -> int:
        pos[v] = (_v_orig[0] + t * delta[0], _v_orig[1] + t * delta[1])
        try:
            return count_crossings_between(pos, edges_incident_v, edges_other)
        finally:
            pos[v] = _v_orig

    if crossings_at(1.0) == 0:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(max_bisect):
        mid = (lo + hi) * 0.5
        if crossings_at(mid) == 0:
            lo = mid
        else:
            hi = mid
    return lo


# ImPrEd 계열 국소 정제: 엣지길이->중앙값 인력 + 정점 반발, 이동은 새 교차 안 만드는 크기로
# 클리핑. 각도 균일화 항은 없음(반발력으로 간접).
def impred_refine(
    pos: dict[int, tuple[float, float]],
    edges: list[tuple[int, int]],
    iterations: int = 60,
    step_frac: float = 0.15,
) -> dict[int, tuple[float, float]]:
    edges = [(min(u, v), max(u, v)) for u, v in edges if u != v]
    nodes = list(pos.keys())
    if len(edges) < 2 or len(nodes) < 2:
        return dict(pos)

    lengths = [math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1]) for u, v in edges]
    target_len = float(np.median(lengths)) if lengths else 1.0
    if target_len <= 1e-12:
        target_len = 1.0
    step = target_len * step_frac

    incident = {v: [] for v in nodes}
    for u, v in edges:
        incident[u].append((u, v))
        incident[v].append((u, v))
    edges_other_of = {v: [e for e in edges if v not in e] for v in nodes}

    cur = dict(pos)
    for _ in range(iterations):
        any_move = False
        for v in nodes:
            fx, fy = 0.0, 0.0
            vx, vy = cur[v]
            for u in cur:
                if u == v:
                    continue
                ux, uy = cur[u]
                dx, dy = vx - ux, vy - uy
                dist = math.hypot(dx, dy)
                if dist < 1e-9:
                    fx += (np.random.default_rng(hash((v, u)) & 0xffffffff).uniform(-1, 1))
                    fy += (np.random.default_rng(hash((u, v)) & 0xffffffff).uniform(-1, 1))
                    continue
                is_edge = (min(u, v), max(u, v)) in edges
                if is_edge:
                    # 인력: 목표보다 길면 v를 u 쪽으로 당겨 줄이고, 짧으면 밀어 늘린다.
                    fmag = (target_len - dist) / target_len
                    fx += fmag * dx / dist
                    fy += fmag * dy / dist
                else:
                    # 반발: 가까운 비인접 정점끼리만(붕괴 방지), 목표 길이 안쪽에서만 작용
                    if dist < target_len:
                        fmag = (target_len - dist) / target_len
                        fx += fmag * dx / dist
                        fy += fmag * dy / dist

            fnorm = math.hypot(fx, fy)
            if fnorm < 1e-9:
                continue
            delta = (fx / fnorm * step, fy / fnorm * step)
            t = _safe_move_fraction(v, delta, cur, incident[v], edges_other_of[v])
            if t > 1e-6:
                cur[v] = (cur[v][0] + t * delta[0], cur[v][1] + t * delta[1])
                any_move = True

        if not any_move:
            break

    return cur


# ---------------------------------------------------------------------------
# 주 경로: 토폴로지-인식 초기배치
# ---------------------------------------------------------------------------

# 라플라시안 고유벡터 기반으로 초기 큐빗 좌표를 계산
def _spectral_initial_coordinate(adj: dict[int, set[int]], edges: list[tuple[int, int]], rng: np.random.Generator) -> np.ndarray:
    n = len(adj)
    laplacian = np.zeros((n, n), dtype=float)
    for u, v in edges:
        laplacian[u, u] += 1.0
        laplacian[v, v] += 1.0
        laplacian[u, v] = -1.0
        laplacian[v, u] = -1.0

    try:
        _, vector = np.linalg.eigh(laplacian)
        if n >= 3: pos = vector[:, 1:3].astype(float)
        else: pos = np.column_stack((vector[:, 1], np.zeros(n)))
    except np.linalg.LinAlgError:
        angle = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
        pos = np.column_stack((np.cos(angle), np.sin(angle)))

    pos += rng.normal(0, 1e-4, size=pos.shape)

    pos = pos - np.mean(pos, axis=0)
    span = np.ptp(pos, axis=0)
    span[span < 1e-9] = 1.0
    return pos / span


# 큐빗 번호가 8개 단위 링을 이루는 옥타곤 격자 구간을 탐지
def _detect_octagon_rings(edges: list[tuple[int, int]], n: int,
                          ring_size: int = 8, min_hit: int = 7) -> list[list[int]]:
    if n < ring_size or n % ring_size != 0: return []
    edge_set = set(edges)
    rings = []
    for k in range(n // ring_size):
        base = k * ring_size
        ring = list(range(base, base + ring_size))
        ring_edges = {(min(ring[i], ring[(i + 1) % ring_size]), max(ring[i], ring[(i + 1) % ring_size]))
                     for i in range(ring_size)}
        if len(ring_edges & edge_set) >= min_hit:
            rings.append(ring)
    return rings


# 탐지된 링들이 전체 큐빗의 충분한 비율을 커버하는지 확인
def _ring_coverage_sufficient(rings: list[list[int]], n: int, min_coverage: float = 0.7) -> bool:
    if not rings: return False
    return sum(len(r) for r in rings) >= min_coverage * n


# 큐빗 -> 링 유닛을 하나로 뭉친 축약 그래프 노드 id 매핑을 생성
def _ring_super_map(n: int, rings: list[list[int]]) -> tuple[dict[int, int], int]:
    qubit_to_super: dict[int, int] = {}
    for ridx, ring in enumerate(rings):
        for q in ring:
            qubit_to_super[q] = ridx
    next_id = len(rings)
    for q in range(n):
        if q not in qubit_to_super:
            qubit_to_super[q] = next_id
            next_id += 1
    return qubit_to_super, next_id


# k각형 유닛의 로컬 정점 오프셋 템플릿을 생성
def _polygon_template(k: int, aspect: float, phase: float = 0.0) -> np.ndarray:
    angles = np.linspace(0.0, 2.0 * math.pi, k, endpoint=False) + phase
    pts = np.column_stack((np.cos(angles), np.sin(angles) * aspect))
    width = float(np.ptp(pts[:, 0]))
    return pts / max(width, 1e-9)


# 축약 그래프의 가중치 최단경로 행렬을 계산
def _weighted_shortest_path(coarse_adj: dict[int, set[int]], n_super: int,
                            edge_weight: dict[tuple[int, int], float] | None = None) -> np.ndarray:
    big = float(n_super) * 1e6
    dist = np.full((n_super, n_super), big, dtype=float)
    np.fill_diagonal(dist, 0.0)
    for u in range(n_super):
        for v in coarse_adj[u]:
            w = float((edge_weight or {}).get((min(u, v), max(u, v)), 1.0))
            dist[u, v] = min(dist[u, v], w)
    for mid in range(n_super):
        dist = np.minimum(dist, dist[:, mid:mid + 1] + dist[mid:mid + 1, :])
    if np.any(dist >= big):
        raise ValueError("coarse ring graph is disconnected — _primary_layout's "
                          "len(components)>1 gate should have routed this input away from the ring-aware path")
    return dist


# 축약 그래프(링 단위)를 MDS와 스프링 완화로 배치
def _coarse_ring_layout(coarse_adj: dict[int, set[int]], coarse_edges: list[tuple[int, int]],
                        n_super: int, k_coarse: float, rng: np.random.Generator, iterations: int = 150,
                        edge_weight: dict[tuple[int, int], float] | None = None) -> np.ndarray:
    if n_super <= 1: return np.zeros((max(n_super, 1), 2))
    if not coarse_edges:
        angle = np.linspace(0.0, 2.0 * math.pi, n_super, endpoint=False)
        return np.column_stack((np.cos(angle), np.sin(angle))) * k_coarse

    dist = _weighted_shortest_path(coarse_adj, n_super, edge_weight)
    if edge_weight is None:
        dist = dist * k_coarse
    d2 = dist * dist
    centering = np.eye(n_super) - np.ones((n_super, n_super)) / n_super
    b = -0.5 * centering @ d2 @ centering
    w, v = np.linalg.eigh(b)
    order = np.argsort(w)[::-1]
    w, v = w[order], v[:, order]
    w2 = np.clip(w[:2], 0.0, None)
    pos = v[:, :2] * np.sqrt(w2)

    if not np.any(np.abs(pos) > 1e-9):
        angle = np.linspace(0.0, 2.0 * math.pi, n_super, endpoint=False)
        pos = np.column_stack((np.cos(angle), np.sin(angle))) * k_coarse

    pos = pos + rng.normal(0, 1e-4, size=pos.shape)

    edge_array = np.array(coarse_edges, dtype=int)
    target = np.array([float((edge_weight or {}).get((min(u, v), max(u, v)), k_coarse))
                       for u, v in coarse_edges], dtype=float)

    temperature = k_coarse * 0.5
    for step in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]
        dd = np.linalg.norm(delta, axis=2) + 1e-9
        repulsion = ((k_coarse * k_coarse * 0.3) / dd)[:, :, None] * delta / dd[:, :, None]
        disp = np.sum(repulsion, axis=1)

        edge_delta = pos[edge_array[:, 0]] - pos[edge_array[:, 1]]
        edge_dist = np.linalg.norm(edge_delta, axis=1) + 1e-9
        attraction = (((edge_dist - target) * edge_dist / k_coarse)[:, None] * edge_delta / edge_dist[:, None])
        np.add.at(disp, edge_array[:, 0], -attraction)
        np.add.at(disp, edge_array[:, 1], attraction)

        disp_norm = np.linalg.norm(disp, axis=1) + 1e-9
        pos += (disp / disp_norm[:, None]) * np.minimum(disp_norm, temperature)[:, None]
        pos -= np.mean(pos, axis=0)
        temperature *= 1.0 - ((step + 1) / iterations) * 0.05

    return pos


# 옥타곤 링 구조를 강체 단위로 배치한 뒤 각 링에 로컬 모양을 얹어 초기 좌표를 생성
def _ring_aware_initial_coordinate(adj: dict[int, set[int]], edges: list[tuple[int, int]],
                                   rings: list[list[int]], edge_target: dict[tuple[int, int], float],
                                   base: float, rng: np.random.Generator, aspect: float = 1.6) -> np.ndarray:
    n = len(adj)
    qubit_to_super, n_super = _ring_super_map(n, rings)

    coarse_edge_set = set()
    for u, v in edges:
        su, sv = qubit_to_super[u], qubit_to_super[v]
        if su != sv:
            coarse_edge_set.add((min(su, sv), max(su, sv)))
    coarse_edges = sorted(coarse_edge_set)
    coarse_adj: dict[int, set[int]] = {i: set() for i in range(n_super)}
    for u, v in coarse_edges:
        coarse_adj[u].add(v); coarse_adj[v].add(u)

    def _ring_scale(ring: list[int]) -> tuple[np.ndarray, float]:
        k = len(ring)
        tmpl = _polygon_template(k, aspect)
        chord = float(np.linalg.norm(tmpl[1] - tmpl[0]))
        targets = [edge_target.get((min(ring[i], ring[(i + 1) % k]), max(ring[i], ring[(i + 1) % k])), base)
                  for i in range(k)]
        avg_target = float(np.mean(targets)) if targets else base
        return tmpl, avg_target / max(chord, 1e-9)

    ring_templates: dict[int, tuple[np.ndarray, float]] = {}
    vertex_radius: dict[int, float] = {}
    ring_extents = []
    for ridx, ring in enumerate(rings):
        tmpl, scale = _ring_scale(ring)
        ring_templates[ridx] = (tmpl, scale)
        ring_extents.append(float(np.max(np.ptp(tmpl * scale, axis=0))))
        for i, q in enumerate(ring):
            vertex_radius[q] = float(np.linalg.norm(tmpl[i] * scale))

    k_coarse = max(float(np.median(ring_extents)) * 2.2 if ring_extents else base * 2.2, base * 1.5)

    bridge_edges: dict[tuple[int, int], list[tuple[int, int, float]]] = defaultdict(list)
    for u, v in edges:
        su, sv = qubit_to_super[u], qubit_to_super[v]
        if su != sv:
            key = (min(su, sv), max(su, sv))
            node_lo, node_hi = (u, v) if su < sv else (v, u)
            bridge_edges[key].append((node_lo, node_hi, edge_target.get((min(u, v), max(u, v)), base)))

    coarse_edge_weight: dict[tuple[int, int], float] = {}
    for e in coarse_edges:
        items = bridge_edges.get(e, [(None, None, base)])
        r_lo = float(np.mean([vertex_radius.get(n0, 0.0) for n0, n1, t in items]))
        r_hi = float(np.mean([vertex_radius.get(n1, 0.0) for n0, n1, t in items]))
        gap = float(np.mean([t for _, _, t in items])) / math.sqrt(max(len(items), 1))
        weight = (r_lo + r_hi) * 1.3 + gap
        if weight <= r_lo + r_hi:
            raise ValueError(
                f"coarse edge weight({weight}) <= 충돌 하한선({r_lo + r_hi}) for ring pair {e} "
                "— 계수(1.3, sqrt 나눗셈)를 건드렸다면 이 하한선 위로 유지되는지 확인할 것"
            )
        coarse_edge_weight[e] = weight

    coarse_pos = _coarse_ring_layout(coarse_adj, coarse_edges, n_super, k_coarse, rng,
                                      edge_weight=coarse_edge_weight)

    pos = np.zeros((n, 2), dtype=float)
    for ridx, ring in enumerate(rings):
        tmpl, scale = ring_templates[ridx]
        center = coarse_pos[ridx]
        for i, q in enumerate(ring):
            pos[q] = center + tmpl[i] * scale

    for q in range(n):
        su = qubit_to_super[q]
        if su >= len(rings):
            pos[q] = coarse_pos[su] + rng.normal(0, base * 0.05, size=2)

    pos += rng.normal(0, 1e-4, size=pos.shape)
    return pos


# 그래프가 사이클 없는 스패닝 트리인지 확인하고 루트(최대 차수 노드)를 반환
def _detect_tree_root(edges: list[tuple[int, int]], n: int, components: list[list[int]]) -> int | None:
    if n < 4 or len(components) != 1 or len(edges) != n - 1:
        return None
    adj: dict[int, set[int]] = {i: set() for i in range(n)}
    for u, v in edges:
        adj[u].add(v); adj[v].add(u)
    return max(range(n), key=lambda q: (len(adj[q]), -q))


# 트리를 BFS로 레벨화해 방사형(radial) 초기 좌표를 생성
def _tree_aware_initial_coordinate(edges: list[tuple[int, int]], root: int,
                                   edge_target: dict[tuple[int, int], float], base: float,
                                   rng: np.random.Generator, n: int) -> np.ndarray:
    tree_adj: dict[int, set[int]] = {i: set() for i in range(n)}
    for u, v in edges:
        tree_adj[u].add(v); tree_adj[v].add(u)

    parent: dict[int, int | None] = {root: None}
    radius: dict[int, float] = {root: 0.0}
    order = [root]
    queue = deque([root])
    while queue:
        node = queue.popleft()
        for nbr in sorted(tree_adj[node]):
            if nbr not in parent:
                parent[nbr] = node
                edge_len = edge_target.get((min(node, nbr), max(node, nbr)), base)
                radius[nbr] = radius[node] + edge_len
                order.append(nbr)
                queue.append(nbr)

    children: dict[int, list[int]] = defaultdict(list)
    for node in order[1:]:
        children[parent[node]].append(node)

    leaf_count: dict[int, int] = {}
    for node in reversed(order):
        kids = children.get(node, [])
        leaf_count[node] = sum(leaf_count[c] for c in kids) if kids else 1

    angle_lo = {root: 0.0}
    angle_hi = {root: 2.0 * math.pi}
    for node in order:
        kids = children.get(node, [])
        if not kids:
            continue
        lo, hi = angle_lo[node], angle_hi[node]
        total_leaves = leaf_count[node]
        cursor = lo
        for c in kids:
            share = (hi - lo) * (leaf_count[c] / total_leaves)
            angle_lo[c] = cursor
            angle_hi[c] = cursor + share
            cursor += share

    pos = np.zeros((n, 2), dtype=float)
    for node in order[1:]:
        theta = 0.5 * (angle_lo[node] + angle_hi[node])
        r = radius[node]
        pos[node] = (r * math.cos(theta), r * math.sin(theta))

    pos += rng.normal(0, 1e-4, size=pos.shape)
    return pos


# 스프링(힘 기반) 모델로 큐빗 좌표를 반복적으로 완화
def _spring_refine_coordinate(pos: np.ndarray, edges: list[tuple[int, int]],
                              edge_target: dict[tuple[int, int], float] | None = None, iterations: int = 350,
                              max_temperature: float = 0.35) -> np.ndarray:
    n = len(pos)
    area = 4.0
    temperature = max_temperature
    k = math.sqrt(area / n)
    edge_array = np.array(edges, dtype=int)
    target = np.array([float((edge_target or {}).get((min(u, v), max(u, v)), k)) for u, v in edges], dtype=float)

    for step in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]
        dist = np.linalg.norm(delta, axis=2) + 1e-9
        repulsion = (k * k / dist)[:, :, None] * delta / dist[:, :, None]
        disp = np.sum(repulsion, axis=1)

        edge_delta = pos[edge_array[:, 0]] - pos[edge_array[:, 1]]
        edge_dist = np.linalg.norm(edge_delta, axis=1) + 1e-9
        attraction = (((edge_dist - target) * edge_dist / k)[:, None] * edge_delta / edge_dist[:, None])

        np.add.at(disp, edge_array[:, 0], -attraction)
        np.add.at(disp, edge_array[:, 1], attraction)

        disp_norm = np.linalg.norm(disp, axis=1) + 1e-9
        pos += (disp / disp_norm[:, None]) * np.minimum(disp_norm, temperature)[:, None]
        pos -= np.mean(pos, axis=0)
        temperature *= 1.0 - ((step + 1) / iterations) * 0.015

    span = np.ptp(pos, axis=0)
    span[span < 1e-9] = 1.0
    pos = pos / span
    return pos


# ---------------------------------------------------------------------------
# 포트 배정 — 각도 기반 순환 순서 보존
# ---------------------------------------------------------------------------
#
# 예전 그리디 버전(2026-09-16까지)은 커플러를 하나씩 순서대로 훑으며 각 큐빗에 "그때 비어
# 있는 후보 포트"를 배정했다. 이 방식은 같은 큐빗에 물린 커플러 여러 개가 서로 다른 순서로
# 처리되면서 포트를 "꼬아서" 배정할 수 있다 — 예: 이웃이 각도순으로 A,B,C,D인데 포트를
# A→p_top_right, B→p_bottom_right, C→p_top_left, D→p_bottom_left처럼 배정하면 A-큐빗과
# C-큐빗으로 가는 두 선이 큐빗 바로 앞에서 서로 교차한다. count_crossings()는 큐빗 "중심"
# 기준이고 정점을 공유하는 엣지 쌍은 애초에 비교 대상에서 뺀다(같은 점에서 만나는 두 선은
# "교차"가 아니라는 일반적인 그래프 교차 정의를 따른 것) — 그래서 이 포트-레벨 꼬임은
# FP의 교차-0 검증을 통과한 채로 남는다. 렌더링(utils/rendering.py)이 포트-포트 직선을
# 그리기 시작하면서 grid_25 9건, xtree_53 8건으로 실측됐다(전부 정점 공유 쌍).
#
# 고친 방식: 같은 큐빗에서 나가는 선들은 "이웃을 각도순으로 정렬한 순서"와 "배정된 포트를
# 각도순으로 정렬한 순서"가 일치하기만 하면 서로 교차할 수 없다 — 포트 4개는 전부 큐빗
# 코너(중심에서 최대 ~283um)에 몰려 있고 이웃은 보통 수백~수천 um 밖에 있어서, 포트에서
# 뻗어나가는 선의 방향은 사실상 "그 포트의 각도" 자체로 근사되기 때문이다(코너 오프셋이
# 이웃까지 거리에 비해 작을수록 이 근사가 좋아진다 — 이 저장소 규모에서는 충분히 좋다).
# 그래서 이웃을 atan2로 정렬하고, 포트 4개도 atan2로 정렬한 뒤, 순서를 그대로 옮겨 붙이면
# (필요하면 회전만 시켜서) 같은 큐빗에서 나가는 엣지끼리는 교차가 원천적으로 불가능하다.

_PORT_NAMES = ("p_top_left", "p_top_right", "p_bottom_left", "p_bottom_right")


def _port_angle(qubit: Qubit, port_name: str) -> float:
    # ports 프로퍼티는 qubit.x/y에 오프셋을 더한 절대좌표를 돌려주지만, 각도 계산에선 그
    # 오프셋(off_x, off_y)만 의미가 있고 qubit.x/y는 빼는 과정에서 상쇄된다 — 그래서 이
    # 함수는 qubit.x/y가 아직 최종값이 아니어도(_place()에서 assign_ports가 new_qubits
    # 생성보다 먼저 호출됨) 항상 안전하다.
    px, py = qubit.ports[port_name]
    return math.atan2(py - qubit.y, px - qubit.x)


def _angular_diff(a: float, b: float) -> float:
    d = abs(a - b) % (2.0 * math.pi)
    return min(d, 2.0 * math.pi - d)


# 큐빗 하나의 이웃들을 각도순 순환 순서를 보존하며 포트에 배정한다 (배정 근거는 위 섹션
# 주석 참고). 포트 4개를 각도순으로 고정해 두고, 이웃 목록(각도순 정렬, 길이 k<=4)을 그
# 순서에 "회전 오프셋" 하나로 통째로 옮겨 붙인다 — 회전 오프셋은 4가지뿐이고(포트가 4개),
# 그중 이웃-포트 각도 차이 합이 최소인 걸 고른다.
#
# 이 "연속 구간 회전" 방식이 차수(k)별로 놓치는 배정이 있는지 확인해봤다:
#   - k=4(포트 전부 사용): 두 원형 4-순열 사이에서 순서를 보존하는 전단사는 정확히 회전
#     4가지뿐이다(수학적으로 전부). 놓치는 배정 없음.
#   - k=3(포트 1개 미사용): 4개 중 1개를 빼면 남는 3개는 원형에서 항상 "연속"이다(4개짜리
#     원을 1점 제거하면 나머지는 무조건 이어져 있음) — 그래서 회전 4가지가 곧
#     C(4,3)=4가지 부분집합 전부와 일치한다. 놓치는 배정 없음.
#   - k<=2: 회전 방식은 "인접한 포트 쌍"만 만든다 — 반대쪽 코너끼리(예: p_top_right·
#     p_bottom_left) 쓰는 게 더 각도 오차가 작은 경우를 놓칠 수 있다. 다만 어떤 조합이든
#     "순서만 보존하면" 교차는 안 생기므로(교차 금지 조건은 인접을 요구하지 않는다) 이건
#     순수히 미관(각도 오차) 문제지 교차 문제가 아니다. 4가지 경우 전부를 동일한 코드
#     경로로 처리하는 단순함이 이 작은 손해보다 낫다고 판단했다.
def _assign_qubit_ports(qubit: Qubit, neighbor_ids: list[int],
                        pos: dict[int, tuple[float, float]]) -> dict[int, str]:
    if not neighbor_ids:
        return {}

    qx, qy = pos[qubit.id]
    ports_sorted = sorted(_PORT_NAMES, key=lambda p: _port_angle(qubit, p))

    def neighbor_angle(nid: int) -> float:
        return math.atan2(pos[nid][1] - qy, pos[nid][0] - qx)

    neighbors_sorted = sorted(neighbor_ids, key=neighbor_angle)
    k = len(neighbors_sorted)

    best_offset, best_cost = 0, None
    for r in range(4):
        cost = sum(
            _angular_diff(neighbor_angle(nid), _port_angle(qubit, ports_sorted[(r + i) % 4]))
            for i, nid in enumerate(neighbors_sorted)
        )
        if best_cost is None or cost < best_cost:
            best_offset, best_cost = r, cost

    return {nid: ports_sorted[(best_offset + i) % 4] for i, nid in enumerate(neighbors_sorted)}


# 큐빗별로 이웃 전체를 한 번에(_assign_qubit_ports) 배정한 뒤, 커플러의 양끝은 각자의
# 결과를 그냥 읽기만 한다. q1과 q2는 서로 "합의"할 필요가 없다 — 포트는 큐빗마다 독립된
# 4개짜리 로컬 자원이고(p_top_left 같은 이름이 q1과 q2에서 같은 이름이어도 물리적으로는
# 전혀 다른 좌표), q1의 배정은 q1 자신의 이웃 각도만으로 정해지고 q2가 자기 쪽 포트를
# 무엇으로 고르든 q1의 순환순서 보존에 전혀 영향을 주지 않는다(그 반대도 마찬가지) — 그래서
# "한쪽은 순서를 지켰는데 다른 쪽에서 깨지는" 상황 자체가 구조적으로 불가능하다. 이건 어디까지나
# "같은 큐빗을 공유하는 엣지끼리" 교차하지 않는다는 로컬 보장이라는 점은 유의: 큐빗을 공유하지
# 않는 두 커플러의 포트-포트 직선이 칩 반대편 어딘가에서 교차하는 건 이 로컬 보장 범위 밖이고,
# 그건 여전히 큐빗-중심 기준 임베딩(위쪽 count_crossings 기반 폴백 로직)이 책임진다 —
# count_port_crossings()가 그 경계를 그대로 반영해 두 값을 나눠서 센다.
def assign_ports(
    couplers: dict[tuple[int, int], Coupler],
    pos: dict[int, tuple[float, float]],
    qubits: dict[int, Qubit],
    processor_name: str = "?",
) -> dict[tuple[int, int], tuple[str, str]]:
    adj: dict[int, list[int]] = defaultdict(list)
    for q1, q2 in couplers:
        adj[q1].append(q2)
        adj[q2].append(q1)

    for qid, neighbor_ids in adj.items():
        if len(neighbor_ids) > 4:
            raise PlacementInfeasibleError(
                f"{processor_name}: qubit {qid}의 차수가 {len(neighbor_ids)}로 4를 초과합니다 "
                "(4포트 설계로는 배치 불가)."
            )

    port_of = {qid: _assign_qubit_ports(qubits[qid], neighbor_ids, pos) for qid, neighbor_ids in adj.items()}

    return {(q1, q2): (port_of[q1][q2], port_of[q2][q1]) for q1, q2 in sorted(couplers)}


def _coupler_port_line(
    qubits: dict[int, Qubit], couplers: dict[tuple[int, int], Coupler],
    port_assignment: dict[tuple[int, int], tuple[str, str]], key: tuple[int, int],
) -> tuple[tuple[float, float], tuple[float, float]]:
    c = couplers[key]
    q1, q2 = qubits[c.q1], qubits[c.q2]
    assignment = port_assignment.get(key)
    if assignment is None:
        return (q1.x, q1.y), (q2.x, q2.y)
    p1_name, p2_name = assignment
    return q1.ports[p1_name], q2.ports[p2_name]


# 실제 배선이 나가는 포트-포트 직선 기준 교차 개수. count_crossings()(큐빗 중심 기준)와
# 달리 정점(큐빗)을 공유하는 엣지 쌍도 제외하지 않는다 — assign_ports의 버그는 정확히 그
# 경우(같은 큐빗의 서로 다른 코너에서 출발하는 두 선)에서 나므로, 그 쌍을 빼면 이 검증
# 함수 자체가 원래 문제를 못 본다. 반환값을 (같은 큐빗을 공유하는 쌍의 교차, 공유하지 않는
# 쌍의 교차)로 나누는 이유는 두 값의 "보장 주체"가 다르기 때문이다 — 전자는 assign_ports
# 하나만으로 항상 0이어야 하고(수학적 불변식, 위 _assign_qubit_ports 주석의 k=3/4 분석
# 참고), 후자는 큐빗-중심 기준 임베딩이 포트 오프셋 아래에서도 여유가 있었는지에 달려 있어
# FP가 원래부터 보장하던 것의 근사치일 뿐 새 불변식이 아니다. O(E^2)라 E~150(eagle) 기준
# 만 쌍 남짓 — 매 FP 호출마다 돌려도 무시할 수준이라 별도 최적화는 하지 않는다.
def count_port_crossings(
    qubits: dict[int, Qubit], couplers: dict[tuple[int, int], Coupler],
    port_assignment: dict[tuple[int, int], tuple[str, str]],
) -> tuple[int, int]:
    keys = sorted(couplers)
    lines = {key: _coupler_port_line(qubits, couplers, port_assignment, key) for key in keys}

    shared_qubit_crossings = 0
    other_crossings = 0
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            k1, k2 = keys[i], keys[j]
            p1, p2 = lines[k1]
            p3, p4 = lines[k2]
            if not segments_cross(p1, p2, p3, p4):
                continue
            if set(k1) & set(k2):
                shared_qubit_crossings += 1
            else:
                other_crossings += 1
    return shared_qubit_crossings, other_crossings


# ---------------------------------------------------------------------------
# 배치 후보 영역 격자 스냅
# ---------------------------------------------------------------------------

# region()이 돌려주는 연속좌표 AABB를 segment_size_um 배수 격자(전역 원점 (0,0) 기준 —
# GlobalPlacement가 쓰는 것과 동일한 격자)로 바깥쪽으로 넓혀 스냅한다. floor/ceil이라
# 원래 박스를 항상 완전히 포함하며(좁아지지 않음), 스냅된 경계는 격자선과 정확히 일치해
# GP가 "박스에 완전히 포함된 셀"을 셀 때 정렬 손실이 0이 된다 — 왜 이게 필요한지는
# Floorplan._place의 coupler_regions 계산부 주석 참고.
def _snap_box_to_grid(
    box: tuple[float, float, float, float], cell: float,
) -> tuple[float, float, float, float]:
    x0, x1, y0, y1 = box
    return (
        math.floor(x0 / cell) * cell,
        math.ceil(x1 / cell) * cell,
        math.floor(y0 / cell) * cell,
        math.ceil(y1 / cell) * cell,
    )


# ---------------------------------------------------------------------------
# Floorplan
# ---------------------------------------------------------------------------

class Floorplan:
    def __init__(self, params):
        self.params = params
        self.skipped: list[tuple[str, str]] = []

    # 칩 목록을 배치한다. PlacementInfeasibleError가 난 칩은 skipped에 (이름, 사유)로
    # 기록하고 건너뛴다 — 한 칩 실패로 전체가 죽지 않는다.
    def run(self, states: list[ChipState]) -> list[ChipState]:
        out: list[ChipState] = []
        for state in states:
            try:
                out.append(self._place(state))
            except PlacementInfeasibleError as e:
                self.skipped.append((state.processor_name, str(e)))
        return out

    def _place(self, state: ChipState) -> ChipState:
        n = state.num_qubits
        edges = sort_edges(state.cmap)
        adj = build_adj(edges, n)

        # --- 평면성 판정 + 비평면 처리 -------------------------------------------
        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(edges)
        is_planar, emb = nx.check_planarity(G)
        cut_edges: list[tuple[int, int]] = []
        guarantee_graph = G
        if not is_planar:
            mode = getattr(self.params, "floorplan_planarity_mode", "crossover")
            guarantee_graph, emb, cut_edges = handle_non_planar(G, mode)
            # strict 모드는 handle_non_planar 안에서 NonPlanarError를 던지므로 여기 도달 안 함

        coupler_lengths = {edge: state.couplers[edge].l for edge in edges}

        # --- 주 경로: spectral/tree/ring 초기배치 + spring 정제 -------------------
        pos = self._primary_layout(adj, edges, state.chip_width, state.chip_height, coupler_lengths)

        # --- 보장 폴백 -----------------------------------------------------------
        # 주 경로는 교차/간격을 사후검증만 하는 필터다. 교차가 남았거나 최소 큐빗 간격을
        # 못 지킨 배치는 신뢰하지 않고 결정론적 보장 경로로 대체한다.
        guarantee_edges = [(min(u, v), max(u, v)) for u, v in guarantee_graph.edges()]
        qubit_w, qubit_h = float(self.params.qubit_width), float(self.params.qubit_height)
        min_sep = float(getattr(self.params, "min_qubit_spacing_um", max(qubit_w, qubit_h)))
        fp_margin = max(qubit_w, min(state.chip_width, state.chip_height) * 0.02)

        # min_sep은 비중첩 하한이다(params.json min_qubit_spacing_um 설명 참고) — 물리적으로
        # 안 겹친다/안다는 이분법적 기준이라 설계상 완화할 여지가 없다. 예전엔 트리거 쪽에
        # -1.0um, 최종 판정 쪽에 *0.9(10%!) 관용치가 있었는데, 후자는 실제 겹침을 최대 10%
        # 까지 통과시키는 셈이라 이 정의와 모순됐다(2026-09-16 aspen_11 조사에서 발견 — 당시엔
        # 이 두 관용치 덕에 실패가 스킵이 아니라 "간신히 통과"로 둔갑할 뻔했다). 아래
        # _SPACING_FP_EPS_UM은 그 관용치가 아니라 부동소수 반올림 오차만 흡수하는 용도다.
        _SPACING_FP_EPS_UM = 1e-6

        fallback_used = False
        if (count_crossings(pos, guarantee_edges) > 0
                or min_pairwise_distance(pos) < min_sep - _SPACING_FP_EPS_UM):
            fallback_used = True
            guaranteed_pos = guaranteed_embedding(guarantee_graph)
            # 엣지 길이 균일화(교차 불변) — Chrobak-Payne 격자 배치는 종횡비가 극단적이라
            # 스케일 전에 최대한 compact하게 만들어야 die에 담긴다.
            guaranteed_pos = impred_refine(guaranteed_pos, guarantee_edges, iterations=300)
            if count_crossings(guaranteed_pos, guarantee_edges) != 0:
                raise RuntimeError(
                    "guaranteed embedding + impred_refine left crossings before scaling "
                    "— should be impossible for a planar graph; internal bug"
                )

            pos_arr = np.array([guaranteed_pos[i] for i in range(n)], dtype=float)
            pos = self._scale_chip(pos_arr, state.chip_width, state.chip_height)
            # bbox 정규화(균일 상사변환: scale·shrink 스칼라 하나)라 정수 직선배치의 proper
            # 교차를 뒤집을 수 없다 — 여기서 >0이면 진짜 내부 버그.
            if count_crossings(pos, guarantee_edges) != 0:
                raise RuntimeError(
                    f"uniform _scale_chip introduced a crossing on a straight-line integer "
                    f"embedding ({state.processor_name}) — internal bug"
                )

            # 교차 0 보장됨. 임베딩 종횡비 탓에 최소 간격 미달이면 die 안에서 균일 확대
            # (교차 보존)로 맞춘다. 확대해도 못 맞추면 이 die로는 불가능 → 스킵.
            pos = uniform_expand_to_spacing(pos, min_sep, state.chip_width, state.chip_height, fp_margin)
            if count_crossings(pos, guarantee_edges) != 0:
                raise RuntimeError(
                    "uniform_expand_to_spacing introduced a crossing — should be impossible "
                    "for a uniform scale-preserving transform; internal bug"
                )
            achieved = min_pairwise_distance(pos)
            if achieved < min_sep - _SPACING_FP_EPS_UM:
                raise PlacementInfeasibleError(
                    f"{state.processor_name}: zero-crossing straight-line embedding does not "
                    f"fit the {state.chip_width:.0f}x{state.chip_height:.0f}um die without "
                    f"qubit-qubit overlap at {min_sep:.1f}um non-overlap floor "
                    f"(best achievable {achieved:.2f}um)"
                )

        # --- 포트 배정 -------------------------------------------------------------
        # count_port_crossings()로 실제 좌표를 검증할 것이므로, 자리표시자(x=y=0.0)가 아니라
        # 최종 좌표가 반영된 new_qubits를 먼저 만들어 assign_ports/검증 둘 다에 쓴다
        # (assign_ports 자체는 pos만 보고 포트 각도는 위치 무관이라 순서가 안 바뀌어도
        # 맞지만, 검증은 절대좌표가 맞아야 의미가 있다).
        new_qubits = {i: replace(state.qubits[i], x=pos[i][0], y=pos[i][1]) for i in range(n)}
        port_assignment = assign_ports(state.couplers, pos, new_qubits, processor_name=state.processor_name)

        shared_xing, other_xing = count_port_crossings(new_qubits, state.couplers, port_assignment)
        if shared_xing > 0:
            raise RuntimeError(
                f"{state.processor_name}: 같은 큐빗을 공유하는 커플러 쌍 사이에 포트-포트 "
                f"직선 교차가 {shared_xing}건 남았습니다 — assign_ports의 각도 기반 순환순서 "
                "배정은 차수 3~4에서 이 경우를 수학적으로 0으로 만들어야 하므로 내부 버그입니다."
            )
        if other_xing > 0:
            logging.warning(
                "[FP] %s: 큐빗을 공유하지 않는 커플러 쌍 사이에 포트-포트 직선 교차가 %d건 "
                "있습니다 — 큐빗-중심 기준 임베딩은 교차 0이지만, 포트가 중심에서 코너로 "
                "치우친 만큼 실제 배선 기준으로는 아슬아슬했던 자리입니다. assign_ports는 "
                "포트가 큐빗별 로컬 자원이라 이 경우를 고칠 권한이 없습니다 — "
                "min_qubit_spacing_um을 늘리거나 GP가 처리해야 합니다.",
                state.processor_name, other_xing,
            )

        # --- 커플러 배치 후보 영역 확정 --------------------------------------------
        # GP가 이 안에서만 세그먼트를 배치한다(core/globalplacement.py) — 즉 여기서부터는
        # Coupler.region()의 "힌트"가 아니라 GP의 하드 탐색 범위다. 포트 배정이 끝난
        # 뒤에만 계산 가능하다(region()이 port1/port2 이름을 받는다).
        #
        # region()이 돌려주는 연속좌표 박스를 그대로 저장하지 않고 segment_size_um 격자에
        # 바깥쪽으로 스냅한다(_snap_box_to_grid). 이유: region()은 면적을
        # required_wire_area(=l*meander_spacing_um) 근처에 거의 여유 없이 딱 맞춘다(ceil
        # 반올림 정도 차이) — 그런데 GP는 그 안에 segment_size_um 정사각형을 전역 격자에
        # 맞춰서만 놓을 수 있어서, 박스 경계가 격자선과 어긋나 있으면 한 칸 가까이를 정렬
        # 손실로 날린다. 2026-09-17 meander_spacing_um을 120->40으로 낮춘 뒤 실측해보니
        # (박스가 3배 작아져 원래도 여유가 빠듯했다) 이 정렬 손실 때문에 falcon 기준
        # 커플러의 46%가 "연속 면적상으론 자리가 있는데" 실패했다 — 박스 자체를 격자에
        # 스냅해 손실을 원천 제거했다(바깥쪽으로만 키우므로 GP 하드 제약이 약해지는 방향이지
        # 좁아지는 방향이 아니다: 실제 배선 여유가 늘어나는 것이지 물리적으로 잘못된 방향이
        # 아니다).
        #
        # region()이 None을 돌려줄 수 있다(2026-09-20, 가용 면적 반복 확장이 die 밖으로
        # 나가야 하거나 반복 한도 안에 못 끝난 경우 — Coupler.region() docstring 참고).
        # 그 경우 스냅을 건너뛰고 coupler_regions[key]=None으로 남긴다 — GP의 기존
        # "박스 없음" 실패 경로(core/globalplacement.py의 _place_chain)가 이 커플러 하나만
        # 실패로 세고 칩 전체는 계속 진행한다.
        cell = float(self.params.segment_size_um)
        coupler_regions: dict[tuple[int, int], tuple[float, float, float, float] | None] = {}
        n_region_failed = 0
        for key in state.couplers:
            box = state.couplers[key].region(
                new_qubits[key[0]], new_qubits[key[1]], *port_assignment[key],
                new_qubits, state.chip_width, state.chip_height,
            )
            if box is None:
                n_region_failed += 1
                coupler_regions[key] = None
            else:
                coupler_regions[key] = _snap_box_to_grid(box, cell)
        if n_region_failed:
            logging.warning(
                "[FP] %s: 커플러 %d개는 가용 면적 확장이 수렴하지 못해 박스를 못 정했습니다 "
                "(die 경계 초과 또는 반복 한도 초과) — GP에서 배치 실패로 집계됩니다.",
                state.processor_name, n_region_failed,
            )

        return replace(
            state,
            qubits=new_qubits,
            port_assignment=port_assignment,
            coupler_regions=coupler_regions,
            embedding=emb,
            cut_edges=tuple(cut_edges),
            fallback_used=fallback_used,
        )

    # 토폴로지-인식 초기 배치와 스프링 정제를 거쳐 큐빗 좌표를 계산 (주 경로 — 교차/간격
    # 미보장, _place()가 사후검증 후 필요하면 보장 폴백으로 대체한다). 01_mainref와 달리
    # 커플러 박스 기반 정제 단계(_coupler_refinement/_eliminate_overlaps/
    # _eliminate_crossings/widen_thin_coupler_spans)는 전부 뺐다 — FP는 세그먼트/커플러
    # 경계를 모른다(그건 GP의 일). edge_target만 커플러 길이(Coupler.l)로 살짝 비대칭을
    # 줘서 긴 공진기가 필요한 큐빗 쌍이 스프링 레이아웃에서 자연히 더 떨어지게 한다.
    def _primary_layout(
        self, adj: dict[int, set[int]], coupling_edges: list[tuple[int, int]],
        chip_width: float, chip_height: float,
        coupler_lengths: dict[tuple[int, int], float] | None = None,
        seed: int | None = None,
    ) -> dict[int, tuple[float, float]]:
        real_edges = set(coupling_edges)
        edges = list(coupling_edges)

        seen: set[int] = set()
        components: list[list[int]] = []
        for start in adj:
            if start in seen: continue
            stack = [start]
            seen.add(start)
            component = []
            while stack:
                node = stack.pop()
                component.append(node)
                for nbr in adj[node]:
                    if nbr not in seen:
                        seen.add(nbr)
                        stack.append(nbr)
            components.append(sorted(component))
        components.sort(key=lambda c: (-len(c), c[0]))

        if len(components) > 1:
            reps = [max(component, key=lambda q: (len(adj[q]), -q)) for component in components]
            virtual_edges = [(min(reps[i], reps[i + 1]), max(reps[i], reps[i + 1])) for i in range(len(reps) - 1)]
            edges = sorted(set(edges + virtual_edges))

        n = len(adj)
        if n <= 1 or not edges:
            raise ValueError(
                f"_primary_layout: n={n} qubits, {len(edges)} edges — trivial/degenerate "
                "coupling graph, no spectral/tree/ring layout is meaningful."
            )

        n_edges = max(len(edges), 1)
        base = math.sqrt(4.0 / n_edges)

        edge_target = {}
        if not coupler_lengths:
            edge_target = {(min(u, v), max(u, v)): base for u, v in edges}
        else:
            physical_length = [coupler_lengths[e] for e in real_edges if e in coupler_lengths]
            median_length = float(np.median(physical_length)) if physical_length else 1.0
            for u, v in edges:
                edge = (min(u, v), max(u, v))
                if edge in real_edges and edge in coupler_lengths:
                    ratio = coupler_lengths[edge] / max(median_length, 1e-9)
                    edge_target[edge] = base * float(np.clip(0.85 + 0.28 * ratio, 0.9, 1.75))
                else:
                    edge_target[edge] = base

        bst_pos = None
        bst_score = None
        base_seed = 0 if seed is None else seed
        base_candidates = int(getattr(self.params, "floorplan_candidate_count", 48))

        n_scale = max(n, 1)
        effort_scale = math.sqrt(n_scale / 27.0)
        candidate_cnt = max(8, int(base_candidates * effort_scale))
        spring_iters = max(200, int(420 * effort_scale))

        octagon_rings = _detect_octagon_rings(list(real_edges), n)
        if not _ring_coverage_sufficient(octagon_rings, n):
            octagon_rings = []
        if len(components) > 1:
            octagon_rings = []
        tree_root = None if octagon_rings else _detect_tree_root(list(real_edges), n, components)
        ring_aspect = float(getattr(self.params, "octagon_ring_aspect", 1.6))
        structured_spring_iters = spring_iters
        structured_max_temperature = 0.35
        if octagon_rings or tree_root is not None:
            candidate_cnt = max(8, int(candidate_cnt * 0.5))
            structured_spring_iters = max(80, int(spring_iters * 0.4))
            structured_max_temperature = 0.12

        offset = 0
        for offset in range(max(candidate_cnt, 1)):
            rng = np.random.default_rng(base_seed + offset)
            if octagon_rings:
                pos = _ring_aware_initial_coordinate(adj, edges, octagon_rings, edge_target, base, rng,
                                                      aspect=ring_aspect)
                pos = _spring_refine_coordinate(pos, edges, edge_target=edge_target,
                                                 iterations=structured_spring_iters,
                                                 max_temperature=structured_max_temperature)
            elif tree_root is not None:
                pos = _tree_aware_initial_coordinate(edges, tree_root, edge_target, base, rng, n)
                pos = _spring_refine_coordinate(pos, edges, edge_target=edge_target,
                                                 iterations=structured_spring_iters,
                                                 max_temperature=structured_max_temperature)
            else:
                pos = _spectral_initial_coordinate(adj, edges, rng)
                pos = _spring_refine_coordinate(pos, edges, edge_target=edge_target, iterations=spring_iters)

            crossing_cnt = _cnt_edge_crossing(pos, sorted(real_edges))
            bbox_area = float(np.prod(np.maximum(np.ptp(pos, axis=0), 1e-9)))
            score = crossing_cnt * 1_000_000.0 + bbox_area
            if bst_score is None or score < bst_score:
                bst_score, bst_pos = score, pos.copy()
                if crossing_cnt == 0:
                    break

        return self._scale_chip(bst_pos, chip_width, chip_height)

    # bbox 정규화만 한다(01_mainref의 커플러 target_span 기반 분기는 제외) — FP는 커플러
    # 스펙/세그먼트를 다루지 않으므로 그 분기가 필요로 하는 정보 자체가 없다. 배치를
    # 중심으로 모은 뒤, die 안(마진 제외)에 딱 맞도록 균일 확대/축소한다.
    def _scale_chip(self, pos: np.ndarray, chip_width: float, chip_height: float) -> dict[int, tuple[float, float]]:
        qubit_width, qubit_height = float(self.params.qubit_width), float(self.params.qubit_height)
        margin_x = max(qubit_width, chip_width * 0.02)
        margin_y = max(qubit_height, chip_height * 0.02)

        centered = pos - np.mean(pos, axis=0)
        half_w = max(chip_width / 2.0 - margin_x, 1.0)
        half_h = max(chip_height / 2.0 - margin_y, 1.0)
        extent = np.max(np.abs(centered), axis=0)
        extent[extent < 1e-9] = 1.0
        shrink = min(half_w / extent[0], half_h / extent[1])
        centered = centered * shrink

        cx, cy = chip_width / 2.0, chip_height / 2.0
        x = cx + centered[:, 0]
        y = cy + centered[:, 1]
        return {i: (float(x[i]), float(y[i])) for i in range(len(pos))}