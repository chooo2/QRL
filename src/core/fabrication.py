"""Fabrication mapping checks for routed couplers.

This stage maps routed centerlines to the CPW footprint that Qiskit-Metal would
draw around them: trace metal plus clearance cut.  It is intentionally
diagnostic first: the state is returned unchanged, while detailed violations are
written as CSV so LG/RT can be tuned against fabrication-level geometry instead
of only centerlines.
"""
from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, replace
from pathlib import Path

from shapely.geometry import LineString, box
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union

from core.state import ChipState, Coupler

_EPS_AREA_UM2 = 1e-6
_EPS_LEN_UM = 1e-6


def _leg_cut_geometries(
    waypoints: list[tuple[float, float]], half_width_um: float,
) -> list[BaseGeometry]:
    out: list[BaseGeometry] = []
    for p0, p1 in zip(waypoints, waypoints[1:]):
        if ((p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2) ** 0.5 <= _EPS_LEN_UM:
            continue
        out.append(
            LineString([p0, p1]).buffer(half_width_um, cap_style=2, join_style=2)
        )
    return out


@dataclass(frozen=True)
class _Footprint:
    key: tuple[int, int]
    coupler: Coupler
    line: LineString
    trace: BaseGeometry
    cut: BaseGeometry


class FabricationMapping:
    def __init__(self, params, output_dir: str):
        self.params = params
        self.output_dir = Path(output_dir)
        self.trace_width_um = float(getattr(params, "cpw_trace_width_um", 10.0))
        self.trace_gap_um = float(getattr(params, "cpw_trace_gap_um", 6.0))
        self.cut_half_width_um = self.trace_width_um / 2.0 + self.trace_gap_um
        self.fillet_radius_um = float(getattr(params, "cpw_fillet_radius_um", 30.0))
        self.min_coupler_spacing_um = float(
            getattr(params, "cpw_min_coupler_spacing_um", 200.0)
        )
        self.reject_spacing_violations = bool(
            getattr(params, "fm_reject_spacing_violations", True)
        )
        self.reports: dict[str, dict[str, int]] = {}
        self.rejected_routes: dict[str, set[tuple[int, int]]] = {}
        self._detail_rows: list[dict[str, object]] = []
        self._summary_rows: list[dict[str, object]] = []

    def run(self, states: list[ChipState]) -> list[ChipState]:
        self.output_dir.mkdir(parents=True, exist_ok=True)

        out: list[ChipState] = []
        for state in states:
            summary, rejected = self._check_state(state, self._detail_rows)
            mapped_state = self._drop_rejected_routes(state, rejected)
            self.reports[state.processor_name] = summary
            self.rejected_routes[state.processor_name] = rejected
            self._summary_rows.append({"processor": state.processor_name, **summary})
            logging.info(
                "[FM] %s: routed=%d cc_overlap=%d qc_overlap=%d region_escape=%d "
                "die_escape=%d short_leg=%d cc_spacing=%d rejected=%d",
                state.processor_name,
                summary["routed"],
                summary["cc_overlap"],
                summary["qc_overlap"],
                summary["region_escape"],
                summary["die_escape"],
                summary["short_leg"],
                summary["cc_spacing_violation"],
                summary["rejected_routes"],
            )
            out.append(mapped_state)

        self._write_csv(self.output_dir / "fabrication_report.csv", self._detail_rows)
        self._write_csv(self.output_dir / "fabrication_summary.csv", self._summary_rows)
        return out

    def _check_state(
        self, state: ChipState, detail_rows: list[dict[str, object]]
    ) -> tuple[dict[str, int | float], set[tuple[int, int]]]:
        fps = self._footprints(state)
        summary = {
            "couplers": len(state.couplers),
            "routed": len(fps),
            "unrouted": len(state.couplers) - len(fps),
            "safe_routed": 0,
            "cc_overlap": 0,
            "qc_overlap": 0,
            "own_qubit_launch_overlap": 0,
            "region_escape": 0,
            "die_escape": 0,
            "short_leg": 0,
            "cc_spacing_violation": 0,
            "self_spacing_violation": 0,
            "min_cc_spacing_um": "",
            "min_self_spacing_um": "",
            "rejected_routes": 0,
        }
        rejected: set[tuple[int, int]] = set()
        min_spacing: float | None = None
        min_self_spacing: float | None = None

        die = box(0.0, 0.0, state.chip_width, state.chip_height)
        qubit_boxes = {
            qid: box(q.x - q.w / 2.0, q.y - q.h / 2.0, q.x + q.w / 2.0, q.y + q.h / 2.0)
            for qid, q in state.qubits.items()
        }

        for fp in fps:
            key_s = self._key_s(fp.key)
            outside_die = fp.cut.difference(die)
            if outside_die.area > _EPS_AREA_UM2:
                summary["die_escape"] += 1
                rejected.add(fp.key)
                detail_rows.append(self._row(state, key_s, "die_escape", area=outside_die.area))

            region = state.coupler_regions.get(fp.key)
            if region is not None:
                rx0, rx1, ry0, ry1 = region
                allowed_region = box(rx0, ry0, rx1, ry1).buffer(
                    self.cut_half_width_um, cap_style=2, join_style=2
                )
                own_qubit_geoms = [
                    qubit_boxes[qid] for qid in fp.key if qid in qubit_boxes
                ]
                if own_qubit_geoms:
                    allowed_region = unary_union([allowed_region, *own_qubit_geoms])
                outside_region = fp.cut.difference(allowed_region)
                if outside_region.area > _EPS_AREA_UM2:
                    summary["region_escape"] += 1
                    rejected.add(fp.key)
                    detail_rows.append(
                        self._row(state, key_s, "region_escape", area=outside_region.area)
                    )

            for leg_idx, (p0, p1) in enumerate(zip(fp.coupler.waypoints, fp.coupler.waypoints[1:])):
                leg_len = ((p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2) ** 0.5
                if leg_len + _EPS_LEN_UM < 2.0 * self.fillet_radius_um:
                    summary["short_leg"] += 1
                    rejected.add(fp.key)
                    detail_rows.append(
                        self._row(state, key_s, "short_leg", index=leg_idx, length_um=leg_len)
                    )

            leg_cuts = _leg_cut_geometries(fp.coupler.waypoints, self.cut_half_width_um)
            for i, leg_a in enumerate(leg_cuts):
                for j in range(i + 2, len(leg_cuts)):
                    if i == 0 and j == len(leg_cuts) - 1:
                        continue
                    leg_b = leg_cuts[j]
                    overlap = leg_a.intersection(leg_b)
                    if overlap.area > _EPS_AREA_UM2:
                        spacing = 0.0
                    else:
                        spacing = leg_a.distance(leg_b)
                    if min_self_spacing is None or spacing < min_self_spacing:
                        min_self_spacing = spacing
                    if spacing + _EPS_LEN_UM < self.min_coupler_spacing_um:
                        summary["self_spacing_violation"] += 1
                        rejected.add(fp.key)
                        detail_rows.append(
                            self._row(
                                state,
                                key_s,
                                "self_spacing_violation",
                                length_um=spacing,
                                index=i,
                                other_index=j,
                            )
                        )
                        break
                if fp.key in rejected:
                    break

            for qid, qgeom in qubit_boxes.items():
                overlap = fp.cut.intersection(qgeom)
                if overlap.area <= _EPS_AREA_UM2:
                    continue
                if fp.coupler.is_own_qubit(qid):
                    summary["own_qubit_launch_overlap"] += 1
                    violation = "own_qubit_launch_overlap"
                else:
                    summary["qc_overlap"] += 1
                    rejected.add(fp.key)
                    violation = "qc_overlap"
                detail_rows.append(
                    self._row(state, key_s, violation, qubit=qid, area=overlap.area)
                )

        for i, a in enumerate(fps):
            for b in fps[i + 1:]:
                overlap = a.cut.intersection(b.cut)
                shared_qubits = set(a.key) & set(b.key)
                a_space = a.cut
                b_space = b.cut
                if shared_qubits:
                    shared_geoms = [
                        qubit_boxes[qid] for qid in shared_qubits if qid in qubit_boxes
                    ]
                    if shared_geoms:
                        shared_union = unary_union(shared_geoms)
                        overlap = overlap.difference(shared_union)
                        a_space = a_space.difference(shared_union)
                        b_space = b_space.difference(shared_union)
                if overlap.area <= _EPS_AREA_UM2:
                    spacing = a_space.distance(b_space)
                    if min_spacing is None or spacing < min_spacing:
                        min_spacing = spacing
                    if spacing + _EPS_LEN_UM < self.min_coupler_spacing_um:
                        summary["cc_spacing_violation"] += 1
                        if self.reject_spacing_violations:
                            rejected.add(b.key)
                        detail_rows.append(
                            self._row(
                                state,
                                f"{self._key_s(a.key)}|{self._key_s(b.key)}",
                                "cc_spacing_violation",
                                length_um=spacing,
                                shared_qubit=bool(shared_qubits),
                            )
                        )
                    continue
                summary["cc_overlap"] += 1
                rejected.add(b.key)
                detail_rows.append(
                    self._row(
                        state,
                        f"{self._key_s(a.key)}|{self._key_s(b.key)}",
                        "cc_overlap",
                        area=overlap.area,
                        shared_qubit=bool(shared_qubits),
                    )
                )

        summary["min_cc_spacing_um"] = "" if min_spacing is None else round(min_spacing, 6)
        summary["min_self_spacing_um"] = (
            "" if min_self_spacing is None else round(min_self_spacing, 6)
        )
        summary["rejected_routes"] = len(rejected)
        summary["safe_routed"] = len(fps) - len(rejected)
        return summary, rejected

    @staticmethod
    def _drop_rejected_routes(
        state: ChipState, rejected: set[tuple[int, int]]
    ) -> ChipState:
        if not rejected:
            return state
        new_couplers = {
            key: replace(coupler, waypoints=[])
            if key in rejected else coupler
            for key, coupler in state.couplers.items()
        }
        return replace(state, couplers=new_couplers)

    def _footprints(self, state: ChipState) -> list[_Footprint]:
        half_cut = self.trace_width_um / 2.0 + self.trace_gap_um
        half_trace = self.trace_width_um / 2.0
        out: list[_Footprint] = []
        for key in sorted(state.couplers):
            coupler = state.couplers[key]
            if len(coupler.waypoints) < 2:
                continue
            line = LineString(coupler.waypoints)
            if line.length <= _EPS_LEN_UM:
                continue
            out.append(
                _Footprint(
                    key=key,
                    coupler=coupler,
                    line=line,
                    trace=line.buffer(half_trace, cap_style=2, join_style=2),
                    cut=line.buffer(half_cut, cap_style=2, join_style=2),
                )
            )
        return out

    @staticmethod
    def _key_s(key: tuple[int, int]) -> str:
        return f"{key[0]}-{key[1]}"

    @staticmethod
    def _row(
        state: ChipState,
        coupler: str,
        violation: str,
        *,
        area: float | None = None,
        length_um: float | None = None,
        index: int | None = None,
        other_index: int | None = None,
        qubit: int | None = None,
        shared_qubit: bool | None = None,
    ) -> dict[str, object]:
        return {
            "processor": state.processor_name,
            "coupler": coupler,
            "violation": violation,
            "area_um2": "" if area is None else f"{area:.6f}",
            "length_um": "" if length_um is None else f"{length_um:.6f}",
            "index": "" if index is None else index,
            "other_index": "" if other_index is None else other_index,
            "qubit": "" if qubit is None else qubit,
            "shared_qubit": "" if shared_qubit is None else int(shared_qubit),
        }

    @staticmethod
    def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
        if rows:
            fields = list(rows[0].keys())
        else:
            fields = [
                "processor",
                "coupler",
                "violation",
                "area_um2",
                "length_um",
                "index",
                "other_index",
                "qubit",
                "shared_qubit",
            ]
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
