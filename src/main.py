import logging
import time
import sys
import os
import csv
import numpy as np
from dataclasses import replace
from pathlib import Path

# ./core
from core.state import ChipState
from core.floorplan import Floorplan
from core.globalplacement import GlobalPlacement
from core.legalization import Legalization
from core.detailedplacement import DetailedPlacement
from core.router import Router
from core.fabrication import FabricationMapping

# ./eval
from func.evaluation import Evaluation

# ./utils
from utils.log import setup as setup_logging, log_chip_states
from utils.params import Params
from utils.parsing import Parser
from utils.rendering import Rendering, rendering_qiskit_metal

# Default
FP      = "output/0_FP"
GP      = "output/1_GP"
LG      = "output/2_LG"
DP      = "output/3_DP"
RT      = "output/4_RT"
FM      = "output/6_FABRICATION"

QM_FP   = "output/5_QISKIT_METAL/0_FP"
QM_GP   = "output/5_QISKIT_METAL/1_GP"
QM_LG   = "output/5_QISKIT_METAL/2_LG"
QM_DP   = "output/5_QISKIT_METAL/3_DP"
QM_RT   = "output/5_QISKIT_METAL/4_RT"

# SAVE_LOG = True
SAVE_LOG = False

def main(params):
    """ Preprocessing """
    stages = ("FP", "GP", "LG", "DP", "RT", "FM")
    runtimes: dict[str, dict[str, float]] = {}

    def run_stage(stage_name, stage, states):
        """Run one stage per benchmark so its runtime can be reported separately."""
        output = []
        for state in states:
            tt = time.perf_counter()
            result = stage.run([state])
            runtimes.setdefault(state.processor_name, {})[stage_name] = (
                time.perf_counter() - tt
            )
            output.extend(result)
        return output

    parser = Parser()
    bench = parser.load_configs(params.processors)
    logging.info("Processors: %s" % [c['processor'] for c in bench])
    init_state: list[ChipState] = [ChipState.processor_config(c, params) for c in bench]
    # log_chip_states("Preprocessing", init_state)


    """ Floorplanning """
    fp = Floorplan(params)
    fp_state = run_stage("FP", fp, init_state)
    # log_chip_states("FP", fp_state)
    for name, reason in fp.skipped:
        logging.warning("[FP] skipped %s: %s", name, reason)


    """ Global Placement"""
### 전반적인 배치 최적화 --> Planarity, Crosstalk, Wirelength, Congestion
    gp = GlobalPlacement(params)
    gp_state = run_stage("GP", gp, fp_state)
    # log_chip_states("GP", gp_state)


    """ Legalization """
### Overlap 해소
    lg = Legalization(params)
    lg_state = run_stage("LG", lg, gp_state)
    # log_chip_states("LG", lg_state)


    """ Detailed Placement """
### Hotspot이 가장 심한 위치 파악 후, 미세 조정
    dp = DetailedPlacement(params)
    dp_state = run_stage("DP", dp, lg_state)
    # log_chip_states("DP", dp_state)


    """ Routing """
### Routing, Finish단계
    rt = Router(params)
    rt_state = run_stage("RT", rt, dp_state)


    """ Fabrication Mapping """
### Qiskit-Metal/제조 footprint 기준 DRC 진단
    fm = FabricationMapping(params, FM)
    fm_state = run_stage("FM", fm, rt_state)
    _write_failure_report(FM, gp, rt, rt_state, fm)


    """ Basic Rendering """
    Rendering(FP, fp_state)
    Rendering(GP, gp_state)
    Rendering(LG, lg_state)
    Rendering(DP, dp_state)
    Rendering(RT, fm_state)

    """ Qiskit-Metal Rendering """
    # rendering_qiskit_metal(QM_FP, fp_state)
    # rendering_qiskit_metal(QM_GP, gp_state)
    # rendering_qiskit_metal(QM_LG, lg_state)
    # rendering_qiskit_metal(QM_DP, dp_state)
    rendering_qiskit_metal(QM_RT, fm_state)


    """ Evaluation """
    eval = Evaluation(params)

    eval_step = [
        ("FloorPlanning",       fp_state),
        ("GlobalPlacement",     gp_state),
        ("Legalization",        lg_state),
        ("DetailedPlacement",   dp_state),
        ("Routing",             fm_state)
    ]

    evaluation_results = {}

    for stage, states in eval_step:
        evaluation_results[stage] = {}

        for state in states:
            metric = eval.compute(state)
            logging.info("[Evaluation][%s][%s]", stage, state.processor_name)
            evaluation_results[stage][state.processor_name] = metric

    # Evaluation 결과 출력
    benchmark_names = [c["processor"] for c in bench]

    eval.print_evaluation_table(evaluation_results, benchmark_names)
##    eval.save_csv(evaluation_results, benchmark_names)

    print("\nRuntime by stage (sec)")
    header = ["Stage", *benchmark_names, "Total"]
    stage_width = 6
    benchmark_width = 8
    total_width = 8
    widths = [stage_width, *([benchmark_width] * len(benchmark_names)), total_width]

    print(" | ".join(
        f"{column:>{width}}"
        for column, width in zip(header, widths)
    ))

    print("-+-".join("-" * width for width in widths))

    for stage in stages:
        values = [
            runtimes.get(benchmark, {}).get(stage)
            for benchmark in benchmark_names
        ]

        total = sum(value for value in values if value is not None)

        formatted = [
            f"{value:{width}.3f}" if value is not None
            else f"{'-':>{width}}"
            for value, width in zip(values, widths[1:-1])
        ]

        print(" | ".join([
            f"{stage:>{stage_width}}",
            *formatted,
            f"{total:{total_width}.3f}",
        ]))




def _write_failure_report(
    output_dir, gp: GlobalPlacement, rt: Router, states: list[ChipState],
    fm: FabricationMapping | None = None,
):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, "routing_failure_report.csv")
    state_by_name = {state.processor_name: state for state in states}
    rows = []

    for processor, failures in sorted(gp.placement_failures.items()):
        state = state_by_name.get(processor)
        for key, reason in failures:
            coupler = state.couplers[key] if state is not None else None
            rows.append(_failure_row(
                processor, "GP", key, reason, coupler,
                gp.failure_stats.get(processor, {}).get(key, {}),
            ))

    for processor, failures in sorted(rt.route_failures.items()):
        state = state_by_name.get(processor)
        for key, reason in failures:
            coupler = state.couplers[key] if state is not None else None
            rows.append(_failure_row(processor, "RT", key, reason, coupler))

    if fm is not None:
        for processor, rejected in sorted(fm.rejected_routes.items()):
            state = state_by_name.get(processor)
            for key in sorted(rejected):
                coupler = state.couplers[key] if state is not None else None
                rows.append(_failure_row(processor, "FM", key, "fabrication_rule_rejected", coupler))

    fieldnames = [
        "processor",
        "stage",
        "q1",
        "q2",
        "reason",
        "target_length_um",
        "route_length_um",
        "length_error_pct",
        "segments",
        "required_segments",
        "region_cells",
        "free_cells_at_failure",
        "occupied_cells_at_failure",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logging.info("[FailureReport] 저장: %s (%d rows)", path, len(rows))


def _failure_row(
    processor: str, stage: str, key: tuple[int, int], reason: str, coupler,
    stats: dict | None = None,
) -> dict[str, object]:
    stats = stats or {}
    route_length = coupler.route_length if coupler is not None else 0.0
    target_length = coupler.l if coupler is not None else 0.0
    err_pct = (
        (route_length - target_length) / target_length * 100.0
        if target_length > 0.0 and route_length > 0.0 else ""
    )
    return {
        "processor": processor,
        "stage": stage,
        "q1": key[0],
        "q2": key[1],
        "reason": reason,
        "target_length_um": f"{target_length:.6f}" if target_length else "",
        "route_length_um": f"{route_length:.6f}" if route_length else "",
        "length_error_pct": f"{err_pct:.6f}" if err_pct != "" else "",
        "segments": len(coupler.segments) if coupler is not None else "",
        "required_segments": stats.get("required_segments", ""),
        "region_cells": stats.get("region_cells", ""),
        "free_cells_at_failure": stats.get("free_cells_at_failure", ""),
        "occupied_cells_at_failure": stats.get("occupied_cells_at_failure", ""),
    }


if __name__ == "__main__":
    params = Params()

    if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
        params.printHelp()
        exit()

    if len(sys.argv) == 2:
        params.load(sys.argv[1])
    elif len(sys.argv) > 2:
        logging.error("Input parameters is required in json format")
        params.printHelp()
        exit()

    setup_logging("PlanarQ", SAVE_LOG)

    os.environ["OMP_NUM_THREADS"] = "%d" % (params.num_threads)
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

    input_dir = Path(__file__).resolve().parent.parent / "input"
    input_dir.mkdir(exist_ok=True)

    main(params)
