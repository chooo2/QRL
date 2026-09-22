import logging
import time
import sys
import os
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

QM_FP   = "output/5_QISKIT_METAL/0_FP"
QM_GP   = "output/5_QISKIT_METAL/1_GP"
QM_LG   = "output/5_QISKIT_METAL/2_LG"
QM_DP   = "output/5_QISKIT_METAL/3_DP"
QM_RT   = "output/5_QISKIT_METAL/4_RT"

# SAVE_LOG = True
SAVE_LOG = False

def main(params):
    """ Preprocessing """
    run: dict[str, float] = {}
    parser = Parser()
    bench = parser.load_configs(params.processors)
    logging.info("Processors: %s" % [c['processor'] for c in bench])
    init_state: list[ChipState] = [ChipState.processor_config(c, params) for c in bench]
    # log_chip_states("Preprocessing", init_state)


    """ Floorplanning """
    tt = time.perf_counter()
    fp = Floorplan(params)
    fp_state = fp.run(init_state)
    run["FP"] = time.perf_counter() - tt
    # log_chip_states("FP", fp_state)
    for name, reason in fp.skipped:
        logging.warning("[FP] skipped %s: %s", name, reason)


    """ Global Placement"""
### 전반적인 배치 최적화 --> Planarity, Crosstalk, Wirelength, Congestion
    tt = time.perf_counter()
    gp = GlobalPlacement(params)
    gp_state = gp.run(fp_state)
    run["GP"] = time.perf_counter() - tt
    # log_chip_states("GP", gp_state)


    """ Legalization """
### Overlap 해소
    tt = time.perf_counter()
    lg = Legalization(params)
    lg_state = lg.run(gp_state)
    run["LG"] = time.perf_counter() - tt
    # log_chip_states("LG", lg_state)


    """ Detailed Placement """
### Hotspot이 가장 심한 위치 파악 후, 미세 조정
    tt = time.perf_counter()
    dp = DetailedPlacement(params)
    dp_state = dp.run(lg_state)
    run["DP"] = time.perf_counter() - tt
    # log_chip_states("DP", dp_state)


    """ Routing """
### Routing, Finish단계
    tt = time.perf_counter()
    rt = Router(params)
    rt_state = rt.run(dp_state)
    run["Finish"] = time.perf_counter() - tt


    """ Basic Rendering """
    Rendering(FP, fp_state)
    Rendering(GP, gp_state)
    Rendering(LG, lg_state)
    Rendering(DP, dp_state)
    Rendering(RT, rt_state)

    """ Qiskit-Metal Rendering """
    rendering_qiskit_metal(QM_FP, fp_state)
    rendering_qiskit_metal(QM_GP, gp_state)
    rendering_qiskit_metal(QM_LG, lg_state)
    rendering_qiskit_metal(QM_DP, dp_state)
    rendering_qiskit_metal(QM_RT, rt_state)


    """ Evaluation """
    eval = Evaluation(params)
    eval_step = [
        ("FloorPlanning",       fp_state),
        ("GlobalPlacement",     gp_state),
        ("Legalization",        lg_state),
        ("DetailedPlacement",   dp_state),
    ]
    for stage, states in eval_step:
        for state in states:
            metric = eval.compute(state)
            logging.info("[Evaluation][%s][%s]", stage, state.processor_name)
            metric.log_result()

    run["total"] = sum(run.values())
    for p, r in run.items(): logging.info("%-20s : %.3f sec", p, r)






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
