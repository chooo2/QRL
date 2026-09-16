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
    log_chip_states("Preprocessing", init_state)


    """ Floorplanning """
    fp_state: list[ChipState] = []
    tt = time.perf_counter()
    Floorplan(init_state, fp_state)
    run["FP"] = time.perf_counter() - tt
    log_chip_states("FP", fp_state)
    qm_fp_state: list[ChipState] = []


    """ Global Placement"""
### 전반적인 배치 최적화 --> Planarity, Crosstalk, Wirelength, Congestion
    gp_state: list[ChipState] = []
    tt = time.perf_counter()
    GlobalPlacement(fp_state, gp_state)
    run["GP"] = time.perf_counter() - tt
    log_chip_states("GP", gp_state)
    qm_gp_state: list[ChipState] = [] # 나중에 qiskit-metal 호환으로 convert시켜야함


    """ Legalization """
### Overlap 해소
    lg_state: list[ChipState] = []
    tt = time.perf_counter()
    Legalization(gp_state, lg_state)
    run["LG"] = time.perf_counter() - tt
    log_chip_states("LG", lg_state)
    qm_lg_state: list[ChipState] = [] # 나중에 qiskit-metal 호환으로 convert시켜야함


    """ Detailed Placement """
### Hotspot이 가장 심한 위치 파악 후, 미세 조정
    dp_state: list[ChipState] = []
    tt = time.perf_counter()
    DetailedPlacement(lg_state, dp_state)
    run["DP"] = time.perf_counter() - tt
    log_chip_states("DP", dp_state)
    qm_dp_state: list[ChipState] = [] # 나중에 qiskit-metal 호환으로 convert시켜야함


    """ Routing """
### Routing, Finish단계
    rt_state: list[ChipState] = []
    tt = time.perf_counter()
    Router(dp_state, rt_state)
    run["Finish"] = time.perf_counter() - tt
    qm_rt_state: list[ChipState] = [] # 나중에 qiskit-metal 호환으로 convert시켜야함


    """ Basic Rendering """
    Rendering(FP, fp_state)
    Rendering(GP, gp_state)
    Rendering(LG, lg_state)
    Rendering(DP, dp_state)
    Rendering(RT, rt_state)

    """ Qiskit-Metal Rendering """
    rendering_qiskit_metal(QM_FP, qm_fp_state)
    rendering_qiskit_metal(QM_GP, qm_gp_state)
    rendering_qiskit_metal(QM_LG, qm_lg_state)
    rendering_qiskit_metal(QM_DP, qm_dp_state)
    rendering_qiskit_metal(QM_RT, qm_rt_state)


    """ Evaluation """
    eval = Evaluation(params)
    eval_step = [
        ("FloorPlanning",       fp_state),
        ("GlobalPlacement",     gp_state),
        ("Legalization",        lg_state),
        ("DetailedPlacement",   dp_state),
    ]

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
