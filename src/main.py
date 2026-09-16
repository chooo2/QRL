import logging
import time
import sys
import os
import torch
import numpy as np
from dataclasses import replace
from pathlib import Path

# ./core
from core.state import ChipState
from core.globalplacement import GlobalPlacement
from core.legalization import Legalization
from core.detailedplacement import DetailedPlacement

# ./eval
from func.evaluation import Evaluation

# ./utils
from utils.log import setup as setup_logging
from utils.params import Params
from utils.parsing import Parser
from utils.rendering import Rendering, rendering_qiskit_metal

# Default
GP      = "output/0_GP"
LG      = "output/1_LG"
DP      = "output/2_DP"
FINISH  = "output/3_FINISH"
QM      = "output/4_QISKIT_METAL"

SAVE_LOG = False

def main(params, device):
    """ Preprocessing """
    run: dict[str, float] = {}
    parser = Parser()
    bench = parser.load_configs(params.processors)
    logging.info("Processors: %s" % [c['topology'] for c in bench])


    """ Global Placement"""
### 전반적인 배치 최적화 --> Planarity, Crosstalk, Wirelength, Congestion
    gp_state: list[ChipState] = []

    tt = time.perf_counter()
    ## Global placement algorithm
    run["GP"] = time.perf_counter() - tt

    qm_gp_state: list[ChipState] = []



    """ Legalization """
### Overlap 해소
    lg_state: list[ChipState] = []

    tt = time.perf_counter()
    ## Legalization algorithm
    run["LG"] = time.perf_counter() - tt

    qm_lg_state: list[ChipState] = []


    """ Detailed Placement """
### Hotspot이 가장 심한 위치 파악 후, 미세 조정
    dp_state: list[ChipState] = []

    tt = time.perf_counter()
    ## Detailed placement algorithm
    run["DP"] = time.perf_counter() - tt

    qm_dp_state: list[ChipState] = []



    """ Routing """
### Routing, Finish단계
    rt_state: list[ChipState] = []

    tt = time.perf_counter()
    ## Routing algorithm
    run["Finish"] = time.perf_counter() - tt

    qm_rt_state: list[ChipState] = []

    """ Basic Rendering """
    Rendering(GP, gp_state)
    Rendering(LG, lg_state)
    Rendering(DP, dp_state)
    Rendering(FINISH, rt_state)

    """ Qiskit-Metal Rendering """
    ## only rendering


    """ Evaluation """
    eval = Evaluation(params)
    eval_step = [
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

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    input_dir = Path(__file__).resolve().parent.parent / "input"
    input_dir.mkdir(exist_ok=True)

    main(params, device)
