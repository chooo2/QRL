import logging
import sys
from datetime import datetime
from pathlib import Path

_LOG_DIR = Path(__file__).resolve().parents[2] / "log"
_FMT = "[%(levelname)-4s] %(message)s"

def setup(name, save_log):
    formatter = logging.Formatter(_FMT)

    root = logging.getLogger()
    root.name = name
    root.setLevel(logging.INFO)
    root.handlers.clear()

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(formatter)
    root.addHandler(sh)

    if not save_log:
        return None

    _LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = _LOG_DIR / f"{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    logging.info("Log file: %s", log_file)
    return log_file

## 각 단계 결과(ChipState) 확인부
def log_chip_states(stage, states):
    for s in states:
        logging.info("[%s][%s] num_qubits=%d chip=(%s, %s) qubits=%d couplers=%d cmap=%d",
                      stage, s.processor_name, s.num_qubits, s.chip_width, s.chip_height,
                      len(s.qubits), len(s.couplers), len(s.cmap))

        logging.info("Qubit")
        for qid in sorted(s.qubits): logging.info(" %s", s.qubits[qid])

        logging.info("Coupler")
        for cid in sorted(s.couplers):
            c = s.couplers[cid]
            logging.info(" %s l=%.2fum", c, c.l)