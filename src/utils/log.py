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
