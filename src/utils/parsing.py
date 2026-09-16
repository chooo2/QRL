import json
import logging
from pathlib import Path

_log = logging.getLogger(__name__)

# input/ 디렉터리의 프로세서 JSON 설정을 읽어오는 파서
class Parser:
    # 지정 경로의 모든 JSON 파일을 읽어 processor config 목록으로 반환
    def load_json(self, file_path="input"):
        processor_config = []
        base_path = Path(__file__).resolve().parents[2]
        target_path = Path(file_path)
        if not target_path.is_absolute():
            target_path = (base_path / target_path).resolve()

        input_roots = [target_path]
        if target_path == (base_path / "input").resolve():
            input_roots.extend(
                p for p in ((base_path / "input2").resolve(), (base_path / "input3").resolve())
                if p.exists()
            )

        seen_processors = set()
        for file in (file for root in input_roots for file in sorted(root.rglob("*.json"))):
            try:
                with open(file, 'r') as f:
                    data = json.load(f)

                    processor_name = data.get("processor")
                    if processor_name in seen_processors:
                        continue
                    seen_processors.add(processor_name)

                    raw_coupling_map = data.get("coupling_graph")
                    topology_family = data.get("topology") if isinstance(data.get("topology"), str) else None
                    if raw_coupling_map is None:
                        raw_topology = data.get("topology")
                        if isinstance(raw_topology, list):
                            raw_coupling_map = raw_topology
                        else:
                            raw_coupling_map = []

                    re_coupling_map = list(set((min(u, v), max(u, v)) for u, v in raw_coupling_map))
                    re_coupling_map.sort()

                    chip_width_um = data.get("chip_width_um")
                    chip_height_um = data.get("chip_height_um")
                    chip_length_um = data.get("chip_length_um")
                    if chip_length_um is not None:
                        chip_width_um = chip_width_um if chip_width_um is not None else chip_length_um
                        chip_height_um = chip_height_um if chip_height_um is not None else chip_length_um

                    config = {
                        "topology": processor_name,
                        "topology_family": topology_family,
                        "num_qubits": data.get("num_qubits"),
                        "coupling_map": re_coupling_map,
                        "freq_ghz": data.get("frequency_ghz") or data.get("node_freq_ghz"),
                        "anharmonicity_ghz": data.get("anharmonicity_ghz"),
                        "edge_freq_ghz": data.get("edge_frequency_ghz") or data.get("edge_freq_ghz"),
                        "chip_length_um": chip_length_um,
                        "chip_width_um": chip_width_um,
                        "chip_height_um": chip_height_um,
                        "source_path": str(file),
                    }

                    processor_config.append(config)
            except Exception as e:
                _log.error("Error reading %s: %s", file, e)

        return processor_config

    # 이름으로 특정 processor config를 찾아 주요 필드를 언팩해 반환
    def target(self, target_name):
        target_chip = next((config for config in self.load_json() if config["topology"] == target_name), None)
        if target_chip is None:
            _log.warning("Target processor '%s' not found.", target_name)
            return None, None, None, None, None

        target_processor = target_chip.get("topology", "Unknown")
        num_qubits = target_chip.get("num_qubits", 0)
        coupling_map = target_chip.get("coupling_map", [])
        freq_ghz = target_chip.get("freq_ghz", [])

        return target_processor, num_qubits, coupling_map, freq_ghz, target_chip

    # 이름 목록에 해당하는 processor config들을 한 번에 일괄 로드
    def load_configs(self, names):
        all_configs = self.load_json()

        by_name = {c['topology']: c for c in all_configs}

        configs = []
        missing = []

        for name in names:
            if name in by_name:
                configs.append(by_name[name])
            else:
                missing.append(name)

        if missing:
            _log.warning("Processors not found in input directory: %s", missing)

        return configs
