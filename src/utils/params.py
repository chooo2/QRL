import os
import json
from collections import OrderedDict

# config/params.json 파라미터를 로드하고 관리하는 클래스
class Params:
    # config/params.json을 읽어 모든 파라미터를 인스턴스 속성으로 로드
    def __init__(self):
        root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        file = os.path.join(root, 'config', 'params.json')
        self.__dict__ = {}
        self._load_from_file(file)

    # {key: {"description":..., "values":...}} 형식 JSON을 읽어 params_dict/속성에 반영
    def _load_from_file(self, file) -> None:
        with open(file, "r") as f:
            params_dict = json.load(f, object_pairs_hook=OrderedDict)
        for key, value in params_dict.items():
            self.__dict__[key] = value['values']
        self.__dict__['params_dict'] = params_dict

    # __init__과 동일한 추출 로직으로 대체 설정 파일(params2.json 등)을 로드
    def load(self, file) -> None:
        self._load_from_file(file)

    # 모든 파라미터의 이름·설명·현재값을 콘솔에 출력
    def printHelp(self):
        print("Setting parameters:")
        param_info = self.list_parameters()
        for key, info in param_info.items():
            desc = info.get('description')
            values = info.get('values')
            print(f"  {key}: {desc} (values: {values})")
        print()

    # 파라미터별 설명/값/현재값을 담은 dict 반환
    def list_parameters(self):
        result = {}
        params_dict = self.__dict__.get('params_dict', {})
        for key, config in params_dict.items():
            result[key] = {
                'description': config.get('description', ''),
                'values': config.get('values'),
                'current': getattr(self, key)
            }
        return result
