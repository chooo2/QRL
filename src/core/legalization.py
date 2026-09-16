from core.state import ChipState


class Legalization:
    def __init__(self, input_state: list[ChipState], output_state: list[ChipState]):
        self.input_state = input_state
        self.output_state = output_state
