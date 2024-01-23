import stk


class Modifier:
    """Modify a building block."""

    def __init__(self) -> None:
        """Initialise the process."""

    def modify(self, building_block: stk.BuildingBlock) -> stk.BuildingBlock:
        raise NotImplementedError
