import stk


class Modifier:
    """
    Modify a building block.

    """

    def __init__(self, building_block: stk.BuildingBlock):
        """
        Initialise the process.

        """
        self._building_block = building_block

    def modify(self) -> stk.BuildingBlock:
        raise NotImplementedError()
