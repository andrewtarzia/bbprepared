import stk


class Modifier:
    """
    Modify a building block.

    """

    def __init__(self):
        """
        Initialise the process.

        """
        pass

    def modify(self, building_block: stk.BuildingBlock) -> stk.BuildingBlock:
        raise NotImplementedError()
