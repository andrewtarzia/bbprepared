import stk


class Selector:
    """
    Select atom ids in stk molecules by deleters.

    """

    def __init__(self):
        """
        Initialise Selector.

        """
        pass

    def select_atoms(self, molecule: stk.Molecule) -> tuple[int]:
        raise NotImplementedError()
