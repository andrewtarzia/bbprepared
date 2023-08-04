import stk

from .selector import Selector


class BySmartsSelector(Selector):
    """
    Select atom ids in stk molecules by smarts string.

    """

    def __init__(self, smarts: str):
        """
        Initialise Selector.

        """
        self._smarts = smarts

    def select_atoms(self, molecule: stk.Molecule) -> tuple[int]:
        print(molecule)
        print(self._smarts)
        raise NotImplementedError()
