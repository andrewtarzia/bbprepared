import stk

from .selector import Selector


class ByIdSelector(Selector):
    """
    Select atom ids in stk molecules by atom ids.

    """

    def __init__(self, ids: tuple[int]):
        """
        Initialise Selector.

        """
        self._ids = ids

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        return tuple(
            i.get_id() for i in molecule.get_atoms(atom_ids=self._ids)
        )
