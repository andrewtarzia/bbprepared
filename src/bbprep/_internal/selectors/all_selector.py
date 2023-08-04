import stk

from .selector import Selector


class AllSelector(Selector):
    """
    Select all atom ids.

    """

    def select_atoms(self, molecule: stk.Molecule) -> tuple[int]:
        return tuple(i.get_id() for i in molecule.get_atoms())
