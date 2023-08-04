import stk

from .selector import Selector


class AllSelector(Selector):
    """
    Select all atom ids.

    """

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        return tuple(i.get_id() for i in molecule.get_atoms())


class AllNonHSelector(Selector):
    """
    Select all atom ids.

    """

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        return tuple(
            i.get_id()
            for i in molecule.get_atoms()
            if i.get_atomic_number() != 1
        )
