import typing

import stk

from .selector import Selector


class DeletersSelector(Selector):
    """
    Select atom ids in stk molecules by deleters.

    """

    def __init__(self):
        """
        Initialise Selector.

        """
        pass

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        assert molecule.get_num_functional_groups() > 0

        atoms = []
        for fg in molecule.get_functional_groups():
            for id_ in fg.get_deleter_ids():  # type: ignore[attr-defined]
                atoms.append(id_)

        return tuple(atoms)

    def yield_stepwise(
        self,
        molecule: stk.BuildingBlock,
    ) -> typing.Iterator[tuple[int, ...]]:
        assert molecule.get_num_functional_groups() > 0

        for fg in molecule.get_functional_groups():
            atoms = []
            for id_ in fg.get_deleter_ids():  # type: ignore[attr-defined]
                atoms.append(id_)
            yield tuple(atoms)
