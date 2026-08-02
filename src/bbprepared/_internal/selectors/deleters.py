from collections import abc

import stk

from .selector import Selector


class DeletersSelector(Selector):
    """Select atom ids in stk molecules by deleters."""

    def __init__(self) -> None:
        """Initialise Selector."""

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        if molecule.get_num_functional_groups() == 0:
            msg = "Molecule has zero functional groups."
            raise ValueError(msg)

        atoms = [
            id_
            for fg in molecule.get_functional_groups()
            for id_ in fg.get_deleter_ids()  # type: ignore[attr-defined]
        ]

        return tuple(atoms)

    def yield_stepwise(
        self,
        molecule: stk.BuildingBlock,
    ) -> abc.Iterator[tuple[int, ...]]:
        if molecule.get_num_functional_groups() == 0:
            msg = "Molecule has zero functional groups."
            raise ValueError(msg)

        for fg in molecule.get_functional_groups():
            atoms = list(fg.get_deleter_ids())  # type: ignore[attr-defined]
            yield tuple(atoms)
