import typing

import stk

from .selector import Selector


class NotPlacersSelector(Selector):
    """
    Select atom ids in stk molecules by placers.

    """

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        assert molecule.get_num_functional_groups() > 0

        atoms = []
        for fg in molecule.get_functional_groups():
            placers = tuple(fg.get_placer_ids())  # type: ignore[attr-defined]
            for id_ in fg.get_atom_ids():  # type: ignore[attr-defined]
                if id_ not in placers:
                    atoms.append(id_)

        return tuple(atoms)

    def yield_stepwise(
        self,
        molecule: stk.BuildingBlock,
    ) -> typing.Iterator[tuple[int, ...]]:
        assert molecule.get_num_functional_groups() > 0

        for fg in molecule.get_functional_groups():
            atoms = []
            placers = tuple(fg.get_placer_ids())  # type: ignore[attr-defined]
            for id_ in fg.get_atom_ids():  # type: ignore[attr-defined]
                if id_ not in placers:
                    atoms.append(id_)
            yield tuple(atoms)
