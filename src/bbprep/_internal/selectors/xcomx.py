import numpy as np
import stk

from .selector import Selector


class XCOMXSelector(Selector):
    """
    Get position of binder atoms with the molecule centroid in the middle.

    """

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        raise NotImplementedError(
            "This class can only provide positions, because one "
            "position is the centroid of the molecule."
        )

    def get_atomic_positions(
        self,
        molecule: stk.BuildingBlock,
    ) -> tuple[np.ndarray, ...]:
        assert molecule.get_num_functional_groups() == 2

        positions = []
        for fg in molecule.get_functional_groups():
            for id_ in fg.get_bonder_ids():  # type: ignore[attr-defined]
                positions.append(
                    next(molecule.get_atomic_positions(atom_ids=id_))
                )
        positions.insert(1, molecule.get_centroid())

        return tuple(positions)
