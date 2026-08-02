import numpy as np
import numpy.typing as npt
import stk

from .selector import Selector


class XCOMXSelector(Selector):
    """Get binder atom positions with the molecule centroid in the middle."""

    def select_atoms(
        self,
        molecule: stk.BuildingBlock,
    ) -> tuple[int, ...]:
        msg = (
            "This class can only provide positions, because one "
            "position is the centroid of the molecule."
        )
        raise NotImplementedError(msg)

    def get_atomic_positions(
        self,
        molecule: stk.BuildingBlock,
    ) -> tuple[npt.NDArray[np.float64], ...]:
        if molecule.get_num_functional_groups() != 2:  # noqa: PLR2004
            msg = "Molecule does not have 2 functional groups."
            raise ValueError(msg)

        positions = [
            next(molecule.get_atomic_positions(atom_ids=id_))
            for fg in molecule.get_functional_groups()
            for id_ in fg.get_bonder_ids()  # type: ignore[attr-defined]
        ]

        positions.insert(1, molecule.get_centroid())

        return tuple(positions)
