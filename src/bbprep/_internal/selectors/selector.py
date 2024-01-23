from collections import abc

import numpy as np
import numpy.typing as npt
import stk


class Selector:
    """Select atom ids in stk molecules by deleters."""

    def __init__(self) -> None:
        """Initialise Selector."""

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        raise NotImplementedError

    def yield_stepwise(
        self,
        molecule: stk.BuildingBlock,
    ) -> abc.Iterator[tuple[int, ...]]:
        raise NotImplementedError

    def get_atomic_positions(
        self,
        molecule: stk.BuildingBlock,
    ) -> tuple[npt.NDArray[np.float64], ...]:
        return tuple(
            molecule.get_atomic_positions(atom_ids=self.select_atoms(molecule))
        )


class NullSelector(Selector):
    """Selecter that does nothing."""

    def select_atoms(
        self,
        molecule: stk.BuildingBlock,  # noqa: ARG002
    ) -> tuple[int, ...]:
        return ()
