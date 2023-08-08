import typing

import numpy as np
import stk


class Selector:
    """
    Select atom ids in stk molecules by deleters.

    """

    def __init__(self):
        """
        Initialise Selector.

        """
        pass

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        raise NotImplementedError()

    def yield_stepwise(
        self,
        molecule: stk.BuildingBlock,
    ) -> typing.Iterator[tuple[int, ...]]:
        raise NotImplementedError()

    def get_atomic_positions(
        self,
        molecule: stk.BuildingBlock,
    ) -> tuple[np.ndarray, ...]:
        return tuple(
            molecule.get_atomic_positions(atom_ids=self.select_atoms(molecule))
        )


class NullSelector(Selector):
    """
    Selecter that does nothing.

    """

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple:
        return ()
