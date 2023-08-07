import numpy as np
import stk

from bbprep._internal.ensemble.ensemble import Conformer

from .process import Process
from .utilities import angle_between


class MinimiseAngle(Process):
    """
    Get the molecule with the min/max of a target angle.

    """

    def _run_process(
        self,
        conformer: Conformer,
        conformer_id: int,
    ) -> float:
        key = stk.Smiles().get_key(conformer.molecule) + f"__{conformer_id}"
        if key in self._data:
            return self._data[key]

        atom_positions = self._selector.get_atomic_positions(
            conformer.molecule
        )

        try:
            assert len(atom_positions) == 3
        except AssertionError:
            raise AssertionError(
                f"Selector found {len(atom_positions)} atoms, not 3"
            )

        vectors = (
            atom_positions[0] - atom_positions[1],
            atom_positions[2] - atom_positions[1],
        )
        value = np.degrees(angle_between(*vectors))
        self._save_to_data(conformer, conformer_id, value)
        return self._data[key]
