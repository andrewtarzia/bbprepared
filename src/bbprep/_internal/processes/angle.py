import numpy as np
import stk
import stko

from bbprep._internal.ensemble.ensemble import Conformer

from .process import Process


class MinimiseAngle(Process):
    """Get the molecule with the min/max of a target angle."""

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

        expected = 3
        if len(atom_positions) != expected:
            msg = f"Selector found {len(atom_positions)} atoms, not 3"
            raise RuntimeError(msg)

        vectors = (
            atom_positions[0] - atom_positions[1],
            atom_positions[2] - atom_positions[1],
        )
        value = np.degrees(stko.vector_angle(*vectors))
        self._save_to_data(conformer, conformer_id, value)
        return self._data[key]
