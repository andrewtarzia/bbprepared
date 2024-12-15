import stk
import stko

from bbprep._internal.ensemble.ensemble import Conformer

from .process import TargetProcess


class TargetTorsion(TargetProcess):
    """Get the molecule with the closest torsion to the target."""

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

        expected = 4
        if len(atom_positions) != expected:
            msg = f"Selector found {len(atom_positions)} atoms, not 4"
            raise RuntimeError(msg)

        value = stko.calculate_dihedral(
            pt1=atom_positions[0],
            pt2=atom_positions[1],
            pt3=atom_positions[2],
            pt4=atom_positions[3],
        )
        self._save_to_data(conformer, conformer_id, value)
        return self._data[key]
