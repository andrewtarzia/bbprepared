import stk
import stko

from bbprep._internal.ensemble.ensemble import Conformer

from .process import Process


class Planarfy(Process):
    """
    Get the most planar molecule in an ensemble.

    """

    def _run_process(
        self,
        conformer: Conformer,
        conformer_id: int,
    ) -> float:
        key = stk.Smiles().get_key(conformer.molecule) + f"__{conformer_id}"
        if key in self._data:
            return self._data[key]

        pc = stko.PlanarityCalculator()
        pc_results = pc.get_results(
            conformer.molecule,
            plane_atom_ids=self._selector.select_atoms(conformer.molecule),
            deviation_atom_ids=self._selector.select_atoms(conformer.molecule),
        )
        value = pc_results.get_planarity_parameter()
        self._save_to_data(conformer, conformer_id, value)
        return self._data[key]
