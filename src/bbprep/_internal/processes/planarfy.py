import typing

import stk
import stko

from bbprep.generators import ETKDG


class Planarfy:
    """
    Get the most planar molecule in an ensemble.

    """

    def __init__(self, generator: ETKDG, molecule: stk.Molecule):
        """
        Initialise the process.

        """
        self._data: dict[str, float]
        self._molecule = molecule
        self._conformers = {
            conformer_id: conformer
            for conformer_id, conformer in enumerate(
                generator.generate_conformers(molecule)
            )
        }
        self._data = {}

    def _save_to_data(
        self,
        molecule: stk.Molecule,
        conformer_id: int,
        value: float,
    ):
        key = stk.Smiles().get_key(molecule) + f"__{conformer_id}"
        self._data[key] = value

    def _run_process(
        self,
        molecule: stk.Molecule,
        conformer_id: int,
    ) -> float:
        key = stk.Smiles().get_key(molecule) + f"{conformer_id}"
        if key in self._data:
            return self._data[key]

        pc = stko.PlanarityCalculator()
        pc_results = pc.get_results(molecule)
        value = pc_results.get_planarity_parameter()
        self._save_to_data(molecule, conformer_id, value)
        return value

    def get_data(self):
        return self._data

    def calculate_score(
        self,
        molecule: stk.Molecule,
        conformer_id: int,
    ) -> float:
        return self._run_process(molecule, conformer_id)

    def get_all_scores(self) -> typing.Iterable[float]:
        scores = []
        for conformer_id in self._conformers:
            conformer = self._conformers[conformer_id]
            scores.append(self._run_process(conformer, conformer_id))
        return scores

    def get_minimum(self) -> stk.Molecule:
        minimum_score = 1e24
        minimum_molecule = self._molecule.clone()
        for conformer_id in self._conformers:
            conformer = self._conformers[conformer_id]
            score = self._run_process(conformer, conformer_id)
            if score < minimum_score:
                minimum_score = score
                minimum_molecule = conformer.clone()
        return minimum_molecule

    def get_minimum_id(self) -> int:
        minimum_score = 1e24
        minimum_id = -1
        for conformer_id in self._conformers:
            conformer = self._conformers[conformer_id]
            score = self._run_process(conformer, conformer_id)
            if score < minimum_score:
                minimum_score = score
                minimum_id = conformer_id
        return minimum_id
