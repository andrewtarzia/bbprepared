import typing

import stk
import stko

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble


class Planarfy:
    """
    Get the most planar molecule in an ensemble.

    """

    def __init__(self, ensemble: Ensemble):
        """
        Initialise the process.

        """
        self._data: dict[str, float]
        self._ensemble = ensemble
        self._data = {}

    def _save_to_data(
        self,
        conformer: Conformer,
        conformer_id: int,
        value: float,
    ):
        key = stk.Smiles().get_key(conformer.molecule) + f"__{conformer_id}"
        self._data[key] = value

    def _run_process(
        self,
        conformer: Conformer,
        conformer_id: int,
    ) -> float:
        key = stk.Smiles().get_key(conformer.molecule) + f"__{conformer_id}"
        if key in self._data:
            return self._data[key]

        pc = stko.PlanarityCalculator()
        pc_results = pc.get_results(conformer.molecule)
        value = pc_results.get_planarity_parameter()
        self._save_to_data(conformer, conformer_id, value)
        return self._data[key]

    def get_data(self):
        return self._data

    def calculate_score(
        self,
        conformer: Conformer,
        conformer_id: int,
    ) -> float:
        return self._run_process(conformer, conformer_id)

    def get_all_scores(self) -> typing.Iterable[float]:
        scores = []
        for conformer in self._ensemble.yield_conformers():
            scores.append(self._run_process(conformer, conformer.conformer_id))
        return scores

    def get_minimum(self) -> Conformer:
        minimum_score = 1e24
        minimum_conformer = Conformer(
            molecule=self._ensemble.get_base_molecule().clone(),
            conformer_id=-1,
            source=None,
        )
        for conformer in self._ensemble.yield_conformers():
            score = self._run_process(conformer, conformer.conformer_id)
            if score < minimum_score:
                minimum_score = score
                minimum_conformer = Conformer(
                    molecule=conformer.molecule.clone(),
                    conformer_id=conformer.conformer_id,
                    source=conformer.source,
                )
        return minimum_conformer

    def get_minimum_id(self) -> int:
        minimum_score = 1e24
        minimum_id = -1
        for conformer in self._ensemble.yield_conformers():
            score = self._run_process(conformer, conformer.conformer_id)
            if score < minimum_score:
                minimum_score = score
                minimum_id = conformer.conformer_id
        return minimum_id
