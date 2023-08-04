import typing

import stk

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble
from bbprep._internal.selectors.selector import Selector


class Process:
    """
    Get the molecule with some property.

    """

    def __init__(self, ensemble: Ensemble, selector: Selector):
        """
        Initialise the process.

        """
        self._data: dict[str, float]
        self._ensemble = ensemble
        self._selector = selector
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
        raise NotImplementedError()

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


class TargetProcess(Process):
    """
    Get the conformer with some target value.

    """

    def __init__(
        self,
        ensemble: Ensemble,
        selector: Selector,
        target_value: float,
    ):
        """
        Initialise the process.

        """
        self._data: dict[str, float]
        self._ensemble = ensemble
        self._selector = selector
        self._target_value = target_value
        self._data = {}

    def calculate_value(
        self,
        conformer: Conformer,
        conformer_id: int,
    ) -> float:
        return self._run_process(conformer, conformer_id)

    def calculate_score(
        self,
        conformer: Conformer,
        conformer_id: int,
    ) -> float:
        return abs(
            self._run_process(conformer, conformer_id) - self._target_value
        )

    def get_all_scores(self) -> typing.Iterable[float]:
        scores = []
        for conformer in self._ensemble.yield_conformers():
            scores.append(
                abs(
                    self._run_process(conformer, conformer.conformer_id)
                    - self._target_value
                )
            )
        return scores

    def get_best(self) -> Conformer:
        best_score = 1e24
        best_conformer = Conformer(
            molecule=self._ensemble.get_base_molecule().clone(),
            conformer_id=-1,
            source=None,
        )
        for conformer in self._ensemble.yield_conformers():
            score = self.calculate_score(conformer, conformer.conformer_id)
            if score < best_score:
                best_score = score
                best_conformer = Conformer(
                    molecule=conformer.molecule.clone(),
                    conformer_id=conformer.conformer_id,
                    source=conformer.source,
                )
        return best_conformer

    def get_best_id(self) -> int:
        best_score = 1e24
        best_id = -1
        for conformer in self._ensemble.yield_conformers():
            score = self.calculate_score(conformer, conformer.conformer_id)
            if score < best_score:
                best_score = score
                best_id = conformer.conformer_id
        return best_id
