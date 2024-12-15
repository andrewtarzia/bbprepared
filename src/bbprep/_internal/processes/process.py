from collections import abc

import stk

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble
from bbprep._internal.selectors.selector import Selector


class Process:
    """Get the conformer that minimises some property."""

    def __init__(self, ensemble: Ensemble, selector: Selector) -> None:
        """Initialise the process."""
        self._data: dict[str, float]
        self._ensemble = ensemble
        self._selector = selector
        self._data = {}

    def _save_to_data(
        self,
        conformer: Conformer,
        conformer_id: int,
        value: float,
    ) -> None:
        key = stk.Smiles().get_key(conformer.molecule) + f"__{conformer_id}"
        self._data[key] = value

    def _run_process(
        self,
        conformer: Conformer,
        conformer_id: int,
    ) -> float:
        raise NotImplementedError

    def get_data(self) -> dict[str, float]:
        return self._data

    def calculate_score(
        self,
        conformer: Conformer,
        conformer_id: int,
    ) -> float:
        return self._run_process(conformer, conformer_id)

    def calculate_all_scores(self) -> None:
        for conformer in self._ensemble.yield_conformers():
            self.calculate_score(conformer, conformer.conformer_id)

    def get_all_scores(self) -> abc.Iterable[float]:
        return [
            self.calculate_score(conformer, conformer.conformer_id)
            for conformer in self._ensemble.yield_conformers()
        ]

    def get_all_scores_by_id(self) -> dict[int, float]:
        return {
            conformer.conformer_id: self.calculate_score(
                conformer, conformer.conformer_id
            )
            for conformer in self._ensemble.yield_conformers()
        }

    def get_minimum(self) -> Conformer:
        minimum_score = float("inf")
        minimum_conformer = Conformer(
            molecule=self._ensemble.get_base_molecule().clone(),
            conformer_id=-1,
            source=None,
            permutation=None,
        )
        for conformer in self._ensemble.yield_conformers():
            score = self._run_process(conformer, conformer.conformer_id)
            if score < minimum_score:
                minimum_score = score
                minimum_conformer = Conformer(
                    molecule=conformer.molecule.clone(),
                    conformer_id=conformer.conformer_id,
                    source=conformer.source,
                    permutation=conformer.permutation,
                )
        return minimum_conformer

    def get_minimum_id(self) -> int:
        minimum_score = float("inf")
        minimum_id = -1
        for conformer in self._ensemble.yield_conformers():
            score = self._run_process(conformer, conformer.conformer_id)
            if score < minimum_score:
                minimum_score = score
                minimum_id = conformer.conformer_id
        return minimum_id


class TargetProcess(Process):
    """Get the conformer closest to some target property."""

    def __init__(
        self,
        ensemble: Ensemble,
        selector: Selector,
        target_value: float,
    ) -> None:
        """Initialise the process."""
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

    def get_best(self) -> Conformer:
        best_score = 1e24
        best_conformer = Conformer(
            molecule=self._ensemble.get_base_molecule().clone(),
            conformer_id=-1,
            source=None,
            permutation=None,
        )
        for conformer in self._ensemble.yield_conformers():
            score = self.calculate_score(conformer, conformer.conformer_id)
            if score < best_score:
                best_score = score
                best_conformer = Conformer(
                    molecule=conformer.molecule.clone(),
                    conformer_id=conformer.conformer_id,
                    source=conformer.source,
                    permutation=conformer.permutation,
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
