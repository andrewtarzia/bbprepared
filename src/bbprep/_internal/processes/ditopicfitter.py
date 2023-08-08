import numpy as np
import stk

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble
from bbprep._internal.selectors.binders import BindersSelector
from bbprep._internal.selectors.notplacers import NotPlacersSelector
from bbprep._internal.selectors.selector import NullSelector

from .process import Process
from .utilities import angle_between


class DitopicFitter(Process):
    """
    Get the conformer best for binding as ditopic linker.

    """

    def __init__(self, ensemble: Ensemble):
        """
        Initialise the process.

        """
        self._data: dict[str, float]
        self._ensemble = ensemble
        self._selector = NullSelector()
        self._data = {}

    def _run_process(
        self,
        conformer: Conformer,
        conformer_id: int,
    ) -> float:
        key = stk.Smiles().get_key(conformer.molecule) + f"__{conformer_id}"
        if key in self._data:
            return self._data[key]

        try:
            assert conformer.molecule.get_num_functional_groups() == 2
        except AssertionError:
            raise AssertionError(
                f"Found {conformer.molecule.get_num_functional_groups()}"
                " functional groups, not 2."
            )

        notplacers = NotPlacersSelector().yield_stepwise(conformer.molecule)
        binders = BindersSelector().yield_stepwise(conformer.molecule)

        notplacers_centroids = []
        for pl in notplacers:
            notplacers_centroids.append(conformer.molecule.get_centroid(pl))

        binders_centroids = []
        for pl in binders:
            binders_centroids.append(conformer.molecule.get_centroid(pl))

        vectors = []
        for i, j in zip(notplacers_centroids, binders_centroids):
            vectors.append((j - i) / np.linalg.norm((j - i)))

        value = angle_between(*vectors)
        self._save_to_data(conformer, conformer_id, value)
        return self._data[key]
