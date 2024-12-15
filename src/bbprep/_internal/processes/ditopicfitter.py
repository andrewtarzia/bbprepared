import numpy as np
import stk
import stko

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble
from bbprep._internal.selectors.binders import BindersSelector
from bbprep._internal.selectors.notplacers import NotPlacersSelector
from bbprep._internal.selectors.selector import NullSelector

from .process import Process


class DitopicFitter(Process):
    r"""Get the conformer best for binding as ditopic linker.

    The process finds the conformer with the minimum angle between binder-COM
    vectors::

        |- - - O <-- binder
        |
        |
        |
        COM
        |
        |
        |
        |- - - O <-- binder


    """

    def __init__(self, ensemble: Ensemble) -> None:
        """Initialise the process."""
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

        expected = 2
        if conformer.molecule.get_num_functional_groups() != expected:
            msg = (
                f"Found {conformer.molecule.get_num_functional_groups()}"
                " functional groups, not 2."
            )
            raise RuntimeError(msg)

        notplacers = NotPlacersSelector().yield_stepwise(conformer.molecule)
        binders = BindersSelector().yield_stepwise(conformer.molecule)

        notplacers_centroids = [
            conformer.molecule.get_centroid(pl) for pl in notplacers
        ]

        binders_centroids = [
            conformer.molecule.get_centroid(pl) for pl in binders
        ]

        vectors = []
        for i, j in zip(notplacers_centroids, binders_centroids, strict=False):
            vectors.append((j - i) / np.linalg.norm(j - i))

        value = stko.vector_angle(*vectors)
        self._save_to_data(conformer, conformer_id, value)
        return self._data[key]
