import os
import pathlib

import numpy as np
from bbprep import DitopicFitter


def test_ditopicfitter(molecule):
    """
    Test :class:`DitopicFitter`.

    Parameters:

        molecule:
            The molecule to process.

    Returns:

        None : :class:`NoneType`

    """

    ensemble = molecule.generator.generate_conformers(molecule.molecule)

    process = DitopicFitter(ensemble=ensemble)

    min_molecule = process.get_minimum()
    path = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    min_molecule.molecule.write(path / f"ditopic_{molecule.name}_min.mol")

    all_scores = process.get_all_scores()
    print(all_scores)
    print(
        process.get_minimum_id(),
        process.calculate_score(min_molecule, molecule.min_id),
    )

    assert np.isclose(molecule.min_value, min(all_scores), rtol=0, atol=1e-3)
    assert molecule.min_id == process.get_minimum_id()
    assert min(all_scores) == process.calculate_score(
        min_molecule, molecule.min_id
    )
