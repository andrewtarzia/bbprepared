import os
import pathlib

import numpy as np
from bbprep import MinimiseAngle
from bbprep.generators import ETKDG


def test_minimiseangle(molecule):
    """
    Test :class:`MinimiseAngle`.

    Parameters:

        molecule:
            The molecule to planarfy.

    Returns:

        None : :class:`NoneType`

    """

    ensemble = ETKDG(num_confs=30).generate_conformers(molecule.molecule)

    process = MinimiseAngle(ensemble=ensemble, selector=molecule.selector)

    min_molecule = process.get_minimum()
    path = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    min_molecule.molecule.write(path / f"angle_{molecule.name}_min.mol")

    all_scores = process.get_all_scores()
    print(all_scores)

    assert np.isclose(molecule.min_value, min(all_scores), rtol=0, atol=1e-3)
    assert molecule.min_id == process.get_minimum_id()
    assert min(all_scores) == process.calculate_score(
        min_molecule, molecule.min_id
    )
