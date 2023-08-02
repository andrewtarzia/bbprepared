import numpy as np
from bbprep import Planarfy
from bbprep.generators import ETKDG


def test_planarfy(molecule):
    """
    Test :class:`Planarfy`.

    Parameters:

        molecule:
            The molecule to planarfy.

    Returns:

        None : :class:`NoneType`

    """

    process = Planarfy(
        generator=ETKDG(num_confs=10),
        molecule=molecule.molecule,
    )

    min_molecule = process.get_minimum()
    all_scores = process.get_all_scores()

    assert min(all_scores) == process.calculate_score(
        min_molecule, molecule.min_id
    )
    assert np.isclose(molecule.min_value, min(all_scores), rtol=0, atol=1e-3)

    assert molecule.min_id == process.get_minimum_id()
