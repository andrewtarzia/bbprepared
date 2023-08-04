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

    ensemble = ETKDG(num_confs=10).generate_conformers(molecule.molecule)

    process = MinimiseAngle(ensemble=ensemble, selector=molecule.selector)

    min_molecule = process.get_minimum()
    all_scores = process.get_all_scores()
    print(all_scores)

    assert min(all_scores) == process.calculate_score(
        min_molecule, molecule.min_id
    )
    assert np.isclose(molecule.min_value, min(all_scores), rtol=0, atol=1e-3)

    assert molecule.min_id == process.get_minimum_id()
