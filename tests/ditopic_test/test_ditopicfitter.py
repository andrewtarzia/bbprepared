import pathlib

import numpy as np

from bbprep import DitopicFitter

from .case_data import CaseData


def test_ditopicfitter(molecule: CaseData) -> None:
    """Test :class:`DitopicFitter`.

    Parameters:

        molecule:
            The molecule to process.

    Returns:
        None : :class:`NoneType`

    """
    ensemble = molecule.generator.generate_conformers(molecule.molecule)

    process = DitopicFitter(ensemble=ensemble)

    min_molecule = process.get_minimum()
    path = pathlib.Path(__file__).parent
    min_molecule.molecule.write(path / f"ditopic_{molecule.name}_min.mol")

    all_scores = process.get_all_scores()

    assert np.isclose(molecule.min_value, min(all_scores), rtol=0, atol=1e-3)
    assert molecule.min_id == process.get_minimum_id()
    assert min(all_scores) == process.calculate_score(
        min_molecule, molecule.min_id
    )
