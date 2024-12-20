import pathlib

import numpy as np

from bbprep import TargetTorsion

from .case_data import CaseData


def test_targettorsion(molecule: CaseData) -> None:
    """Test :class:`TargetTorsion`.

    Parameters:

        molecule:
            The molecule.

    Returns:
        None : :class:`NoneType`

    """
    ensemble = molecule.generator.generate_conformers(molecule.molecule)

    process = TargetTorsion(
        ensemble=ensemble,
        selector=molecule.selector,
        target_value=molecule.target_value,
    )

    best_molecule = process.get_best()
    path = pathlib.Path(__file__).parent
    best_molecule.molecule.write(path / f"torsion_{molecule.name}_best.mol")

    all_scores = process.get_all_scores()

    assert molecule.best_id == process.get_best_id()
    assert min(all_scores) == process.calculate_score(
        best_molecule, molecule.best_id
    )
    assert np.isclose(
        molecule.best_value,
        process.calculate_value(
            best_molecule,
            molecule.best_id,
        ),
        rtol=0,
        atol=5e-1,
    )
