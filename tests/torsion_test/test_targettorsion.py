import os
import pathlib

from bbprep import TargetTorsion
from bbprep.generators import ETKDG


def test_targettorsion(molecule):
    """
    Test :class:`TargetTorsion`.

    Parameters:

        molecule:
            The molecule.

    Returns:

        None : :class:`NoneType`

    """

    ensemble = ETKDG(num_confs=30).generate_conformers(molecule.molecule)

    process = TargetTorsion(
        ensemble=ensemble,
        selector=molecule.selector,
        target_value=molecule.target_value,
    )

    best_molecule = process.get_best()
    path = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    best_molecule.molecule.write(path / f"torsion_{molecule.name}_best.mol")

    all_scores = process.get_all_scores()
    print(all_scores)

    assert molecule.best_id == process.get_best_id()
    assert min(all_scores) == process.calculate_score(
        best_molecule, molecule.best_id
    )
