import stko

from .case_data import CaseData


def test_geometryscanner(molecule: CaseData) -> None:
    """Test :class:`GeometryScanner`.

    Parameters:

        molecule:
            The molecule to scan.

    Returns:
        None : :class:`NoneType`

    """
    ensemble = molecule.generator.generate_conformers(molecule.molecule)
    print(ensemble)
    assert ensemble.get_num_conformers() == molecule.num_confs
    test_energies = tuple(
        [
            stko.MMFFEnergy(ignore_inter_interactions=False).get_energy(
                conformer.molecule
            )
            for conformer in ensemble.yield_conformers()
        ]
    )

    print(
        test_energies.index(min(test_energies)),
        min(test_energies),
        molecule.min_energy,
    )
    assert (
        min(test_energies),
        test_energies.index(min(test_energies)),
    ) == molecule.min_energy
    print(
        test_energies.index(max(test_energies)),
        max(test_energies),
        molecule.max_energy,
    )
    assert (
        max(test_energies),
        test_energies.index(max(test_energies)),
    ) == molecule.max_energy
    print(test_energies[4], molecule.energy_5)
    assert (test_energies[4], 4) == molecule.energy_5
