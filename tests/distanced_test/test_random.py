from bbprep import RandomFGs

from .case_data import CaseData


def test_randomfgs(molecule: CaseData) -> None:
    """Test :class:`RandomFGs`.

    Parameters:

        molecule:
            The molecule to modify.

    Returns:
        None : :class:`NoneType`

    """
    original_fgs = tuple(molecule.molecule.get_functional_groups())
    for seed in (123, 985, 23):
        modified = RandomFGs().modify(
            building_block=molecule.molecule,
            seed=seed,
            desired_functional_groups=molecule.desired_functional_groups,
        )
        modified_fgs = set(modified.get_functional_groups())
        assert len(modified_fgs) == molecule.desired_functional_groups
        for idx, fg in enumerate(original_fgs):
            if fg in modified_fgs:
                assert idx in molecule.random_ids[seed]
