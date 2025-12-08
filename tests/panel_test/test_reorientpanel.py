import pathlib

import numpy as np
import stk

from .case_data import CaseData


def test_reorientpanel(molecule: CaseData) -> None:
    """Test :class:`ReorientC2Panel` and :class:`ReorientC1Panel`.

    Parameters:

        molecule:
            The molecule to modify.

    Returns:
        None : :class:`NoneType`

    """
    path = pathlib.Path(__file__).parent
    target_coords = (
        np.array(
            (
                np.array([1, 1, 0]),
                np.array([1, -1, 0]),
                np.array([-1, -1, 0]),
                np.array([-1, 1, 0]),
            )
        )
        * molecule.molecule.get_maximum_diameter()
        / 2
    )
    modified = molecule.orientmethod.modify(molecule.molecule)
    modified.write(path / f"panel_{molecule.name}_final.mol")

    fg_positions = [
        modified.get_centroid(atom_ids=fg.get_placer_ids())
        for fg in modified.get_functional_groups()
    ]

    # Write out the structure.
    xyz_string = stk.XyzWriter().to_string(molecule.molecule).split("\n")
    xyz_string += stk.XyzWriter().to_string(modified).split("\n")[2:]
    xyz_string[0] = f"{int(xyz_string[0]) * 2 + 4}"
    fg_distances = []
    for pos in target_coords:
        xyz_string.append(
            f"He {round(pos[0], 2)} {round(pos[1], 2)} {round(pos[2], 2)}"
        )
        fg_distances.append([np.linalg.norm(i - pos) for i in fg_positions])

    with open(path / f"panel_{molecule.name}.xyz", "w") as f:  # noqa: PTH123
        f.write("\n".join(xyz_string))

    modified_fgs = list(modified.get_functional_groups())
    original_fgs = list(molecule.molecule.get_functional_groups())

    mapping = {}
    for idx, fg in enumerate(modified_fgs):
        old_idx = original_fgs.index(fg)
        mapping[idx] = old_idx
    assert mapping == molecule.mapping

    # Test that the reassigned functional groups match expectations.
    assert tuple(i.index(min(i)) for i in fg_distances) == molecule.fg_reorder
