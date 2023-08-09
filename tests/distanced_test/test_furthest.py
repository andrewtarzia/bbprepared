import os
import pathlib

import stk
from bbprep import FurthestFGs


def test_furtherfgs(molecule):
    """
    Test :class:`FurthestFGs`.

    Parameters:

        molecule:
            The molecule to modify.

    Returns:

        None : :class:`NoneType`

    """

    modified = FurthestFGs(molecule.molecule).modify(
        molecule.desired_functional_groups
    )
    modified_fgs = set(modified.get_functional_groups())
    original_fgs = tuple(molecule.molecule.get_functional_groups())
    fg_positions = []
    for idx, fg in enumerate(original_fgs):
        if fg in modified_fgs:
            assert idx in molecule.furthest_ids
            fg_positions.append(
                molecule.molecule.get_centroid(atom_ids=fg.get_placer_ids())
            )

    xyz_string = stk.XyzWriter().to_string(molecule.molecule).split("\n")
    xyz_string[
        0
    ] = f"{int(xyz_string[0]) + molecule.desired_functional_groups}"
    for fgpos in fg_positions:
        xyz_string.append(
            f"He {round(fgpos[0], 2)} {round(fgpos[1], 2)} "
            f"{round(fgpos[2], 2)}"
        )

    path = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    with open(path / f"dist_{molecule.name}_fur.xyz", "w") as f:
        f.write("\n".join(xyz_string))