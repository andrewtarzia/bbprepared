import itertools

import numpy as np
import stk

from .modifier import Modifier


class DistancedFGs(Modifier):
    """
    Modify a building block.

    """

    def __init__(self):
        """
        Initialise the process.

        """
        self._reverse: bool

    def modify(  # type: ignore[override]
        self,
        building_block: stk.BuildingBlock,
        desired_functional_groups: int,
    ) -> stk.BuildingBlock:
        selected_fgs: list[stk.FunctionalGroup]
        if (
            building_block.get_num_functional_groups()
            == desired_functional_groups
        ):
            return building_block.clone()

        try:
            assert (
                building_block.get_num_functional_groups()
                > desired_functional_groups
            )
        except AssertionError:
            raise AssertionError(
                f"{building_block} has more functional groups than"
                f" asked for ({desired_functional_groups})."
            )

        fg_centroids = [
            (
                fg,
                building_block.get_centroid(atom_ids=fg.get_placer_ids()),
            )
            for fg in building_block.get_functional_groups()
        ]
        fg_distances = sorted(
            [
                (i[0], j[0], np.linalg.norm(i[1] - j[1]))
                for i, j in itertools.combinations(fg_centroids, 2)
            ],
            key=lambda x: x[2],
            reverse=self._reverse,
        )

        selected_fgs = []
        for idx, fg_pair in enumerate(fg_distances):
            if len(selected_fgs) == 0:
                selected_fgs.append(fg_distances[0][0])
                selected_fgs.append(fg_distances[0][1])
            else:
                if fg_distances[idx][0] in set(selected_fgs):
                    selected_fgs.append(fg_distances[idx][1])
                elif fg_distances[idx][1] in set(selected_fgs):
                    selected_fgs.append(fg_distances[idx][0])

            if len(selected_fgs) == desired_functional_groups:
                break

        return building_block.with_functional_groups(selected_fgs)


class FurthestFGs(DistancedFGs):
    """
    Modify a building block.

    """

    def __init__(self):
        """
        Initialise the process.

        """
        self._reverse = True


class ClosestFGs(DistancedFGs):
    """
    Modify a building block.

    """

    def __init__(self):
        """
        Initialise the process.

        """
        self._reverse = False
