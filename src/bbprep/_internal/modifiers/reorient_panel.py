import numpy as np
import stk
from stk._internal.topology_graphs.utilities import _FunctionalGroupSorter
from stk._internal.utilities.utilities import (
    get_acute_vector,
    get_plane_normal,
)

from .modifier import Modifier


class PanelBuildingBlock(stk.BuildingBlock):
    def get_bonders(self):
        return {
            i: tuple(fg.get_bonder_ids())[0]
            for i, fg in enumerate(self.get_functional_groups())
        }

    def show_long_axis(self, long_axis, path):
        string = stk.XyzWriter().to_string(self)

        # Long axis is along y direction.
        x_pos = [0, long_axis[0]]
        y_pos = [0, long_axis[1]]
        z_pos = [0, long_axis[2]]

        for x, y, z in zip(x_pos, y_pos, z_pos):
            string += f"Ar {x} {y} {z}\n"

        string = string.split("\n")
        string[0] = str(int(string[0]) + len(x_pos))
        string = "\n".join(string)

        with open(path, "w") as f:
            f.write(string)

    def get_long_axis(self):
        atom_ids = self.get_bonders()
        if len(atom_ids) != 4:
            raise ValueError(f"{self} has too many functional groups.")

        sorted_fg_ids = list(_FunctionalGroupSorter(self).get_items())
        long_axis = None
        max_pair_distance = 0
        for i, fg_id in enumerate(sorted_fg_ids):
            id1_centroid = self.get_centroid(
                atom_ids=(atom_ids[fg_id]),
            )
            if i != 3:
                id2_centroid = self.get_centroid(
                    atom_ids=(atom_ids[sorted_fg_ids[i + 1]]),
                )
            else:
                id2_centroid = self.get_centroid(
                    atom_ids=(atom_ids[sorted_fg_ids[0]]),
                )
            pair_axis = id1_centroid - id2_centroid
            pair_distance = np.linalg.norm(pair_axis)
            if pair_distance > max_pair_distance:
                long_axis = pair_axis.copy()
                max_pair_distance = pair_distance

        return long_axis

    def get_concave_direction(self):
        atom_ids = self.get_bonders()
        if len(atom_ids) != 4:
            raise ValueError(f"{self} has too many functional groups.")

        sorted_fg_ids = list(_FunctionalGroupSorter(self).get_items())
        pair_list = list(
            zip(sorted_fg_ids, sorted_fg_ids[1:] + sorted_fg_ids[:1])
        )
        max_pair_distance = 0
        for fg_id1, fg_id2 in pair_list:
            id1_centroid = self.get_centroid(
                atom_ids=(atom_ids[fg_id1]),
            )
            id2_centroid = self.get_centroid(
                atom_ids=(atom_ids[fg_id2]),
            )
            pair_axis = id1_centroid - id2_centroid
            pair_distance = np.linalg.norm(pair_axis)
            if pair_distance > max_pair_distance:
                max_pair_distance = pair_distance
                long_axis_centroid = np.divide(id1_centroid + id2_centroid, 2)
                for fg_id3, fg_id4 in pair_list:
                    if fg_id3 not in (fg_id1, fg_id2) and fg_id4 not in (
                        fg_id1,
                        fg_id2,
                    ):
                        id3_centroid = self.get_centroid(
                            atom_ids=(atom_ids[fg_id3]),
                        )
                        id4_centroid = self.get_centroid(
                            atom_ids=(atom_ids[fg_id4]),
                        )

                        short_axis_centroid = np.divide(
                            id4_centroid + id3_centroid, 2
                        )

        concave_direction = short_axis_centroid - long_axis_centroid
        return concave_direction


class ReorientPanel(Modifier):
    def reassign_functional_groups(
        self,
        building_block: PanelBuildingBlock,
    ) -> PanelBuildingBlock:
        # Set functional group ordering based on long axis.
        fg_centroids = tuple(
            building_block.get_centroid(
                atom_ids=fg.get_placer_ids(),
            )
            for fg in building_block.get_functional_groups()
        )
        plus_minus_fg_id = tuple(
            i
            for i, cent in enumerate(fg_centroids)
            if cent[0] > 0 and cent[1] < 0
        )[0]
        fg1_id = plus_minus_fg_id
        fg2_id, fg3_id, fg4_id = tuple(
            i
            for i in range(building_block.get_num_functional_groups())
            if i != fg1_id
        )
        new_fgs = tuple(building_block.get_functional_groups())
        building_block = building_block.with_functional_groups(
            functional_groups=(
                new_fgs[fg1_id],
                new_fgs[fg2_id],
                new_fgs[fg3_id],
                new_fgs[fg4_id],
            )
        )

        return building_block


class ReorientC2Panel(ReorientPanel):
    """
    Modify a building block.

    """

    def modify(self, building_block: stk.BuildingBlock) -> PanelBuildingBlock:
        building_block = PanelBuildingBlock.init_from_molecule(
            molecule=building_block,
            functional_groups=building_block.get_functional_groups(),
        )
        try:
            assert building_block.get_num_functional_groups() == 4
        except AssertionError:
            raise AssertionError(
                f"{building_block} does not have 4 functional groups."
            )

        target_coords = (
            np.array(
                (
                    np.array([1, 1, 0]),
                    np.array([1, -1, 0]),
                    np.array([-1, -1, 0]),
                    np.array([-1, 1, 0]),
                )
            )
            * building_block.get_maximum_diameter()
            / 2
        )
        centroid_pos = np.array([0, 0, 0])
        building_block = building_block.with_centroid(
            position=centroid_pos,
            atom_ids=building_block.get_placer_ids(),
        )

        edge_centroid = sum(target_coords) / len(target_coords)
        edge_normal = get_acute_vector(
            reference=edge_centroid,
            vector=get_plane_normal(
                points=np.array(target_coords),
            ),
        )

        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        edge_position = target_coords[0]
        building_block = building_block.with_rotation_to_minimize_angle(
            start=fg_bonder_centroid - centroid_pos,
            target=edge_position - edge_centroid,
            axis=edge_normal,
            origin=centroid_pos,
        )

        # Flatten wrt to xy plane.
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = get_acute_vector(
            reference=core_centroid - centroid_pos,
            vector=normal,
        )
        building_block = building_block.with_rotation_between_vectors(
            start=normal,
            target=np.array([0, 0, 1]),
            origin=centroid_pos,
        )

        # Align long axis of molecule (defined by deleter atoms) with
        # y axis.
        long_axis_vector = building_block.get_long_axis()
        edge_centroid = sum(target_coords) / len(target_coords)
        edge_normal = get_acute_vector(
            reference=edge_centroid,
            vector=get_plane_normal(
                points=np.array(target_coords),
            ),
        )
        building_block = building_block.with_rotation_to_minimize_angle(
            start=long_axis_vector,
            target=np.array([1, 0, 0]),
            axis=edge_normal,
            origin=centroid_pos,
        )

        return self.reassign_functional_groups(building_block)


class ReorientC1Panel(ReorientPanel):
    """
    Modify a building block.

    """

    def modify(self, building_block: stk.BuildingBlock) -> PanelBuildingBlock:
        building_block = PanelBuildingBlock.init_from_molecule(
            molecule=building_block,
            functional_groups=building_block.get_functional_groups(),
        )
        try:
            assert building_block.get_num_functional_groups() == 4
        except AssertionError:
            raise AssertionError(
                f"{building_block} does not have 4 functional groups."
            )

        target_coords = (
            np.array(
                (
                    np.array([1, 1, 0]),
                    np.array([1, -1, 0]),
                    np.array([-1, -1, 0]),
                    np.array([-1, 1, 0]),
                )
            )
            * building_block.get_maximum_diameter()
            / 2
        )
        centroid_pos = np.array([0, 0, 0])
        edge_centroid = sum(target_coords) / len(target_coords)

        building_block = building_block.with_centroid(
            position=centroid_pos,
            atom_ids=building_block.get_placer_ids(),
        )

        edge_normal = get_acute_vector(
            reference=edge_centroid,
            vector=get_plane_normal(
                points=np.array(target_coords),
            ),
        )

        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        edge_position = target_coords[0]
        building_block = building_block.with_rotation_to_minimize_angle(
            start=fg_bonder_centroid - centroid_pos,
            target=edge_position - edge_centroid,
            axis=edge_normal,
            origin=centroid_pos,
        )

        # Flatten wrt to xy plane.
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = get_acute_vector(
            reference=core_centroid - centroid_pos,
            vector=normal,
        )
        building_block = building_block.with_rotation_between_vectors(
            start=normal,
            target=np.array([0, 0, 1]),
            origin=centroid_pos,
        )

        # Align long axis of molecule with the x axis.
        building_block = building_block.with_rotation_to_minimize_angle(
            start=building_block.get_long_axis(),
            target=np.array([1, 0, 0]),
            axis=edge_normal,
            origin=centroid_pos,
        )
        # Align concave axis of molecule with the y axis.
        building_block = building_block.with_rotation_to_minimize_angle(
            start=building_block.get_concave_direction(),
            target=np.array([0, 1, 0]),
            axis=edge_normal,
            origin=centroid_pos,
        )

        return self.reassign_functional_groups(building_block)
