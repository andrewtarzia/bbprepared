import numpy as np
import stk

from .modifier import Modifier


class RandomFGs(Modifier):
    """Modify a building block."""

    def modify(  # type: ignore[override]
        self,
        building_block: stk.BuildingBlock,
        desired_functional_groups: int,
        seed: int = 1000,
    ) -> stk.BuildingBlock:
        selected_fgs: list[stk.FunctionalGroup]
        if (
            building_block.get_num_functional_groups()
            == desired_functional_groups
        ):
            return building_block.clone()

        if (
            building_block.get_num_functional_groups()
            < desired_functional_groups
        ):
            msg = (
                f"{building_block} has less functional groups than"
                f" asked for ({desired_functional_groups})."
            )
            raise RuntimeError(msg)

        generator = np.random.default_rng(seed)

        existing_fgs = list(building_block.get_functional_groups())

        selected_fgs = list(
            generator.choice(
                np.array(existing_fgs),
                size=desired_functional_groups,
                replace=False,
            )
        )

        return building_block.with_functional_groups(selected_fgs)
