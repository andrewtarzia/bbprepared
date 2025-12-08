from collections import abc

import stk


class CaseData:
    """A test case."""

    def __init__(
        self,
        molecule: stk.BuildingBlock,
        desired_functional_groups: int,
        closest_ids: abc.Sequence[int],
        furthest_ids: abc.Sequence[int],
        random_ids: dict[int, abc.Sequence[int]],
        name: str,
    ) -> None:
        self.molecule = molecule
        self.desired_functional_groups = desired_functional_groups
        self.closest_ids = closest_ids
        self.furthest_ids = furthest_ids
        self.random_ids = random_ids
        self.name = name
