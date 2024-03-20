import stk


class CaseData:
    """A test case."""

    def __init__(
        self,
        molecule: stk.Molecule,
        desired_functional_groups: int,
        closest_ids: tuple[int, int],
        furthest_ids: tuple[int, int],
        random_ids: dict[int, tuple[int, int]],
        name: str,
    ) -> None:
        self.molecule = molecule
        self.desired_functional_groups = desired_functional_groups
        self.closest_ids = closest_ids
        self.furthest_ids = furthest_ids
        self.random_ids = random_ids
        self.name = name
