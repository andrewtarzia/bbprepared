import stk

import bbprep


class CaseData:
    """A test case."""

    def __init__(
        self,
        molecule: stk.BuildingBlock,
        generator: bbprep.generators.Generator,
        num_confs: int,
        min_energy: tuple[float, int],
        max_energy: tuple[float, int],
        energy_5: tuple[float, int],
        name: str,
    ) -> None:
        self.molecule = molecule
        self.generator = generator
        self.num_confs = num_confs
        self.min_energy = min_energy
        self.max_energy = max_energy
        self.energy_5 = energy_5
        self.name = name
