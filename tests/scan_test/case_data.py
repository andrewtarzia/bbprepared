import bbprep
import stk


class CaseData:
    """A test case."""

    def __init__(
        self,
        molecule: stk.Molecule,
        generator: bbprep.generators.Generator,
        num_confs: int,
        min_energy: tuple,
        max_energy: tuple,
        energy_5: tuple,
        name: str,
    ) -> None:
        self.molecule = molecule
        self.generator = generator
        self.num_confs = num_confs
        self.min_energy = min_energy
        self.max_energy = max_energy
        self.energy_5 = energy_5
        self.name = name
