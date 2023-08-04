import bbprep
import stk


class CaseData:
    """
    A test case.

    Attributes:

    """

    def __init__(
        self,
        molecule: stk.Molecule,
        generator: bbprep.generators.Generator,
        selector: bbprep.selectors.Selector,
        min_id: int,
        min_value: float,
        name: str,
    ) -> None:
        self.molecule = molecule
        self.generator = generator
        self.selector = selector
        self.min_id = min_id
        self.min_value = min_value
        self.name = name
