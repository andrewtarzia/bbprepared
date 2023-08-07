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
        target_value: float,
        best_id: int,
        best_value: float,
        name: str,
    ) -> None:
        self.molecule = molecule
        self.generator = generator
        self.selector = selector
        self.target_value = target_value
        self.best_id = best_id
        self.best_value = best_value
        self.name = name
