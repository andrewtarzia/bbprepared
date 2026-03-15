import stk

import bbprepared


class CaseData:
    """A test case."""

    def __init__(
        self,
        molecule: stk.BuildingBlock,
        generator: bbprepared.generators.Generator,
        selector: bbprepared.selectors.Selector,
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
