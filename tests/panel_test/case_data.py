import stk
from bbprep import Modifier


class CaseData:
    """
    A test case.

    Attributes:

    """

    def __init__(
        self,
        molecule: stk.Molecule,
        orientmethod: Modifier,
        fg_reorder: tuple[int, int, int, int],
        mapping: dict[int, int],
        name: str,
    ) -> None:
        self.molecule = molecule
        self.orientmethod = orientmethod
        self.fg_reorder = fg_reorder
        self.mapping = mapping
        self.name = name
