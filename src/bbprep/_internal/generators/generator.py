import stk

from bbprep._internal.ensemble.ensemble import Ensemble


class Generator:
    """
    Generate an ensemble of molecules.

    """

    def generate_conformers(
        self,
        molecule: stk.BuildingBlock,
    ) -> Ensemble:
        raise NotImplementedError()
