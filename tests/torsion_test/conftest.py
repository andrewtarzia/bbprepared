import bbprep
import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock(smiles="C1=CC=NC(=C1)C2=CC=CC=N2"),
            target_value=180,
            best_value=180,
            best_id=15,
            selector=bbprep.selectors.BySmartsSelector(
                smarts="[#7][#6][#6][#7]",
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(smiles="C[C@@H](C(=O)O)N"),
            target_value=0,
            best_value=0,
            best_id=3,
            selector=bbprep.selectors.BySmartsSelector(
                smarts="[#7][#6][#6]=[#8]",
            ),
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
