import bbprep
import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3",
            ),
            min_value=0.8743,
            selector=bbprep.selectors.AllSelector(),
            min_id=0,
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
