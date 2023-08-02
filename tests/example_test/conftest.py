import pytest
from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            name=name,
        ),
    )
)
def example(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
