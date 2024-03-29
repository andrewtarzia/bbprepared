import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=CC(=CN=C1)C2=NC=C(C=C2)C3=CC=NC=C3",
                functional_groups=stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]",
                    bonders=(1,),
                    deleters=(),
                ),
            ),
            desired_functional_groups=2,
            closest_ids=(0, 1),
            furthest_ids=(0, 2),
            random_ids={123: (0, 2), 985: (1, 0), 23: (0, 2)},
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C(C(Br)Br)Br",
                functional_groups=stk.BromoFactory(),
            ),
            desired_functional_groups=2,
            closest_ids=(1, 2),
            furthest_ids=(0, 1),
            random_ids={123: (0, 2), 985: (1, 0), 23: (0, 2)},
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=CC(=C(C=C1Br)C2=CC(=CC(=C2)Br)Br)Br",
                functional_groups=stk.BromoFactory(),
            ),
            desired_functional_groups=3,
            closest_ids=(0, 2, 3),
            furthest_ids=(1, 2, 3),
            random_ids={123: (0, 2, 3), 985: (1, 0, 2), 23: (0, 1, 2)},
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1(C(C(C(C(C(C1Br)Br)Br)Br)Br)Br)Br",
                functional_groups=stk.BromoFactory(),
            ),
            desired_functional_groups=2,
            closest_ids=(0, 1),
            furthest_ids=(1, 5),
            random_ids={123: (0, 4), 985: (1, 3), 23: (0, 4)},
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
