import pytest
import stk

import bbprepared

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1C=CN=CC=1C1=CC=CC(C2=CC=CN=C2)=C1",
                functional_groups=stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]",
                    bonders=(1,),
                    deleters=(),
                ),
            ),
            generator=bbprepared.generators.ETKDG(num_confs=30),
            min_value=0.09384645484464693,
            min_id=22,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1C=CN=CC=1C1=CC=C(C2C=NC=CC=2)C=C1",
                functional_groups=stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]",
                    bonders=(1,),
                    deleters=(),
                ),
            ),
            generator=bbprepared.generators.ETKDG(num_confs=30),
            min_value=0.8128047398711786,
            min_id=13,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1(C=NC=CC=1)C1C=CC=C(C#CC2=CC=CC3=C2C=CN=C3)C=1",
                functional_groups=stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]",
                    bonders=(1,),
                    deleters=(),
                ),
            ),
            generator=bbprepared.generators.ETKDG(num_confs=30),
            min_value=0.810207345088923,
            min_id=5,
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",  # type: ignore[attr-defined]
    )
