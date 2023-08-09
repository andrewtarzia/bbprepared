import numpy as np
import pytest
import stk
from bbprep import ReorientC1Panel, ReorientC2Panel

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=C2C=C(C(=CC2=CC(=C1Br)Br)Br)Br",
                functional_groups=stk.BromoFactory(),
            ).with_rotation_about_axis(
                angle=20,
                axis=np.array((1, 2, 3)),
                origin=np.array((0, 0, 0)),
            ),
            orientmethod=ReorientC2Panel(),
            fg_reorder=(1, 0, 3, 2),
            mapping={0: 0, 1: 1, 2: 2, 3: 3},
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=C2C=C(C(=CC2=CC(=C1Br)Br)Br)Br",
                functional_groups=stk.BromoFactory(),
            ).with_rotation_about_axis(
                angle=20,
                axis=np.array((1, 2, 3)),
                origin=np.array((0, 0, 0)),
            ),
            orientmethod=ReorientC1Panel(),
            fg_reorder=(1, 0, 3, 2),
            mapping={0: 0, 1: 1, 2: 2, 3: 3},
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=C2C=C(C(=CC2=CC(=C1Br)I)Br)I",
                functional_groups=(
                    stk.BromoFactory(),
                    stk.IodoFactory(),
                ),
            ).with_rotation_about_axis(
                angle=50,
                axis=np.array((1, 2, 3)),
                origin=np.array((0, 0, 0)),
            ),
            orientmethod=ReorientC2Panel(),
            fg_reorder=(1, 0, 2, 3),
            mapping={0: 2, 1: 0, 2: 1, 3: 3},
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles=(
                    "C1=CC(=CC=C1N(C2=CC=C(C=C2)Br)C3=CC=C(C=C3)Br)N(C4"
                    "=CC=C(C=C4)Br)C5=CC=C(C=C5)Br"
                ),
                functional_groups=stk.BromoFactory(),
            ).with_rotation_about_axis(
                angle=50,
                axis=np.array((1, 2, 3)),
                origin=np.array((0, 0, 0)),
            ),
            orientmethod=ReorientC2Panel(),
            fg_reorder=(1, 0, 2, 3),
            mapping={0: 1, 1: 0, 2: 2, 3: 3},
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles=(
                    "C1=CC(=CC=C1N(C2=CC=C(C=C2)Br)C3=CC=C(S3)N(C4=CC=C"
                    "(C=C4)Br)C5=CC=C(C=C5)Br)Br"
                ),
                functional_groups=stk.BromoFactory(),
            ).with_rotation_about_axis(
                angle=80,
                axis=np.array((1, 2, 3)),
                origin=np.array((0, 0, 0)),
            ),
            orientmethod=ReorientC1Panel(),
            fg_reorder=(1, 0, 3, 2),
            mapping={0: 1, 1: 0, 2: 2, 3: 3},
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles=(
                    "BrC1C=CC(=CC=1)N(C#CI)C1=CC=C(N(C#CI)C2C=CC(=CC=2)"
                    "Br)S1"
                ),
                functional_groups=(
                    stk.BromoFactory(),
                    stk.IodoFactory(),
                ),
            ).with_rotation_about_axis(
                angle=80,
                axis=np.array((1, 2, 3)),
                origin=np.array((0, 0, 0)),
            ),
            orientmethod=ReorientC1Panel(),
            fg_reorder=(3, 0, 1, 2),
            mapping={0: 1, 1: 0, 2: 2, 3: 3},
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
