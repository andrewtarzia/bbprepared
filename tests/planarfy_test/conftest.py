import bbprep
import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles=(
                    "C1=CC=C2C=C(C=CC2=C1)C3=CC(=CC=C3)C4=CC=CC5=CC=CC=C54"
                ),
            ),
            generator=bbprep.generators.ETKDG(num_confs=10),
            selector=bbprep.selectors.AllSelector(),
            min_value=0.8743,
            min_id=0,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles=(
                    "C1C(Br)=C2C(C=C(C3C=CC=C(C4C5C(=CC=C(Br)C=5)C=C(Br"
                    ")C=4)C=3)C=C2)=CC=1Br"
                ),
                functional_groups=stk.BromoFactory(),
            ),
            generator=bbprep.generators.ETKDG(num_confs=10),
            selector=bbprep.selectors.BindersSelector(),
            min_value=0.2443,
            min_id=0,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="c1ccccc1",
            ),
            generator=bbprep.generators.ETKDG(num_confs=10),
            selector=bbprep.selectors.AllSelector(),
            min_value=0.0,
            min_id=0,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles=(
                    "C1=CC=C2C=C(C=CC2=C1)C3=CC(=CC=C3)C4=CC=CC5=CC=CC=C54"
                ),
            ),
            generator=bbprep.generators.TorsionScanner(
                target_torsions=bbprep.generators.TargetTorsion(
                    smarts="[#6][#6]-!@[#6][#6]",
                    expected_num_atoms=4,
                    torsion_ids=(0, 1, 2, 3),
                ),
                angle_range=range(0, 362, 40),
            ),
            selector=bbprep.selectors.AllSelector(),
            min_value=0.419,
            min_id=87,
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
