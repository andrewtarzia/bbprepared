import bbprep
import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="CCC",
            ),
            generator=bbprep.generators.ETKDG(num_confs=30),
            min_value=104.235,
            selector=bbprep.selectors.AllNonHSelector(),
            min_id=6,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3",
                functional_groups=stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]",
                    bonders=(1,),
                    deleters=(),
                ),
            ),
            generator=bbprep.generators.ETKDG(num_confs=30),
            min_value=119.237,
            selector=bbprep.selectors.XCOMXSelector(),
            min_id=26,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=CC(=CN=C1)C2=CC=C(C=C2)C3=CN=CC=C3",
                functional_groups=stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]",
                    bonders=(1,),
                    deleters=(),
                ),
            ),
            min_value=149.311,
            min_id=41,
            selector=bbprep.selectors.XCOMXSelector(),
            generator=bbprep.generators.TorsionScanner(
                target_torsions=(
                    bbprep.generators.TargetTorsion(
                        smarts="[#7X2]@[#6X3]@[#6X3H0]-!@[#6X3H0]@[#6X3]",
                        expected_num_atoms=5,
                        torsion_ids=(1, 2, 3, 4),
                    ),
                ),
                angle_range=range(0, 362, 40),
            ),
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
