import pytest
import stk

import bbprepared

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3"
            ),
            generator=bbprepared.generators.GeometryScanner(
                target_ranges=(
                    bbprepared.generators.AngleRange(
                        smarts="[#7X2]@[#6X3]@[#6X3H0]-!@[#6X3H0]@[#6X3]",
                        expected_num_atoms=5,
                        scanned_ids=(2, 3, 4),
                        scanned_range=[-5.0, 0, 5.0, 10.0],
                    ),
                    bbprepared.generators.BondRange(
                        smarts="[#7X2]@[#6X3]",
                        expected_num_atoms=2,
                        scanned_ids=(0, 1),
                        scanned_range=[-0.2, 0],
                    ),
                ),
            ),
            num_confs=64,
            min_energy=(72.82838602908487, 23),
            max_energy=(118.0264414866404, 60),
            energy_5=(111.22211705068317, 4),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(smiles="C1=CC=C(C=C1)C2=CC=CC=C2"),
            generator=bbprepared.generators.GeometryScanner(
                target_ranges=(
                    bbprepared.generators.BondRange(
                        smarts="[#6X3H0]-!@[#6X3H0]",
                        expected_num_atoms=2,
                        scanned_ids=(0, 1),
                        scanned_range=[-0.5, -0.3, -0.2, 0, 0.2, 0.3, 0.5],
                    ),
                ),
            ),
            num_confs=7,
            min_energy=(39.329228141023066, 3),
            max_energy=(225.68026567784406, 0),
            energy_5=(46.815866634648344, 4),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=CC(=CC(=C1)Br)Br",
                functional_groups=stk.BromoFactory(deleters=(), bonders=(1,)),
            ),
            generator=bbprepared.generators.SelectorDistanceScanner(
                selector=bbprepared.selectors.BindersSelector(),
                scanned_changes=[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2.0],
            ),
            num_confs=9,
            min_energy=(10.901384927843285, 4),
            max_energy=(170.89862539078675, 8),
            energy_5=(10.901384927843285, 4),
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",  # type: ignore[attr-defined]
    )
