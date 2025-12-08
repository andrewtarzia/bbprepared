import pytest
import stk

import bbprep

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3"
            ),
            generator=bbprep.generators.GeometryScanner(
                target_ranges=(
                    bbprep.generators.AngleRange(
                        smarts="[#7X2]@[#6X3]@[#6X3H0]-!@[#6X3H0]@[#6X3]",
                        expected_num_atoms=5,
                        scanned_ids=(2, 3, 4),
                        scanned_range=[-5.0, 0, 5.0, 10.0],
                    ),
                    bbprep.generators.BondRange(
                        smarts="[#7X2]@[#6X3]",
                        expected_num_atoms=2,
                        scanned_ids=(0, 1),
                        scanned_range=[-0.2, 0],
                    ),
                ),
            ),
            num_confs=64,
            min_energy=(72.82838605023622, 23),
            max_energy=(118.02645660051499, 60),
            energy_5=(111.35399540923027, 4),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(smiles="C1=CC=C(C=C1)C2=CC=CC=C2"),
            generator=bbprep.generators.GeometryScanner(
                target_ranges=(
                    bbprep.generators.BondRange(
                        smarts="[#6X3H0]-!@[#6X3H0]",
                        expected_num_atoms=2,
                        scanned_ids=(0, 1),
                        scanned_range=[-0.5, -0.3, -0.2, 0, 0.2, 0.3, 0.5],
                    ),
                ),
            ),
            num_confs=7,
            min_energy=(39.81232878272989, 3),
            max_energy=(225.68003897535004, 0),
            energy_5=(48.52349318923524, 4),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                smiles="C1=CC(=CC(=C1)Br)Br",
                functional_groups=stk.BromoFactory(deleters=(), bonders=(1,)),
            ),
            generator=bbprep.generators.SelectorDistanceScanner(
                selector=bbprep.selectors.BindersSelector(),
                scanned_changes=[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2.0],
            ),
            num_confs=9,
            min_energy=(10.90138492752934, 4),
            max_energy=(170.89861660707226, 8),
            energy_5=(10.90138492752934, 4),
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",  # type: ignore[attr-defined]
    )
