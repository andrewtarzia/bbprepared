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
                selected_indices=(0, 1, 2, 3),
            ),
            generator=bbprep.generators.ETKDG(num_confs=20),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(smiles="C[C@@H](C(=O)O)N"),
            target_value=0,
            best_value=0,
            best_id=3,
            selector=bbprep.selectors.BySmartsSelector(
                smarts="[#7][#6][#6]=[#8]",
                selected_indices=(0, 1, 2, 3),
            ),
            generator=bbprep.generators.ETKDG(num_confs=20),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(smiles="C1=CC=NC(=C1)C2=CC=CC=N2"),
            target_value=120,
            best_value=120,
            best_id=3,
            selector=bbprep.selectors.BySmartsSelector(
                smarts="[#7][#6][#6][#7]",
                selected_indices=(0, 1, 2, 3),
            ),
            generator=bbprep.generators.TorsionScanner(
                target_torsions=(
                    bbprep.generators.TargetTorsion(
                        smarts="[#7][#6][#6][#7]",
                        expected_num_atoms=4,
                        torsion_ids=(0, 1, 2, 3),
                    ),
                ),
                angle_range=range(0, 362, 40),
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(smiles="C1=CC=NC(=C1)C2=CC=CC=N2"),
            target_value=-120,
            best_value=-120,
            best_id=6,
            selector=bbprep.selectors.BySmartsSelector(
                smarts="[#6][#7][#6][#6][#7][#6]",
                selected_indices=(1, 2, 3, 4),
            ),
            generator=bbprep.generators.TorsionScanner(
                target_torsions=bbprep.generators.TargetTorsion(
                    smarts="[#7][#6][#6][#7]",
                    expected_num_atoms=4,
                    torsion_ids=(0, 1, 2, 3),
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
