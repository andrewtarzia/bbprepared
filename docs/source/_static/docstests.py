"""So I can get the values for the doc tests."""

import logging

import numpy as np
import rmsd
import stk
import stko

import bbprepared

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def log_print(value) -> None:  # noqa: ANN001
    """Handy function."""
    logger.info(value)


def recipe1() -> None:
    logger.info("recipe1.rst")

    building_block = stk.BuildingBlock(
        smiles="C1=CC=NC(=C1)C=NBr",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#35]",
                bonders=(1,),
                deleters=(),
            ),
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )

    ensemble = bbprepared.generators.ETKDG(num_confs=100).generate_conformers(
        building_block
    )

    # Select the C-C-N-Br torsion.
    process1 = bbprepared.TargetTorsion(
        ensemble=ensemble,
        selector=bbprepared.selectors.BySmartsSelector(
            smarts="[#6]~[#6]~[#7]~[#35]",
            selected_indices=(0, 1, 2, 3),
        ),
        target_value=180,
    )
    p1_by_id = process1.get_all_scores_by_id()

    # Select the N-C-C-N torsion.
    process2 = bbprepared.TargetTorsion(
        ensemble=ensemble,
        selector=bbprepared.selectors.BySmartsSelector(
            smarts="[#7]~[#6]~[#6]~[#7]",
            selected_indices=(0, 1, 2, 3),
        ),
        target_value=0,
    )
    p2_by_id = process2.get_all_scores_by_id()

    best_score = float("inf")
    best_conformer = bbprepared.Conformer(
        molecule=ensemble.get_base_molecule().clone(),
        conformer_id=-1,
        source=None,
        permutation=None,
    )
    for conformer in ensemble.yield_conformers():
        p1score = p1_by_id[conformer.conformer_id]
        p2score = p2_by_id[conformer.conformer_id]
        sum_score = p1score + p2score
        if sum_score < best_score:
            best_conformer = bbprepared.Conformer(
                molecule=conformer.molecule.clone(),
                conformer_id=conformer.conformer_id,
                source=conformer.source,
                permutation=conformer.permutation,
            )
            best_score = sum_score

    # Get the best conformer as an stk.BuildingBlock.
    log_print(best_score)
    log_print(best_conformer.conformer_id)


def recipe2() -> None:
    logger.info("recipe2.rst")

    building_block = stk.BuildingBlock(
        smiles="C1=CC=NC(=C1)C=NBr",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#35]",
                bonders=(1,),
                deleters=(),
            ),
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )

    # This uses the rdkit conformer generation.
    ensemble = bbprepared.generators.ETKDG(num_confs=100).generate_conformers(
        building_block
    )

    calculator = bbprepared.EnergyCalculator(
        name="MMFFEnergy",
        function=stko.MMFFEnergy().get_energy,
    )

    # Iterate over ensemble without optimisation.
    minimum_conformer = ensemble.get_lowest_energy_conformer(
        calculator=calculator,
    )
    minimum_score_no_opt = minimum_conformer.score

    # With optimisation...
    optimiser = bbprepared.Optimiser(
        name="MMFF",
        function=stko.MMFF().optimize,
    )
    new_ensemble = ensemble.optimise_conformers(
        optimiser=optimiser,
    )
    minimum_conformer = new_ensemble.get_lowest_energy_conformer(
        calculator=calculator,
    )
    minimum_score_opt = minimum_conformer.score

    log_print(minimum_score_no_opt)
    log_print(minimum_score_opt)
    log_print(minimum_conformer.conformer_id)


def recipe4() -> None:
    logger.info("recipe4.rst")

    bb = stk.BuildingBlock(
        smiles="C1C=CN=CC=1C1=CC=C(C2C=NC=CC=2)C=C1",
        functional_groups=stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    )

    generator = bbprepared.generators.ETKDG(num_confs=30)
    ensemble = generator.generate_conformers(bb)

    process = bbprepared.DitopicFitter(ensemble=ensemble)
    min_molecule = process.get_minimum()

    # As a stk.BuildingBlock
    bb = min_molecule.molecule
    log_print(
        process.calculate_score(
            conformer=min_molecule,
            conformer_id=min_molecule.conformer_id,
        )
    )


def recipe5() -> None:
    logger.info("recipe5.rst")

    smiles = "C(C#CC1C=CC=C(Br)C=1)1C=CC=C(Br)C=1"

    building_block = stk.BuildingBlock(
        smiles=smiles,
        functional_groups=(stk.BromoFactory(),),
    )

    # Rotate around alkynes.
    smarts_to_rotate = "[#6X3][#6X3][#6X2H0]#!@[#6X2H0][#6X3][#6X3]"
    selected_indices = (0, 1, 4, 5)
    expected_num_atoms = 6

    selector = bbprepared.selectors.BySmartsSelector(
        smarts=smarts_to_rotate,
        selected_indices=selected_indices,
    )

    # Scan with MMFF.
    generator = bbprepared.generators.TorsionScanner(
        target_torsions=(
            bbprepared.generators.TorsionRange(
                smarts=smarts_to_rotate,
                expected_num_atoms=expected_num_atoms,
                scanned_ids=selected_indices,
                scanned_range=range(0, 359, 20),
            ),
        ),
    )
    ensemble = generator.generate_conformers(building_block)

    energies = []
    torsions = []
    for conformer in ensemble.yield_conformers():
        molecule = conformer.molecule.with_centroid(np.array((0, 0, 0)))

        # Align to conformer 0.
        if conformer.conformer_id == 0:
            first_pos_mat = molecule.get_position_matrix()
            aligned_ = molecule.clone()

        # Else, align.
        else:
            current_pos_mat = molecule.get_position_matrix()
            alignment = rmsd.kabsch(current_pos_mat, first_pos_mat)
            aligned_pos_mat = np.dot(current_pos_mat, alignment)
            aligned_ = molecule.with_position_matrix(aligned_pos_mat)

        energy = (
            stko.MMFFEnergy(ignore_inter_interactions=False).get_energy(
                aligned_
            )
            * 4.184
        )
        energies.append(energy)

        atom_positions = aligned_.get_position_matrix()
        atoms_to_be_constrained = set()
        for found in selector.yield_stepwise(aligned_):
            if not any(i in atoms_to_be_constrained for i in found):
                for i in found:
                    atoms_to_be_constrained.add(i)

                torsion = float(
                    stko.calculate_dihedral(
                        pt1=atom_positions[found[0]],
                        pt2=atom_positions[found[1]],
                        pt3=atom_positions[found[2]],
                        pt4=atom_positions[found[3]],
                    )
                )

                torsions.append(torsion)

    log_print(torsions[0])
    log_print(energies[0])


def generators() -> None:
    logger.info("generators.rst")
    bb = stk.BuildingBlock(smiles="C1=CC=C(C=C1)C2=CC=CC=C2")

    generator = bbprepared.generators.GeometryScanner(
        target_ranges=(
            bbprepared.generators.BondRange(
                smarts="[#6X3H0]-!@[#6X3H0]",
                expected_num_atoms=2,
                scanned_ids=(0, 1),
                scanned_range=[-0.5, -0.3, -0.2, 0, 0.2, 0.3, 0.5],
            ),
        ),
    )
    ensemble = generator.generate_conformers(bb)
    energies = [
        stko.MMFFEnergy(ignore_inter_interactions=False).get_energy(
            conformer.molecule
        )
        for conformer in ensemble.yield_conformers()
    ]

    log_print(energies[3])
    log_print(min(energies))


def processes() -> None:
    logger.info("processes.rst")

    bb = stk.BuildingBlock(
        smiles=("C1=CC=C2C=C(C=CC2=C1)C3=CC(=CC=C3)C4=CC=CC5=CC=CC=C54"),
    )

    selector = bbprepared.selectors.AllSelector()
    # This is not the cheapest approach, but works well in this case!
    generator = bbprepared.generators.TorsionScanner(
        target_torsions=bbprepared.generators.TorsionRange(
            smarts="[#6][#6]-!@[#6][#6]",
            expected_num_atoms=4,
            scanned_ids=(0, 1, 2, 3),
            scanned_range=range(0, 362, 60),
        ),
    )
    ensemble = generator.generate_conformers(bb)
    process = bbprepared.Planarfy(ensemble=ensemble, selector=selector)

    # Get the minimum bbprepared.Conformer.
    min_molecule = process.get_minimum()

    # Can access the stk molecule using:
    # min_molecule.molecule, and the score as:
    min_score = process.calculate_score(min_molecule, process.get_minimum_id())
    log_print(min_score)


def selectors() -> None:
    logger.info("selectors.rst")

    bb = stk.BuildingBlock(
        smiles="C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3",
        functional_groups=stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    )
    count_all = len(bbprepared.selectors.AllSelector().select_atoms(bb))
    count_allnonh = len(
        bbprepared.selectors.AllNonHSelector().select_atoms(bb)
    )
    count_bind = len(bbprepared.selectors.BindersSelector().select_atoms(bb))
    count_dele = len(bbprepared.selectors.DeletersSelector().select_atoms(bb))
    count_null = len(bbprepared.selectors.NullSelector().select_atoms(bb))

    log_print(count_all)
    log_print(count_allnonh)
    log_print(count_bind)
    log_print(count_dele)
    log_print(count_null)

    bb = stk.BuildingBlock(
        smiles="C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3",
        functional_groups=stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    )

    selector = bbprepared.selectors.XCOMXSelector()
    positions = selector.get_atomic_positions(bb)

    log_print(positions[0])
    log_print(positions[1])
    log_print(positions[2])


def main() -> None:
    selectors()
    generators()
    processes()
    recipe1()
    recipe2()
    recipe4()
    recipe5()


if __name__ == "__main__":
    main()
