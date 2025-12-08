Orienting two torsions
======================

Can also be found `here <https://gist.github.com/andrewtarzia/f6541f0e0244c30739ee574a52cfa210>`_.

This comes from an issue in automating the alignment of functional groups
either side of imine bonds, where ETKDG just happens to give the configuration
I do not want.

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb = stk.BuildingBlock(
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

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                bb.get_atoms(),
                bb.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in bb.get_bonds()
        ),
    )

.. testcode:: recipe1-test

    import stk
    import bbprep

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

    ensemble = bbprep.generators.ETKDG(num_confs=100).generate_conformers(
        building_block
    )

Here, we use two distinct torsional selection processes either side of the
imine to get the configuration we want. This benefits from processes occuring
in the same order over the same ensemble.

.. testcode:: recipe1-test

    # Select the C-C-N-Br torsion.
    process1 = bbprep.TargetTorsion(
        ensemble=ensemble,
        selector=bbprep.selectors.BySmartsSelector(
            smarts="[#6]~[#6]~[#7]~[#35]",
            selected_indices=(0, 1, 2, 3),
        ),
        target_value=180,
    )
    p1_by_id = process1.get_all_scores_by_id()

    # Select the N-C-C-N torsion.
    process2 = bbprep.TargetTorsion(
        ensemble=ensemble,
        selector=bbprep.selectors.BySmartsSelector(
            smarts="[#7]~[#6]~[#6]~[#7]",
            selected_indices=(0, 1, 2, 3),
        ),
        target_value=0,
    )
    p2_by_id = process2.get_all_scores_by_id()


Then you can iterate over the ensemble with the selections and compute the
scores.

.. testcode:: recipe1-test

    best_score = float("inf")
    best_conformer = bbprep.Conformer(
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
            best_conformer = bbprep.Conformer(
                molecule=conformer.molecule.clone(),
                conformer_id=conformer.conformer_id,
                source=conformer.source,
                permutation=conformer.permutation,
            )
            best_score = sum_score

    # Get the best conformer as an stk.BuildingBlock.
    best_bb = best_conformer.molecule

.. testcode:: recipe1-test
    :hide:

    import numpy as np
    assert np.isclose(best_score, 0.0, atol=1E-2)
    assert best_conformer.conformer_id == 83

The desired conformation:

.. moldoc::

    import moldoc.molecule as molecule
    import stk
    import bbprep

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

    ensemble = bbprep.generators.ETKDG(num_confs=100).generate_conformers(
        building_block
    )

    # Select the C-C-N-Br torsion.
    process1 = bbprep.TargetTorsion(
        ensemble=ensemble,
        selector=bbprep.selectors.BySmartsSelector(
            smarts="[#6]~[#6]~[#7]~[#35]",
            selected_indices=(0, 1, 2, 3),
        ),
        target_value=180,
    )
    p1_by_id = process1.get_all_scores_by_id()

    # Select the N-C-C-N torsion.
    process2 = bbprep.TargetTorsion(
        ensemble=ensemble,
        selector=bbprep.selectors.BySmartsSelector(
            smarts="[#7]~[#6]~[#6]~[#7]",
            selected_indices=(0, 1, 2, 3),
        ),
        target_value=0,
    )
    p2_by_id = process2.get_all_scores_by_id()

    # Iterate over both selected torsions and merge their scoring function.
    best_score = float("inf")
    best_conformer = bbprep.Conformer(
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
            best_conformer = bbprep.Conformer(
                molecule=conformer.molecule.clone(),
                conformer_id=conformer.conformer_id,
                source=conformer.source,
                permutation=conformer.permutation,
            )
            best_score = sum_score

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                best_conformer.molecule.get_atoms(),
                best_conformer.molecule.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in best_conformer.molecule.get_bonds()
        ),
    )
