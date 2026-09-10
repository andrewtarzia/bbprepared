Set torsions in a BB.
=====================

This is a simple script for finding conformers in an ensemble with a particular
torsion.


We can scan over a torsion.

.. testcode:: recipe3-test

    import stk
    import bbprepared

    bb = stk.BuildingBlock(smiles="C1=CC=NC(=C1)C2=CC=CC=N2")

    selector = bbprepared.selectors.BySmartsSelector(
        smarts="[#6][#7][#6][#6][#7][#6]",
        selected_indices=(1, 2, 3, 4),
    )
    generator = bbprepared.generators.TorsionScanner(
        target_torsions=bbprepared.generators.TorsionRange(
            smarts="[#7][#6][#6][#7]",
            expected_num_atoms=4,
            scanned_ids=(0, 1, 2, 3),
            scanned_range=range(0, 362, 20),
        ),
    )
    ensemble = generator.generate_conformers(bb)

Then we can select specific torsions.

.. testcode:: recipe3-test

    target1 = 120
    process1 = bbprepared.TargetTorsion(
        ensemble=ensemble,
        selector=selector,
        target_value=target1,
    )
    molecule1 = process1.get_best()

    target2 = 60
    process2 = bbprepared.TargetTorsion(
        ensemble=ensemble,
        selector=selector,
        target_value=target2,
    )
    molecule2 = process2.get_best()



.. moldoc::

    import moldoc.molecule as molecule

    import stk
    import bbprepared

    bb = stk.BuildingBlock(smiles="C1=CC=NC(=C1)C2=CC=CC=N2")

    selector = bbprepared.selectors.BySmartsSelector(
        smarts="[#6][#7][#6][#6][#7][#6]",
        selected_indices=(1, 2, 3, 4),
    )
    generator = bbprepared.generators.TorsionScanner(
        target_torsions=bbprepared.generators.TorsionRange(
            smarts="[#7][#6][#6][#7]",
            expected_num_atoms=4,
            scanned_ids=(0, 1, 2, 3),
            scanned_range=range(0, 362, 20),
        ),
    )
    ensemble = generator.generate_conformers(bb)

    target2 = 60
    process2 = bbprepared.TargetTorsion(
        ensemble=ensemble,
        selector=selector,
        target_value=target2,
    )
    molecule2 = process2.get_best()

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                molecule2.molecule.get_atoms(),
                molecule2.molecule.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in molecule2.molecule.get_bonds()
        ),
    )


.. moldoc::

    import moldoc.molecule as molecule

    import stk
    import bbprepared

    bb = stk.BuildingBlock(smiles="C1=CC=NC(=C1)C2=CC=CC=N2")

    selector = bbprepared.selectors.BySmartsSelector(
        smarts="[#6][#7][#6][#6][#7][#6]",
        selected_indices=(1, 2, 3, 4),
    )
    generator = bbprepared.generators.TorsionScanner(
        target_torsions=bbprepared.generators.TorsionRange(
            smarts="[#7][#6][#6][#7]",
            expected_num_atoms=4,
            scanned_ids=(0, 1, 2, 3),
            scanned_range=range(0, 362, 20),
        ),
    )
    ensemble = generator.generate_conformers(bb)
    target1 = 120
    process1 = bbprepared.TargetTorsion(
        ensemble=ensemble,
        selector=selector,
        target_value=target1,
    )
    molecule1 = process1.get_best()

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                molecule1.molecule.get_atoms(),
                molecule1.molecule.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in molecule1.molecule.get_bonds()
        ),
    )





.. testcode:: recipe3-test
    :hide:

    import numpy as np

    assert np.isclose(
        process1.calculate_score(
            conformer=molecule1,
            conformer_id=molecule1.conformer_id,
        ),
        0.10068184660579504,
    )
    assert np.isclose(
        process2.calculate_score(
            conformer=molecule2,
            conformer_id=molecule2.conformer_id,
        ),
        0.10004690858116305,
    )
