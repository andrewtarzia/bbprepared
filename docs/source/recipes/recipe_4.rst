Ditopic functional group alignment.
===================================

This is a simple script for selecting a building block with two functional
groups aligned.


.. testcode:: recipe6-test

    import stk
    import bbprepared

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


.. moldoc::

    import moldoc.molecule as molecule
    import stk
    import bbprepared

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

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                min_molecule.molecule.get_atoms(),
                min_molecule.molecule.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in min_molecule.molecule.get_bonds()
        ),
    )

.. testcode:: recipe6-test
    :hide:

    import numpy as np

    assert np.isclose(
        process.calculate_score(
            conformer=min_molecule,
            conformer_id=min_molecule.conformer_id,
        ),
        0.8249489663822132,
    )
