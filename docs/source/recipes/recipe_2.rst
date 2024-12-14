Getting the lowest energy conformer.
====================================

This is a simple script for iterating over an ensemble and getting the lowest
energy :class:`bbprep.Conformer`.


.. testcode:: recipe2-test

    import stk
    import stko
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

    # This uses the rdkit conformer generation.
    ensemble = bbprep.generators.ETKDG(num_confs=100).generate_conformers(
        building_block
    )

Note that you could couple this with any energy function, especially those
provided in :mod:`stko`.

.. testcode:: recipe2-test

    # Iterate over ensemble.
    minimum_score = 1e24
    minimum_conformer = bbprep.Conformer(
        molecule=ensemble.get_base_molecule().clone(),
        conformer_id=-1,
        source=None,
        permutation=None,
    )
    for conformer in ensemble.yield_conformers():
        # Something here to calculate energy.
        score = stko.MMFFEnergy().get_energy(conformer.molecule)
        if score < minimum_score:
            minimum_score = score
            minimum_conformer = bbprep.Conformer(
                molecule=conformer.molecule.clone(),
                conformer_id=conformer.conformer_id,
                source=conformer.source,
                permutation=conformer.permutation,
            )

.. testcode:: recipe2-test
    :hide:

    assert minimum_score == 43.2237165419001
    assert minimum_conformer.conformer_id == 34

