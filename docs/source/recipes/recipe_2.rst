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
provided in :mod:`stko`. And in the new verstion, it is a one-liner!

.. testcode:: recipe2-test

    calculator = bbprep.EnergyCalculator(
        name="MMFFEnergy",
        function=stko.MMFFEnergy().get_energy,
    )

    # Iterate over ensemble without optimisation.
    minimum_conformer = ensemble.get_lowest_energy_conformer(
        calculator=calculator,
    )
    minimum_score_no_opt = minimum_conformer.score

    # With optimisation...
    optimiser = bbprep.Optimiser(
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

.. testcode:: recipe2-test
    :hide:

    assert minimum_score_no_opt == 44.139049720694786
    assert minimum_score_opt == 39.52500653652586
    assert minimum_conformer.conformer_id == 55
