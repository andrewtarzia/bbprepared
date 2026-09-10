Generators
==========


Classes for generating conformer ensembles using standard rdkit algorithms
or torsion scans.

.. toctree::
  :maxdepth: 1

  Generator <_autosummary/bbprepared.generators.Generator>
  ETKDG <_autosummary/bbprepared.generators.ETKDG>
  TorsionScanner <_autosummary/bbprepared.generators.TorsionScanner>
  XtbTorsionScanner <_autosummary/bbprepared.generators.XtbTorsionScanner>
  GeometryScanner <_autosummary/bbprepared.generators.GeometryScanner>
  SelectorDistanceScanner <_autosummary/bbprepared.generators.SelectorDistanceScanner>


Example:
--------

We can use geometrical ranges (in the ``GeometryScanner`` classes) to build
ensembles along specific collective variables:

.. toctree::
  :maxdepth: 1

  BondRange <_autosummary/bbprepared.generators.BondRange>
  AngleRange <_autosummary/bbprepared.generators.AngleRange>
  TorsionRange <_autosummary/bbprepared.generators.TorsionRange>


.. testcode:: generators-test

    import stk
    import bbprepared
    import stko

    bb = stk.BuildingBlock(smiles="C1=CC=C(C=C1)C2=CC=CC=C2")

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb = stk.BuildingBlock(smiles="C1=CC=C(C=C1)C2=CC=CC=C2")

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

The desired range is actually a change from the initial configuration and can
be set by the user (using SMARTS strings). Any combination of variables can be
used, although you will start hitting a combinatorial nightmare with too
many ranges. Once you have an ensemble, you can do normal analysis.

.. testcode:: generators-test

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

.. testcode:: generators-test
    :hide:

    import numpy as np

    assert np.isclose(energies[3], 39.81232878272989, atol=1E-4)
    assert np.isclose(min(energies), 39.81232878272989, atol=1E-4)

.. moldoc::

    import moldoc.molecule as molecule
    import stk
    import bbprepared

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
    conformers = list(ensemble.yield_conformers())

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                conformers[0].molecule.get_atoms(),
                conformers[0].molecule.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in conformers[0].molecule.get_bonds()
        ),
    )


Using the `xtb scan method <https://xtb-docs.readthedocs.io/en/latest/scan.html>`_,
``bbprepared`` can now scan torsions with improved accuracy in :class:`bbprepared.generators.XtbTorsionScanner`.

You can install ``xtb`` with ``mamba install xtb`` in a `mamba environment`, or
by downloading the source from the site above.

.. code:: python

    import stk
    import bbprepared

    smiles = "C(C#CC1C=CC=C(Br)C=1)1C=CC=C(Br)C=1"

    building_block = stk.BuildingBlock(
        smiles=smiles,
        functional_groups=(stk.BromoFactory(),),
    )

    # Rotate around alkynes.
    smarts_to_rotate = "[#6X3][#6X3][#6X2H0]#!@[#6X2H0][#6X3][#6X3]"
    selected_indices = (0, 1, 4, 5)
    expected_num_atoms = 6

    selector = bbprep.selectors.BySmartsSelector(
        smarts=smarts_to_rotate,
        selected_indices=selected_indices,
    )

    generator = XtbTorsionScanner(
        target_torsions=(
            bbprep.generators.TorsionRange(
                smarts=smarts_to_rotate,
                expected_num_atoms=expected_num_atoms,
                scanned_ids=selected_indices,
                scanned_range=range(0, 359, 20),
            ),
        ),
        # We cannot test this code because xtb is a dependancy.
        xtb_path=xtb_path,
        output_dir=calculation_dir,
    )
    ensemble = generator.generate_conformers(molecule)

    # This gives an ensemble with conformers that have a score attribute
    # equal to the GFN2-xTB energy in Ha.
    for conformer in ensemble.yield_conformers():
        energy = conformer.score * 2625.5 # to kJ/mol.
        # Plot...
