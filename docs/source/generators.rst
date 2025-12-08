Generators
==========


Classes for generating conformer ensembles using standard rdkit algorithms
or torsion scans.

.. toctree::
  :maxdepth: 1

  Generator <_autosummary/bbprep.generators.Generator>
  ETKDG <_autosummary/bbprep.generators.ETKDG>
  TorsionScanner <_autosummary/bbprep.generators.TorsionScanner>
  GeometryScanner <_autosummary/bbprep.generators.GeometryScanner>
  SelectorDistanceScanner <_autosummary/bbprep.generators.SelectorDistanceScanner>


Example:
--------

We can use geometrical ranges (in the ``GeometryScanner`` classes) to build
ensembles along specific collective variables:

.. toctree::
  :maxdepth: 1

  BondRange <_autosummary/bbprep.generators.BondRange>
  AngleRange <_autosummary/bbprep.generators.AngleRange>
  TorsionRange <_autosummary/bbprep.generators.TorsionRange>


.. testcode:: generators-test

    import stk
    import bbprep
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

    generator = bbprep.generators.GeometryScanner(
        target_ranges=(
            bbprep.generators.BondRange(
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
    import bbprep

    bb = stk.BuildingBlock(smiles="C1=CC=C(C=C1)C2=CC=CC=C2")

    generator = bbprep.generators.GeometryScanner(
        target_ranges=(
            bbprep.generators.BondRange(
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
