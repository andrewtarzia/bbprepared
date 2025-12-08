Processes
=========

A series of classes for processing ensembles to get a particular conformer
property.

.. toctree::
  :maxdepth: 1

  Process (to target the minimum) <_autosummary/bbprep.Process>
  TargetProcess (for a specific target value) <_autosummary/bbprep.TargetProcess>
  MinimiseAngle <_autosummary/bbprep.MinimiseAngle>
  DitopicFitter <_autosummary/bbprep.DitopicFitter>
  Planarfy <_autosummary/bbprep.Planarfy>
  TargetTorsion <_autosummary/bbprep.TargetTorsion>


Example:
--------

One example is to make a building block more planar based on a
plane-of-best-fit through desired atoms.

.. testcode:: processes-test

    import stk
    import bbprep

    bb = stk.BuildingBlock(
        smiles=(
            "C1=CC=C2C=C(C=CC2=C1)C3=CC(=CC=C3)C4=CC=CC5=CC=CC=C54"
        ),
    )

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb = stk.BuildingBlock(
        smiles=(
            "C1=CC=C2C=C(C=CC2=C1)C3=CC(=CC=C3)C4=CC=CC5=CC=CC=C54"
        ),
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

The desired atoms are set by the user, using a selector (here, all atoms!),
and the best is chosen from an ensemble built with the ETKDG algorithm.
The `gists <https://gist.github.com/andrewtarzia>`_ have other ways of
generating ensembles.

.. testcode:: processes-test

    selector = bbprep.selectors.AllSelector()
    # This is not the cheapest approach, but works well in this case!
    generator = bbprep.generators.TorsionScanner(
        target_torsions=bbprep.generators.TorsionRange(
            smarts="[#6][#6]-!@[#6][#6]",
            expected_num_atoms=4,
            scanned_ids=(0, 1, 2, 3),
            scanned_range=range(0, 362, 60),
        ),
    )
    ensemble = generator.generate_conformers(bb)
    process = bbprep.Planarfy(ensemble=ensemble, selector=selector)

    # Get the minimum bbprep.Conformer.
    min_molecule = process.get_minimum()

    # Can access the stk molecule using:
    # min_molecule.molecule, and the score as:
    min_score = process.calculate_score(
        min_molecule, process.get_minimum_id()
    )

.. testcode:: processes-test
    :hide:

    assert min_score == 0.6451099779697566

.. moldoc::

    import moldoc.molecule as molecule
    import stk
    import bbprep

    bb = stk.BuildingBlock(
        smiles=(
            "C1=CC=C2C=C(C=CC2=C1)C3=CC(=CC=C3)C4=CC=CC5=CC=CC=C54"
        ),
    )
    selector = bbprep.selectors.AllSelector()
    generator = bbprep.generators.TorsionScanner(
        target_torsions=bbprep.generators.TorsionRange(
            smarts="[#6][#6]-!@[#6][#6]",
            expected_num_atoms=4,
            scanned_ids=(0, 1, 2, 3),
            scanned_range=range(0, 362, 60),
        ),
    )
    ensemble = generator.generate_conformers(bb)
    process = bbprep.Planarfy(ensemble=ensemble, selector=selector)

    # Get the minimum bbprep.Conformer.
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
