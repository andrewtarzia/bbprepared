Modifiers
=========

Classes that modify the properties of existing building blocks, returning
clones, to get a particular property.

.. toctree::
  :maxdepth: 1

  Modifier <_autosummary/bbprep.Modifier>
  FurthestFGs <_autosummary/bbprep.FurthestFGs>
  ClosestFGs <_autosummary/bbprep.ClosestFGs>
  RandomFGs <_autosummary/bbprep.RandomFGs>
  ReorientPanel <_autosummary/bbprep.ReorientPanel>
  ReorientC2Panel <_autosummary/bbprep.ReorientC2Panel>
  ReorientC1Panel <_autosummary/bbprep.ReorientC1Panel>


Example:
--------

These are most useful when you have a building block with more than the desired
number of functional groups, and you want to simply pick the closest, furthest
or random options.

Here, we have 3 ``C-N-C`` groups, but want the furthest 2.

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb = stk.BuildingBlock(
        smiles="C1=CC(=CN=C1)C2=NC=C(C=C2)C3=CC=NC=C3",
        functional_groups=stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
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

.. testcode:: modifiers-test

    import stk
    import bbprep

    bb = stk.BuildingBlock(
        smiles="C1=CC(=CN=C1)C2=NC=C(C=C2)C3=CC=NC=C3",
        functional_groups=stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    )
    modified = bbprep.FurthestFGs().modify(
        building_block=bb,
        desired_functional_groups=2,
    )

.. testcode:: modifiers-test
    :hide:

    new_fg_ids = tuple(
        sorted(
            {j for i in modified.get_functional_groups() for j in i.get_atom_ids()}
        )
    )
    old_fg_ids = tuple(
        sorted({j for i in bb.get_functional_groups() for j in i.get_atom_ids()})
    )

    assert new_fg_ids == (3, 4, 5, 14, 15, 16)
    assert old_fg_ids == (3, 4, 5, 6, 7, 8, 14, 15, 16)
