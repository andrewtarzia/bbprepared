Selectors
=========

Selectors take advantage of :mod:`rdkit` and :mod:`stk` to provide access to
specific subsets of a molecule for manipulation.

.. toctree::
  :maxdepth: 1

  Selector <_autosummary/bbprep.selectors.Selector>
  AllNonHSelector <_autosummary/bbprep.selectors.AllNonHSelector>
  AllSelector <_autosummary/bbprep.selectors.AllSelector>
  BindersSelector <_autosummary/bbprep.selectors.BindersSelector>
  ByIdSelector <_autosummary/bbprep.selectors.ByIdSelector>
  BySmartsSelector <_autosummary/bbprep.selectors.BySmartsSelector>
  DeletersSelector <_autosummary/bbprep.selectors.DeletersSelector>
  NotPlacersSelector <_autosummary/bbprep.selectors.NotPlacersSelector>
  NullSelector <_autosummary/bbprep.selectors.NullSelector>
  XCOMXSelector <_autosummary/bbprep.selectors.XCOMXSelector>

Examples:
---------

To be nonspecific, you can use null or all selectors. And this taps into the
:mod:`stk` functional group interface.

.. testcode:: selector-test

    import stk
    import bbprep

    bb = stk.BuildingBlock(
        smiles="C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3",
        functional_groups=stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    )
    count_all = len(bbprep.selectors.AllSelector().select_atoms(bb))
    count_allnonh = len(bbprep.selectors.AllNonHSelector().select_atoms(bb))
    count_bind = len(bbprep.selectors.BindersSelector().select_atoms(bb))
    count_dele = len(bbprep.selectors.DeletersSelector().select_atoms(bb))
    count_null = len(bbprep.selectors.NullSelector().select_atoms(bb))

.. testcode:: selector-test
    :hide:

    assert count_all == 30
    assert count_allnonh == 18
    assert count_bind == 2
    assert count_dele == 0
    assert count_null == 0

Some selectors provide the atom positions only, when you are selecting based on
a geometrical feature.
``X`` is often the ``binder`` positions in :class:`stk.BuildingBlock`.

.. testcode:: selectors2-test

    import stk
    import bbprep
    import numpy as np

    bb = stk.BuildingBlock(
        smiles="C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3",
        functional_groups=stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    )

    selector = bbprep.selectors.XCOMXSelector()
    positions = selector.get_atomic_positions(bb)

.. testcode:: selectors2-test
    :hide:

    assert np.allclose(
        positions[0], np.array([4.79283128, 0.22960641, 0.7650025])
    )
    assert np.allclose(
        positions[1], np.array([-3.46389584e-15,  2.50910404e-15, 1.04306619e-02])
    )
    assert np.allclose(
        positions[2], np.array([-3.70115826, 1.81330028, -1.16595142])
    )
