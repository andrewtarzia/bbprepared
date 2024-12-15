Containers
==========

:class:`Ensemble` contain :class:`Conformer`, which contain
:class:`stk.BuildingBlock`.

.. toctree::
  :maxdepth: 1

  Conformer <_autosummary/bbprep.Conformer>
  Ensemble <_autosummary/bbprep.Ensemble>

You can optimise and calculate energies for the ensemble based on
``ensembles.py``, using the following classes to contain user-defined methods
or those in :mod:`stko`

.. toctree::
  :maxdepth: 1

  EnergyCalculator <_autosummary/bbprep.EnergyCalculator>
  Optimiser <_autosummary/bbprep.Optimiser>

Example:
--------

.. testcode:: containers-test

    import stk
    import bbprep

    bb = stk.BuildingBlock(smiles="C1=CC=C(C=C1)C2=CC=CC=C2")
    generator = bbprep.generators.ETKDG(num_confs=30)
    ensemble = generator.generate_conformers(bb)
    conformers = list(ensemble.yield_conformers())
