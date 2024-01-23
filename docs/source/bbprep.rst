bbprep
======




BuildingBlocks
--------------

New BuildingBlock classes to replace `stk.BuildingBlock`.

.. toctree::
  :maxdepth: 1

  PanelBuildingBlock <_autosummary/bbprep.PanelBuildingBlock>


Ensemble
--------

Contains Conformers.

.. toctree::
  :maxdepth: 1

  Conformer <_autosummary/bbprep.Conformer>
  Ensemble <_autosummary/bbprep.Ensemble>

Generators
----------

For generating conformer ensembles.

.. toctree::
  :maxdepth: 1

  Generator <_autosummary/bbprep.generators.Generator>
  ETKDG <_autosummary/bbprep.generators.ETKDG>
  TorsionScanner <_autosummary/bbprep.generators.TorsionScanner>
  GeometryScanner <_autosummary/bbprep.generators.GeometryScanner>
  BondRange <_autosummary/bbprep.generators.BondRange>
  AngleRange <_autosummary/bbprep.generators.AngleRange>
  TorsionRange <_autosummary/bbprep.generators.TorsionRange>


Selectors
---------

For selecting atoms in molecules for further analysis and manipulation.

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

Processes
---------

For processing ensembles to get a particular conformer property.

.. toctree::
  :maxdepth: 1

  Process <_autosummary/bbprep.Process>
  TargetProcess <_autosummary/bbprep.TargetProcess>
  MinimiseAngle <_autosummary/bbprep.MinimiseAngle>
  DitopicFitter <_autosummary/bbprep.DitopicFitter>
  Planarfy <_autosummary/bbprep.Planarfy>
  TargetTorsion <_autosummary/bbprep.TargetTorsion>


Modifiers
---------

For modifying existing conformers to get a particular conformer property.

.. toctree::
  :maxdepth: 1

  Modifier <_autosummary/bbprep.Modifier>
  FurthestFGs <_autosummary/bbprep.FurthestFGs>
  ClosestFGs <_autosummary/bbprep.ClosestFGs>
  ReorientPanel <_autosummary/bbprep.ReorientPanel>
  ReorientC2Panel <_autosummary/bbprep.ReorientC2Panel>
  ReorientC1Panel <_autosummary/bbprep.ReorientC1Panel>