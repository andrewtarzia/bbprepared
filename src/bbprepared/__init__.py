"""bbprepared module."""

from importlib.metadata import PackageNotFoundError, version
from importlib.util import find_spec

from bbprepared import generators, selectors
from bbprepared._internal.ensemble.calculators import (
    EnergyCalculator,
    Optimiser,
)
from bbprepared._internal.ensemble.ensemble import Conformer, Ensemble
from bbprepared._internal.modifiers.distanced_functional_groups import (
    ClosestFGs,
    FurthestFGs,
)
from bbprepared._internal.modifiers.modifier import Modifier
from bbprepared._internal.modifiers.random_functional_groups import RandomFGs
from bbprepared._internal.modifiers.reorient_panel import (
    PanelBuildingBlock,
    ReorientC1Panel,
    ReorientC2Panel,
    ReorientPanel,
)
from bbprepared._internal.processes.angle import MinimiseAngle
from bbprepared._internal.processes.ditopicfitter import DitopicFitter
from bbprepared._internal.processes.planarfy import Planarfy
from bbprepared._internal.processes.process import Process, TargetProcess
from bbprepared._internal.processes.torsion import TargetTorsion


def _has(pack: str) -> bool:
    try:
        version(pack)
    except PackageNotFoundError:
        return False
    return True


for full, light in (("stk", "stk-lite"), ("stko", "stko-lite")):
    if _has(full) and _has(light):
        msg = (
            f"{full} and {light} are both installed and conflict.\n"
            f"Fix: pip uninstall -y {full} {light} && "
            "pip install 'bbprepared[full]' "
            "(or 'bbprepared[lite]')"
        )
        raise ImportError(msg)

if find_spec("stk") is None or find_spec("stko") is None:
    msg = (
        "Missing dependencies. Install 'bbprepared[full]' or "
        "'bbprepared[lite]'."
    )
    raise ImportError(msg)

__all__ = [
    "ClosestFGs",
    "Conformer",
    "DitopicFitter",
    "EnergyCalculator",
    "Ensemble",
    "FurthestFGs",
    "MinimiseAngle",
    "Modifier",
    "Optimiser",
    "PanelBuildingBlock",
    "Planarfy",
    "Process",
    "RandomFGs",
    "ReorientC1Panel",
    "ReorientC2Panel",
    "ReorientPanel",
    "TargetProcess",
    "TargetTorsion",
    "generators",
    "selectors",
]
