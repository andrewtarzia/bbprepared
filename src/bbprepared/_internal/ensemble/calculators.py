"""Module for calculator containers."""

from collections import abc
from dataclasses import dataclass

import stk


@dataclass
class EnergyCalculator:
    name: str
    function: abc.Callable[[stk.BuildingBlock], float]


@dataclass
class Optimiser:
    name: str
    function: abc.Callable[[stk.BuildingBlock], stk.BuildingBlock]
