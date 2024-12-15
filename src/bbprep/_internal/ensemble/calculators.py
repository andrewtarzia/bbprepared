"""Module for calculator containers."""

from collections import abc
from dataclasses import dataclass


@dataclass
class EnergyCalculator:
    name: str
    function: abc.Callable


@dataclass
class Optimiser:
    name: str
    function: abc.Callable
