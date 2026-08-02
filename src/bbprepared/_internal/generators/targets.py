from collections import abc
from dataclasses import dataclass


@dataclass
class BondRange:
    smarts: str
    expected_num_atoms: int
    scanned_ids: tuple[int, int]
    scanned_range: abc.Iterable[float]


@dataclass
class AngleRange:
    smarts: str
    expected_num_atoms: int
    scanned_ids: tuple[int, int, int]
    scanned_range: abc.Iterable[float]


@dataclass
class TorsionRange:
    smarts: str
    expected_num_atoms: int
    scanned_ids: tuple[int, int, int, int]
    scanned_range: abc.Iterable[float]
