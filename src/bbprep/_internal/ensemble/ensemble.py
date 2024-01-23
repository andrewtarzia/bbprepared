# Distributed under the terms of the MIT License.

"""Module for ensemble class.

Author: Andrew Tarzia

"""

from collections import abc
from dataclasses import dataclass

import stk


@dataclass
class Conformer:
    molecule: stk.BuildingBlock
    conformer_id: int
    source: str | None = None
    permutation: dict[tuple[int], float] | None = None


class Ensemble:
    def __init__(
        self,
        base_molecule: stk.BuildingBlock,
    ) -> None:
        self._base_molecule = base_molecule
        self._molecule_num_atoms = base_molecule.get_num_atoms()
        self._conformers: dict[int, Conformer] = {}

    def get_num_conformers(self) -> int:
        return len(self._conformers)

    def add_conformer(self, conformer: Conformer) -> None:
        conf_id = conformer.conformer_id

        if conf_id in self._conformers:
            msg = f"{conf_id} already in conformer list"
            raise RuntimeError(msg)

        self._conformers[conf_id] = conformer

    def yield_conformers(self) -> abc.Iterator[Conformer]:
        for conf_id in self._conformers:
            yield self._conformers[conf_id]

    def get_base_molecule(self) -> stk.BuildingBlock:
        return self._base_molecule

    def get_molecule_num_atoms(self) -> int:
        return self._molecule_num_atoms

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"num_confs={self.get_num_conformers()})"
        )

    def __repr__(self) -> str:
        return str(self)
