#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for ensemble class.

Author: Andrew Tarzia

"""

from dataclasses import dataclass

import stk


@dataclass
class Conformer:
    molecule: stk.BuildingBlock
    conformer_id: int | None = None
    source: str | None = None


class Ensemble:
    def __init__(
        self,
        base_molecule: stk.Molecule,
    ):
        self._base_molecule = base_molecule
        self._molecule_num_atoms = base_molecule.get_num_atoms()
        self._conformers: dict[int, Conformer] = {}

    def get_num_conformers(self):
        return len(self._conformers)

    def add_conformer(self, conformer):
        if conformer.conformer_id is None:
            conf_id = self.get_num_conformers()
        else:
            conf_id = conformer.conformer_id

        assert conf_id not in self._conformers

        self._conformers[conf_id] = conformer

    def yield_conformers(self):
        for conf_id in self._conformers:
            yield self._conformers[conf_id]

    def get_base_molecule(self):
        return self._base_molecule

    def get_molecule_num_atoms(self):
        return self._molecule_num_atoms

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"num_confs={self.get_num_conformers()})"
        )

    def __repr__(self) -> str:
        return str(self)
