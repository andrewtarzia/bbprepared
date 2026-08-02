"""Module for ensemble class."""

from collections import abc
from dataclasses import dataclass

import stk

from .calculators import EnergyCalculator, Optimiser


@dataclass
class Conformer:
    molecule: stk.BuildingBlock
    conformer_id: int
    source: str | None = None
    permutation: dict[tuple[int], float] | None = None
    score: int | float | None = None


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

    def get_conformers(self) -> dict[int, Conformer]:
        return self._conformers

    def get_conformer(self, idx: int) -> Conformer:
        try:
            return self._conformers[idx]
        except KeyError as e:
            msg = f"Conformer with ID: {idx} not in ensemble."
            raise KeyError(msg) from e

    def get_base_molecule(self) -> stk.BuildingBlock:
        return self._base_molecule

    def get_molecule_num_atoms(self) -> int:
        return self._molecule_num_atoms

    def optimise_conformers(self, optimiser: Optimiser) -> "Ensemble":
        """Get a new ensemble with optimised conformers."""
        new_ensemble = Ensemble(base_molecule=self._base_molecule)
        for conformer in self.yield_conformers():
            new_ensemble.add_conformer(
                conformer=Conformer(
                    molecule=optimiser.function(conformer.molecule),
                    conformer_id=conformer.conformer_id,
                    source=f"{conformer.source}:{optimiser.name}",
                    permutation=None,
                )
            )
        return new_ensemble

    def get_lowest_energy_conformer(
        self,
        calculator: EnergyCalculator,
    ) -> Conformer:
        """Get the lowest energy conformer based on a calculator.

        The calculator must be a callable that returns a value.
        """
        lowest_energy_conformer: Conformer
        lowest_energy = float("inf")
        for conformer in self.yield_conformers():
            energy = calculator.function(conformer.molecule)
            if energy < lowest_energy:
                lowest_energy = energy
                lowest_energy_conformer = Conformer(
                    molecule=conformer.molecule,
                    conformer_id=conformer.conformer_id,
                    source=conformer.source,
                    permutation=None,
                    score=energy,
                )
        return lowest_energy_conformer

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}(num_confs={self.get_num_conformers()})"
        )

    def __repr__(self) -> str:
        return str(self)
