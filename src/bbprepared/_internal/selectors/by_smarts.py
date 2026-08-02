from collections import abc

import stk
from rdkit.Chem import AllChem

from .selector import Selector


class BySmartsSelector(Selector):
    """Select atom ids in stk molecules by smarts string."""

    def __init__(self, smarts: str, selected_indices: tuple[int, ...]) -> None:
        """Initialise Selector."""
        self._smarts = smarts
        self._selected_indices = selected_indices

    def select_atoms(self, molecule: stk.BuildingBlock) -> tuple[int, ...]:
        rdkit_mol = molecule.to_rdkit_mol()
        AllChem.SanitizeMol(rdkit_mol)  # type: ignore[attr-defined]
        matches = rdkit_mol.GetSubstructMatches(
            query=AllChem.MolFromSmarts(self._smarts),  # type: ignore[attr-defined]
        )
        atoms = []
        for match in matches:
            for idx, atom_id in enumerate(match):
                if idx in self._selected_indices:
                    atoms.append(atom_id)
        return tuple(atoms)

    def yield_stepwise(
        self,
        molecule: stk.BuildingBlock,
    ) -> abc.Iterator[tuple[int, ...]]:
        rdkit_mol = molecule.to_rdkit_mol()
        AllChem.SanitizeMol(rdkit_mol)  # type: ignore[attr-defined]
        matches = rdkit_mol.GetSubstructMatches(
            query=AllChem.MolFromSmarts(self._smarts),  # type: ignore[attr-defined]
        )
        for match in matches:
            atoms = []
            for idx, atom_id in enumerate(match):
                if idx in self._selected_indices:
                    atoms.append(atom_id)
            yield tuple(atoms)
