import stk
from rdkit.Chem import AllChem as rdkit

from .selector import Selector


class BySmartsSelector(Selector):
    """
    Select atom ids in stk molecules by smarts string.

    """

    def __init__(self, smarts: str):
        """
        Initialise Selector.

        """
        self._smarts = smarts

    def select_atoms(self, molecule: stk.Molecule) -> tuple[int]:
        rdkit_mol = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        matches = rdkit_mol.GetSubstructMatches(
            query=rdkit.MolFromSmarts(self._smarts),
        )
        atoms = []
        for match in matches:
            for atom_id in match:
                atoms.append(atom_id)
        return tuple(atoms)
