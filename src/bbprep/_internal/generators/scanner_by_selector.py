import itertools as it
from collections import abc

import numpy as np
import stk
import stko
from rdkit.Chem import AllChem

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble
from bbprep.selectors import Selector

from .generator import Generator


class SelectorDistanceScanner(Generator):
    """Generate conformers by scanning over one selector."""

    def __init__(
        self,
        selector: Selector,
        scanned_changes: abc.Iterable[float],
    ) -> None:
        """Initialise generator."""
        self._selector = selector
        self._scanned_changes = scanned_changes

    def generate_conformers(
        self,
        molecule: stk.BuildingBlock,
    ) -> Ensemble:
        # Optimise the initial ligand structure.
        molecule = stko.MMFF(  # type:ignore[assignment]
            ignore_inter_interactions=False
        ).optimize(
            mol=molecule,
        )
        ensemble = Ensemble(base_molecule=molecule)
        rdkit_molecule = molecule.to_rdkit_mol()
        AllChem.SanitizeMol(rdkit_molecule)  # type: ignore[attr-defined]
        rdkit_properties = AllChem.MMFFGetMoleculeProperties(  # type: ignore[attr-defined]
            rdkit_molecule, mmffVariant="MMFF94s"
        )

        selected_ids = self._selector.select_atoms(molecule)
        atom_positions = self._selector.get_atomic_positions(molecule)

        initial_value = float(
            np.linalg.norm(atom_positions[1] - atom_positions[0])
        )

        key = tuple(i for i in selected_ids)
        matched_changes = {
            key: [
                round(initial_value + test, 2)
                for test in self._scanned_changes
            ]
        }

        test_molecule = molecule.clone()
        keys, values = zip(*matched_changes.items(), strict=False)
        permutations_dicts = [
            dict(zip(keys, v, strict=False)) for v in it.product(*values)
        ]

        for cid, permutation in enumerate(permutations_dicts):
            rdkit_molecule = test_molecule.to_rdkit_mol()
            AllChem.SanitizeMol(rdkit_molecule)  # type: ignore[attr-defined]
            rdkit_properties = AllChem.MMFFGetMoleculeProperties(  # type: ignore[attr-defined]
                rdkit_molecule
            )
            ff = AllChem.MMFFGetMoleculeForceField(  # type: ignore[attr-defined]
                rdkit_molecule,
                rdkit_properties,
            )
            for change in permutation:
                actual_value = permutation[change]

                # Add constraint.
                ff.MMFFAddDistanceConstraint(
                    change[0],
                    change[1],
                    False,  # noqa: FBT003
                    actual_value - 0.01,
                    actual_value + 0.01,
                    1.0e5,
                )

            ff.Minimize(maxIts=500)
            pos_mat = rdkit_molecule.GetConformer(-1).GetPositions()
            test_molecule = molecule.with_position_matrix(pos_mat)
            ensemble.add_conformer(
                conformer=Conformer(
                    molecule=test_molecule,
                    conformer_id=cid,
                    source="distscan",
                    permutation=permutation,
                ),
            )

        return ensemble
