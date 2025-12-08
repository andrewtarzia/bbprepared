import itertools as it

import stk
import stko
from rdkit.Chem import AllChem, rdMolTransforms

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble

from .generator import Generator
from .targets import TorsionRange


class TorsionScanner(Generator):
    """Generate conformers by scanning target torsions."""

    def __init__(
        self,
        target_torsions: TorsionRange | tuple[TorsionRange],
    ) -> None:
        """Initialise generator."""
        if not isinstance(target_torsions, tuple):
            self._target_torsions = (target_torsions,)
        else:
            self._target_torsions = target_torsions

    def generate_conformers(self, molecule: stk.BuildingBlock) -> Ensemble:
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

        matched_torsions = {}
        atoms_to_be_constrained = set()
        for target in self._target_torsions:
            matches = rdkit_molecule.GetSubstructMatches(
                query=AllChem.MolFromSmarts(target.smarts),  # type: ignore[attr-defined]
            )

            for match in matches:
                if len(match) != target.expected_num_atoms:
                    msg = (
                        f"{len(match)} not as expected ("
                        f"{target.expected_num_atoms})"
                    )
                    raise RuntimeError(msg)

                if not any(i in atoms_to_be_constrained for i in match):
                    for i in match:
                        atoms_to_be_constrained.add(i)

                    initial_torsion = rdMolTransforms.GetDihedralDeg(
                        rdkit_molecule.GetConformer(0),
                        match[target.scanned_ids[0]],
                        match[target.scanned_ids[1]],
                        match[target.scanned_ids[2]],
                        match[target.scanned_ids[3]],
                    )
                    key = tuple(match[i] for i in target.scanned_ids)
                    matched_torsions[key] = [
                        round(initial_torsion + angle, 2)
                        for angle in target.scanned_range
                    ]

        test_molecule = molecule.clone()
        keys, values = zip(*matched_torsions.items(), strict=False)
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
            for torsion in permutation:
                actual_angle = permutation[torsion]

                ff.MMFFAddTorsionConstraint(
                    torsion[0],
                    torsion[1],
                    torsion[2],
                    torsion[3],
                    False,  # noqa: FBT003
                    actual_angle - 0.1,
                    actual_angle + 0.1,
                    100.0,
                )

            ff.Minimize(maxIts=500)
            pos_mat = rdkit_molecule.GetConformer(-1).GetPositions()
            test_molecule = molecule.with_position_matrix(pos_mat)
            ensemble.add_conformer(
                conformer=Conformer(
                    molecule=test_molecule,
                    conformer_id=cid,
                    source="torsionscan",
                    permutation=permutation,
                ),
            )

        return ensemble
