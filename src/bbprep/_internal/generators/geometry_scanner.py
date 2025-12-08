import itertools as it
from collections import abc

import stk
import stko
from rdkit.Chem import AllChem, rdMolTransforms

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble

from .generator import Generator
from .targets import AngleRange, BondRange, TorsionRange


class GeometryScanner(Generator):
    """Generate conformers by scanning multiple settable geometries."""

    def __init__(
        self,
        target_ranges: abc.Sequence[BondRange | AngleRange | TorsionRange],
    ) -> None:
        """Initialise generator."""
        if not isinstance(target_ranges, tuple):
            msg = "`target_ranges` should be tuple."
            raise TypeError(msg)

        self._target_ranges = target_ranges

    def generate_conformers(  # noqa: C901, PLR0912
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

        matched_changes = {}
        atoms_to_be_constrained_bonds = set()
        atoms_to_be_constrained_angles = set()
        atoms_to_be_constrained_torsions = set()
        for target in self._target_ranges:
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
                if isinstance(target, AngleRange) and not any(
                    i in atoms_to_be_constrained_angles for i in match
                ):
                    for i in match:
                        atoms_to_be_constrained_angles.add(i)
                    initial_value = rdMolTransforms.GetAngleDeg(
                        rdkit_molecule.GetConformer(0),
                        match[target.scanned_ids[0]],
                        match[target.scanned_ids[1]],
                        match[target.scanned_ids[2]],
                    )
                    key = tuple(match[i] for i in target.scanned_ids)
                    matched_changes[key] = [
                        round(initial_value + test, 2)
                        for test in target.scanned_range
                    ]

                elif isinstance(target, TorsionRange) and not any(
                    i in atoms_to_be_constrained_torsions for i in match
                ):
                    for i in match:
                        atoms_to_be_constrained_torsions.add(i)
                    initial_value = rdMolTransforms.GetDihedralDeg(
                        rdkit_molecule.GetConformer(0),
                        match[target.scanned_ids[0]],
                        match[target.scanned_ids[1]],
                        match[target.scanned_ids[2]],
                        match[target.scanned_ids[3]],
                    )
                    key = tuple(match[i] for i in target.scanned_ids)
                    matched_changes[key] = [
                        round(initial_value + test, 2)
                        for test in target.scanned_range
                    ]

                elif isinstance(target, BondRange) and not any(
                    i in atoms_to_be_constrained_bonds for i in match
                ):
                    for i in match:
                        atoms_to_be_constrained_bonds.add(i)
                    initial_value = rdMolTransforms.GetBondLength(
                        rdkit_molecule.GetConformer(0),
                        match[target.scanned_ids[0]],
                        match[target.scanned_ids[1]],
                    )

                    key = tuple(match[i] for i in target.scanned_ids)
                    matched_changes[key] = [
                        round(initial_value + test, 2)
                        for test in target.scanned_range
                    ]

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

                if len(change) == 2:  # noqa: PLR2004
                    # Add constraint.
                    ff.MMFFAddDistanceConstraint(
                        change[0],
                        change[1],
                        False,  # noqa: FBT003
                        actual_value - 0.01,
                        actual_value + 0.01,
                        1.0e5,
                    )

                elif len(change) == 3:  # noqa: PLR2004
                    # Add constraint.
                    ff.MMFFAddAngleConstraint(
                        change[0],
                        change[1],
                        change[2],
                        False,  # noqa: FBT003
                        actual_value - 0.1,
                        actual_value + 0.1,
                        100.0,
                    )

                elif len(change) == 4:  # noqa: PLR2004
                    # Add constraint.
                    ff.MMFFAddTorsionConstraint(
                        change[0],
                        change[1],
                        change[2],
                        change[3],
                        False,  # noqa: FBT003
                        actual_value - 0.1,
                        actual_value + 0.1,
                        100.0,
                    )

            ff.Minimize(maxIts=500)
            pos_mat = rdkit_molecule.GetConformer(-1).GetPositions()
            test_molecule = molecule.with_position_matrix(pos_mat)
            ensemble.add_conformer(
                conformer=Conformer(
                    molecule=test_molecule,
                    conformer_id=cid,
                    source="geomscan",
                    permutation=permutation,
                ),
            )

        return ensemble
