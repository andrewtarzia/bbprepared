import itertools
import typing
from dataclasses import dataclass

import stk
import stko
from rdkit.Chem import AllChem as rdkit
from rdkit.Chem import rdMolTransforms

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble

from .generator import Generator


@dataclass
class TargetTorsion:
    smarts: str
    expected_num_atoms: int
    torsion_ids: tuple[int, int, int, int]


class TorsionScanner(Generator):
    """
    Generate conformers by scanning target torsions.

    """

    def __init__(
        self,
        target_torsions: TargetTorsion | tuple[TargetTorsion],
        angle_range: typing.Iterable[float],
    ):
        """
        Initialise generator.

        """

        if not isinstance(target_torsions, tuple):
            self._target_torsions = (target_torsions,)
        else:
            self._target_torsions = target_torsions

        self._angle_range = angle_range

    def generate_conformers(
        self,
        molecule: stk.BuildingBlock,
    ) -> Ensemble:
        # Optimise the initial ligand structure.
        molecule = stko.MMFF(ignore_inter_interactions=False).optimize(
            mol=molecule,
        )
        ensemble = Ensemble(base_molecule=molecule)
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        rdkit_properties = rdkit.MMFFGetMoleculeProperties(
            rdkit_molecule, mmffVariant="MMFF94s"
        )

        matched_torsions = {}
        atoms_to_be_constrained = set()
        for target in self._target_torsions:
            matches = rdkit_molecule.GetSubstructMatches(
                query=rdkit.MolFromSmarts(target.smarts),
            )

            for match in matches:
                assert len(match) == target.expected_num_atoms

                if not any(i in atoms_to_be_constrained for i in match):
                    for i in match:
                        atoms_to_be_constrained.add(i)

                    initial_torsion = rdMolTransforms.GetDihedralDeg(
                        rdkit_molecule.GetConformer(0),
                        match[target.torsion_ids[0]],
                        match[target.torsion_ids[1]],
                        match[target.torsion_ids[2]],
                        match[target.torsion_ids[3]],
                    )
                    key = tuple(match[i] for i in target.torsion_ids)
                    matched_torsions[key] = initial_torsion

        cid = 0
        test_molecule = molecule.clone()
        for target_angles in itertools.product(
            self._angle_range, repeat=len(matched_torsions)
        ):
            rdkit_molecule = test_molecule.to_rdkit_mol()
            rdkit.SanitizeMol(rdkit_molecule)
            rdkit_properties = rdkit.MMFFGetMoleculeProperties(rdkit_molecule)
            ff = rdkit.MMFFGetMoleculeForceField(
                rdkit_molecule,
                rdkit_properties,
            )
            for angle, torsion in zip(target_angles, matched_torsions):
                actual_angle = round(matched_torsions[torsion] + angle, 2)
                ff.MMFFAddTorsionConstraint(
                    torsion[0],
                    torsion[1],
                    torsion[2],
                    torsion[3],
                    False,
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
                ),
            )
            cid += 1

        return ensemble

        # XTB CODE:
        # # ff = rdkit.MMFFGetMoleculeForceField(
        # # rdk_mol, rdk_props)
        # # ff.MMFFAddTorsionConstraint(
        # #     torsion_atom_ids[0],
        # #     torsion_atom_ids[1],
        # #     torsion_atom_ids[2],
        # #     torsion_atom_ids[3],
        # #     False,
        # #     angle - 0.2,
        # #     angle + 0.2,
        # #     1000.0,
        # # )
        # # ff.Minimize()

        # # new_conf = rdkit.Conformer(rdk_mol.GetNumAtoms())
        # # for i in range(rdk_mol.GetNumAtoms()):
        # #     new_conf.SetAtomPosition(
        # #         i, (
        # # rdk_mol.GetConformer(-1).GetAtomPosition(i))
        # #     )

        # # Set up directory.
        # xtb_dir = conformer_search_dir / f"{cname}_xtbopt"
        # out_file = xtb_dir / "output.txt"
        # if not os.path.exists(out_file):
        #     if os.path.exists(xtb_dir):
        #         shutil.rmtree(xtb_dir)
        #     os.mkdir(xtb_dir)
        #     init_dir = os.getcwd()
        #     os.chdir(xtb_dir)

        #     # Write constraint.inp file.
        #     # Xtb indexes start at 1. not 0.
        #     atoms1 = (
        #         f"{torsion_atom_ids_1[0]+1},"
        #         f"{torsion_atom_ids_1[1]+1},"
        #         f"{torsion_atom_ids_1[2]+1},"
        #         f"{torsion_atom_ids_1[3]+1}"
        #     )
        #     atoms2 = (
        #         f"{torsion_atom_ids_2[0]+1},"
        #         f"{torsion_atom_ids_2[1]+1},"
        #         f"{torsion_atom_ids_2[2]+1},"
        #         f"{torsion_atom_ids_2[3]+1}"
        #     )
        #     if overlapping:
        #         tcstr1 = f"dihedral: {atoms1},{float(angle1)}\n"
        #         tcstr2 = ""
        #     else:
        #         tcstr1 = f"dihedral: {atoms1},{float(angle1)}\n"
        #         tcstr2 = f"dihedral: {atoms2},{float(angle2)}\n"
        #     with open(xtb_dir / "input.inp", "w") as f:
        #         f.write(
        #             "$constrain\n"
        #             "force constant=1.0\n"
        #             f"{tcstr1}"
        #             f"{tcstr2}"
        #             "$end\n"
        #         )

        #     # Run xTB opt.
        #     conf_mol.write(xtb_dir / "input.xyz")
        #     xtb_path = (
        #         "/home/atarzia/miniconda3/envs/simple_het/bin/" "xtb"
        #     )
        #     cmd = (
        #         "ulimit -s unlimited ; "
        #         f"{xtb_path} "
        #         "input.xyz --gfnff "
        #         "--opt normal --parallel 4 "
        #         "--chrg 0 "
        #         "--input input.inp"
        #     )
        #     with open(out_file, "w") as f:
        #         # Note that sp.call will hold the program until
        #         # completion of the calculation.
        #         sp.call(
        #             cmd,
        #             stdin=sp.PIPE,
        #             stdout=f,
        #             stderr=sp.PIPE,
        #             # Shell is required for complex arguments.
        #             shell=True,
        #         )
        #     os.chdir(init_dir)

        # # Extract structure and energy.
        # # conf_mol = update_from_rdkit_conf(
        # #     stk_mol=opt_mol,
        # #     rdk_mol=rdk_mol,
        # #     conf_id=-1,
        # # )

        # try:
        #     conf_mol = conf_mol.with_structure_from_file(
        #         xtb_dir / "xtbopt.xyz"
        #     )
        #     with open(out_file, "r") as f:
        #         for line in f.readlines():
        #             if "TOTAL ENERGY" in line:
        #                 XtbEnergy = float(line.strip().split()[3])

        # except FileNotFoundError:
        #     logging.info(
        #         f"{xtb_dir} did not finish - ignoring this pair"
        #     )

        #     # This jumps back quite a bit in both angles, to
        #     # avoid just being at a clash point again.
        #     if overlapping:
        #         previous_cname = f"c_{angle1-20}_{int(angle2)}"
        #     else:
        #         previous_cname = f"c_{angle1-20}_{angle2-20}"

        #     previous_xtb_dir = (
        #         conformer_search_dir / f"{previous_cname}_xtbopt"
        #     )
        #     if not os.path.exists(previous_xtb_dir / "xtbopt.xyz"):
        #         if overlapping:
        #             previous_cname = f"c_{angle1-40}_{int(angle2)}"
        #         else:
        #             previous_cname = f"c_{angle1-40}_{angle2-40}"
        #         previous_xtb_dir = (
        #             conformer_search_dir / f"{previous_cname}_xtbopt"
        #         )

        #     if not os.path.exists(previous_xtb_dir / "xtbopt.xyz"):
        #         raise ValueError(
        #             "Had to go back many times, and none found."
        #             " This is broken!"
        #         )
        #     conf_mol = conf_mol.with_structure_from_file(
        #         previous_xtb_dir / "xtbopt.xyz"
        #     )
        #     continue
