import logging
import os
import pathlib
import shutil
import subprocess as sp

import numpy as np
import stk
from rdkit.Chem import AllChem, rdMolTransforms

from bbprepared._internal.ensemble.ensemble import Conformer, Ensemble

from .generator import Generator
from .targets import TorsionRange

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


class XtbTorsionScanner(Generator):
    """Generate conformers by scanning target torsions with xtb."""

    def __init__(  # noqa: PLR0913, PLR0917
        self,
        target_torsions: TorsionRange | tuple[TorsionRange],
        output_dir: pathlib.Path,
        xtb_path: pathlib.Path,
        gfn_version: int = 2,
        opt_level: str = "normal",
        num_cores: int = 1,
        electronic_temperature: float = 300,
        solvent_model: str = "gbsa",
        solvent: str | None = None,
        solvent_grid: str = "normal",
        charge: int = 0,
        num_unpaired_electrons: int = 0,
        unlimited_memory: bool = False,  # noqa: FBT001, FBT002
        concerted: bool = False,  # noqa: FBT001, FBT002
    ) -> None:
        """Initialise generator."""
        if not isinstance(target_torsions, tuple):
            self._target_torsions = (target_torsions,)
        else:
            self._target_torsions = target_torsions
        self._check_path(xtb_path)
        self._xtb_path = xtb_path
        self._gfn_version = str(gfn_version)
        self._output_dir = output_dir.resolve()
        self._opt_level = opt_level
        self._num_cores = str(num_cores)
        self._electronic_temperature = str(electronic_temperature)
        self._solvent = solvent
        self._solvent_model = solvent_model
        self._solvent_grid = solvent_grid
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._unlimited_memory = unlimited_memory
        self._concerted = concerted

        if self._output_dir.exists():
            shutil.rmtree(self._output_dir)
        self._output_dir.mkdir(parents=True)

        if (
            self._concerted
            and {len(list(i.scanned_range)) for i in self._target_torsions}
            != 1
        ):
            msg = (
                "A concerted scan can only carried out if all "
                "constraints are scanned with the same number of steps."
            )
            raise RuntimeError(msg)

    def _check_path(self, path: pathlib.Path) -> None:
        path = pathlib.Path(path)
        if not path.exists():
            msg = f"XTB not found at {path}"
            raise RuntimeError(msg)

    def _run_xtb(self, xyz: pathlib.Path, out_file: pathlib.Path) -> None:
        """Run GFN-xTB."""
        out_file = pathlib.Path(out_file)

        # Modify the memory limit.
        memory = "ulimit -s unlimited ;" if self._unlimited_memory else ""

        # Set optimization level and type.
        # Do optimization.
        optimization = f"--opt {self._opt_level}"

        if self._solvent is not None:
            solvent = f"--{self._solvent_model} {self._solvent} "
        else:
            solvent = ""

        cmd = (
            f"{memory} {self._xtb_path} {xyz} "
            f"--gfn {self._gfn_version} "
            f"{optimization} --parallel {self._num_cores} "
            f"--etemp {self._electronic_temperature} "
            f"{solvent} --chrg {self._charge} "
            f"--uhf {self._num_unpaired_electrons} -I det_control.in"
        )

        with out_file.open("w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(  # noqa: S602
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )

    def _optimiser(
        self,
        mol: stk.BuildingBlock,
        prefix: str,
    ) -> stk.BuildingBlock:
        """Optimize `mol`."""
        logger.info("optimising %s", prefix)
        output_dir = self._output_dir / f"{prefix}_opt"
        if output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True)

        init_dir = pathlib.Path.cwd()
        os.chdir(output_dir)

        try:
            xyz = pathlib.Path(f"input_structure_{prefix}.xyz")
            out_file = pathlib.Path(f"optimization_{prefix}.output")
            mol.write(xyz)

            string = f"$gbsa\n gbsagrid={self._solvent_grid}\n"
            with pathlib.Path("det_control.in").open("w") as f:
                f.write(string)

            self._run_xtb(xyz=xyz, out_file=out_file)
            # Check if the optimization is complete.
            output_xyz = pathlib.Path("xtbopt.xyz")

            if not out_file.exists():
                # No simulation has been run.
                msg = "XTB: Optimization failed to start"
                raise RuntimeError(msg)

            # If convergence is achieved, then .xtboptok should exist.
            if not pathlib.Path(".xtboptok").exists():
                # No simulation has been run.
                msg = "XTB: Optimization failed to finish"
                raise RuntimeError(msg)

            if pathlib.Path("NOT_CONVERGED").exists():
                msg = "XTB: Optimization not converged."
                raise RuntimeError(msg)

            # Optimization is complete.
            # Update mol from xtbopt.xyz.
            mol = mol.with_structure_from_file(output_xyz)

        finally:
            os.chdir(init_dir)

        return mol

    def _scanner(self, mol: stk.BuildingBlock, scan_lines: list[str]) -> None:
        """Optimize `mol`."""
        logger.info("scanning over %s torsions", len(scan_lines))
        output_dir = self._output_dir / "scan"
        if output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True)

        init_dir = pathlib.Path.cwd()
        os.chdir(output_dir)

        try:
            xyz = pathlib.Path("input_structure.xyz")
            out_file = pathlib.Path("scan.output")
            mol.write(xyz)

            string = (
                f"$gbsa\n gbsagrid={self._solvent_grid}\n"
                f"$constrain\n  force constant=0.05\n$scan\n"
            )

            for scanline in scan_lines:
                string += f"{scanline}\n"
            string += "$end\n"

            with pathlib.Path("det_control.in").open("w") as f:
                f.write(string)

            self._run_xtb(xyz=xyz, out_file=out_file)

            if not out_file.exists():
                # No simulation has been run.
                msg = "XTB: Optimization failed to start"
                raise RuntimeError(msg)

            # If convergence is achieved, then .xtboptok should exist.
            if not pathlib.Path(".xtboptok").exists():
                # No simulation has been run.
                msg = "XTB: Optimization failed to finish"
                raise RuntimeError(msg)

        finally:
            os.chdir(init_dir)

    def get_matched_torsions(
        self,
        molecule: stk.BuildingBlock,
    ) -> dict[tuple[int], float]:
        """Get the values of the torsions for a molecule."""
        rdkit_molecule = molecule.to_rdkit_mol()
        AllChem.SanitizeMol(rdkit_molecule)  # type: ignore[attr-defined]

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
                    matched_torsions[key] = round(initial_torsion, 2)
        return matched_torsions

    def generate_conformers(
        self,
        molecule: stk.BuildingBlock,
    ) -> Ensemble:
        # Optimise the initial ligand structure.
        molecule = self._optimiser(molecule, prefix="initial")

        ensemble = Ensemble(base_molecule=molecule)
        rdkit_molecule = molecule.to_rdkit_mol()
        AllChem.SanitizeMol(rdkit_molecule)  # type: ignore[attr-defined]

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

        scan_lines = [
            # IDs are +1, start at 1, not 0.
            f"  dihedral: {i[0] + 1},{i[1] + 1},{i[2] + 1},{i[3] + 1},{j[0]};"
            f" {j[0]},{j[-1]},{len(j)}"
            for i, j in matched_torsions.items()
        ]

        self._scanner(molecule, scan_lines=scan_lines)

        scanlog = self._output_dir / "scan" / "xtbscan.log"
        if not scanlog.exists():
            raise RuntimeError

        with scanlog.open("r") as f:
            lines = list(f.readlines())

        num_atoms = molecule.get_num_atoms()

        for i in np.arange(0, len(lines), num_atoms + 2):
            step = lines[i : i + num_atoms + 2]

            if len(step) == 0:
                continue
            if int(step[0].strip()) != num_atoms:
                raise RuntimeError

            position_matrix = np.array(
                [
                    [
                        float(i.strip().split()[1]),
                        float(i.strip().split()[2]),
                        float(i.strip().split()[3]),
                    ]
                    for i in step[2:]
                ]
            )
            conformer = molecule.with_position_matrix(position_matrix)

            angles = self.get_matched_torsions(conformer)

            ensemble.add_conformer(
                conformer=Conformer(
                    molecule=conformer,
                    conformer_id=int(i),
                    source="xtbtorsionscan",
                    score=float(step[1].split()[1]),
                    permutation=angles,
                ),
            )

        return ensemble
