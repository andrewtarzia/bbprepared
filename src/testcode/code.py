#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build the ligand in this project.

Author: Andrew Tarzia

"""

import matplotlib.pyplot as plt
import logging
import pathlib
import sys
import numpy as np
from openbabel import pybel as pb
import copy
import os
import json
import stk
import uuid
import shutil
import stko
from rdkit.Chem import AllChem as rdkit
from rdkit.Chem.rdmolfiles import (
    MolFromMolBlock,
    MolToMolBlock,
    MolToXYZBlock,
)
import matplotlib as mpl

# import shutil
import itertools
from rdkit.Chem import rdMolTransforms
import subprocess as sp


from utilities import (
    update_from_rdkit_conf,
    calculate_N_centroid_N_angle,
    calculate_NN_distance,
    calculate_NN_BCN_angles,
    calculate_NCCN_dihedral,
    get_furthest_pair_FGs,
    AromaticCNCFactory,
)


class XTBFFMD(stko.XTBFF):
    def __init__(
        self,
        xtb_path,
        output_dir=None,
        num_cores=1,
        charge=0,
        unlimited_memory=False,
        temp=298.15,
        time=50.0,
        dump=50.0,
        step=4.0,
        shake=2,
    ):
        """

        temp: 298.15 # in K
        time: 50.0  # in ps
        dump: 50.0  # in fs
        step: 4.0  # in fs
        shake: 2 # 2 for all, 1 for just H, 0 for none

        """

        self._xtb_path = xtb_path
        self._output_dir = output_dir
        self._num_cores = str(num_cores)
        self._charge = str(charge)
        self._unlimited_memory = unlimited_memory
        self._temp = temp
        self._time = time
        self._dump = dump
        self._step = step
        self._shake = shake

    def _run_xtb(self, xyz, out_file):

        # Modify the memory limit.
        if self._unlimited_memory:
            memory = "ulimit -s unlimited ;"
        else:
            memory = ""

        optimization = "--md "

        cmd = (
            f"{memory} {self._xtb_path} {xyz} "
            f"--gfnff "
            f"{optimization} --parallel {self._num_cores} "
            f"--chrg {self._charge} "
            f"-I det_control.in"
        )

        with open(out_file, "w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )

    def _write_detailed_control(self):
        string = (
            f"$md\n"
            f"   temp={self._temp}\n"
            f"   time={self._time}\n"
            f"   dump={self._dump}\n"
            f"   step={self._step}\n"
            f"   velo=false\n"
            f"   nvt =true\n"
            f"   hmass=4\n"
            f"   shake={self._shake}\n"
            f"   sccacc=2.0\n"
            f"$end\n"
        )

        with open("det_control.in", "w") as f:
            f.write(string)

    def _run_md(self, mol):

        xyz = "input_structure_1.xyz"
        out_file = "md_1.output"
        mol.write(xyz)
        self._write_detailed_control()
        self._run_xtb(xyz=xyz, out_file=out_file)
        # Check if the optimization is complete.
        outputchk = "xtbmdok"
        if not os.path.exists(outputchk):
            raise stko.XTBOptimizerError(f"{outputchk} not made.")

        return mol

    def optimize(self, mol):

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        try:
            mol = self._run_md(mol)
        finally:
            os.chdir(init_dir)

        return mol


def split_xyz_file(num_atoms, xyz_file, suffix=".xyz"):
    """
    Splits xyz trajectory file into xyz files.

    """

    with open(xyz_file, "r") as f:
        lines = f.readlines()

    file_strings = []
    string = []
    for line in lines:
        if f" {num_atoms} " in f" {line.strip()} ":
            if len(string) == 0:
                string.append(line)
            else:
                # New block.
                file_strings.append(string)
                string = [line]
        else:
            string.append(line)
    # Add last set.
    file_strings.append(string)

    if suffix not in xyz_file:
        raise ValueError(
            f"{suffix} is not in {xyz_file}. Therefore, the file will "
            'be overwritten! Change "suffix" argument'
        )

    out_files = []
    for i, fs in enumerate(file_strings):
        file_name = xyz_file.replace(suffix, f"_s{i}.xyz")
        with open(file_name, "w") as f:
            for line in fs:
                f.write(line)
        out_files.append(file_name)

    return out_files


def conformer_generation_xtb(
    molecule,
    name,
    conf_data_file,
    calc_dir,
):
    """
    Build a large conformer ensemble with UFF optimisation.

    """

    logging.info(f"building conformer ensemble of {name}")

    xtbmd_output = os.path.join(calc_dir, f"{name}_xtbmd.mol")
    output_dir = os.path.join(calc_dir, f"{name}_xtbmd")

    if not os.path.exists(xtbmd_output):
        logging.info(f"xtbFF MD of {name}")
        xtb_md = XTBFFMD(
            xtb_path="/home/atarzia/miniconda3/envs/beves/bin/xtb",
            output_dir=output_dir,
            # gfn_version=2,
            num_cores=6,
            charge=0,
            # num_unpaired_electrons=0,
            unlimited_memory=True,
            # solvent=solv,
            # solvent_model="alpb",
            # solvent_grid="verytight",
            temp=500,
            time=200,
            dump=500,
            step=1,
            shake=2,
            # step=3,
        )
        molecule = xtb_md.optimize(mol=molecule)
        logging.info(f"xtb MD traj saved: {output_dir}")
        molecule.write(xtbmd_output)
    else:
        logging.info(f"loading: {xtbmd_output}")
        molecule = molecule.with_structure_from_file(xtbmd_output)

    outfiles = split_xyz_file(
        num_atoms=molecule.get_num_atoms(),
        xyz_file=f"{output_dir}/xtb.trj",
        suffix=".trj",
    )

    lig_conf_data = {}
    num_confs = 0
    for cid, xyzfile in enumerate(outfiles):
        logging.info(f"analysing {cid}")

        # Update stk_mol to conformer geometry.
        new_mol = molecule.with_structure_from_file(xyzfile)

        # Need to define the functional groups.
        new_mol = stk.BuildingBlock.init_from_molecule(
            molecule=new_mol,
            functional_groups=[AromaticCNCFactory()],
        )
        # Only get two FGs.
        new_mol = new_mol.with_functional_groups(
            functional_groups=get_furthest_pair_FGs(new_mol),
        )

        # Save unoptimised structure.
        lig_conf_data[f"{cid}_unopt"] = {
            "NcentroidN_angle": calculate_N_centroid_N_angle(new_mol),
            "NCCN_dihedral": calculate_NCCN_dihedral(new_mol),
            "NN_distance": calculate_NN_distance(new_mol),
            "NN_BCN_angles": calculate_NN_BCN_angles(new_mol),
            "UFFEnergy;kj/mol": stko.UFFEnergy(
                ignore_inter_interactions=False
            ).get_energy(new_mol)
            * 4.184,
            "xtbenergy;kj/mol": stko.XTBEnergy(
                xtb_path="/home/atarzia/miniconda3/envs/beves/bin/xtb",
                unlimited_memory=True,
                output_dir=calc_dir / f"{name}_{cid}_unopteyxtbff",
            ).get_energy(new_mol)
            * 2625.5,
        }

        # Save optimised molecule.
        try:
            new_mol = stko.XTBFF(
                xtb_path="/home/atarzia/miniconda3/envs/beves/bin/xtb",
                unlimited_memory=True,
                output_dir=calc_dir / f"{name}_{cid}_xtbffmdunopt",
            ).optimize(mol=new_mol)

            lig_conf_data[f"{cid}_opt"] = {
                "NcentroidN_angle": calculate_N_centroid_N_angle(
                    new_mol
                ),
                "NCCN_dihedral": calculate_NCCN_dihedral(new_mol),
                "NN_distance": calculate_NN_distance(new_mol),
                "NN_BCN_angles": calculate_NN_BCN_angles(new_mol),
                "UFFEnergy;kj/mol": stko.UFFEnergy(
                    ignore_inter_interactions=False
                ).get_energy(new_mol)
                * 4.184,
                "xtbenergy;kj/mol": stko.XTBEnergy(
                    xtb_path="/home/atarzia/miniconda3/envs/beves/bin/xtb",
                    unlimited_memory=True,
                    output_dir=calc_dir / f"{name}_{cid}_xtbffmdopt",
                ).get_energy(new_mol)
                * 2625.5,
            }
            num_confs += 1
        except stko.XTBConvergenceError:
            continue

    logging.info(f"{num_confs} conformers generated for {name}")

    with open(conf_data_file, "w") as f:
        json.dump(lig_conf_data, f)


def conformer_generation_uff(
    molecule,
    name,
    conf_data_file,
    calc_dir,
):
    """
    Build a large conformer ensemble with UFF optimisation.

    """

    logging.info(f"building conformer ensemble of {name}")

    confs = molecule.to_rdkit_mol()
    etkdg = rdkit.srETKDGv3()
    etkdg.randomSeed = 1000
    etkdg.pruneRmsThresh = 0.05
    cids = rdkit.EmbedMultipleConfs(
        mol=confs,
        numConfs=500,
        params=etkdg,
    )

    lig_conf_data = {}
    num_confs = 0
    for cid in cids:
        conf_opt_file_name = calc_dir / f"{name}_c{cid}_cuff.mol"
        # Update stk_mol to conformer geometry.
        new_mol = update_from_rdkit_conf(
            stk_mol=molecule,
            rdk_mol=confs,
            conf_id=cid,
        )
        # Need to define the functional groups.
        new_mol = stk.BuildingBlock.init_from_molecule(
            molecule=new_mol,
            functional_groups=[AromaticCNCFactory()],
        )
        # Only get two FGs.
        new_mol = new_mol.with_functional_groups(
            functional_groups=get_furthest_pair_FGs(new_mol),
        )

        # Save unoptimised structure.
        lig_conf_data[f"{cid}_unopt"] = {
            "NcentroidN_angle": calculate_N_centroid_N_angle(new_mol),
            "NCCN_dihedral": calculate_NCCN_dihedral(new_mol),
            "NN_distance": calculate_NN_distance(new_mol),
            "NN_BCN_angles": calculate_NN_BCN_angles(new_mol),
            "UFFEnergy;kj/mol": stko.UFFEnergy(
                ignore_inter_interactions=False
            ).get_energy(new_mol)
            * 4.184,
            "xtbenergy;kj/mol": stko.XTBEnergy(
                xtb_path="/home/atarzia/miniconda3/envs/beves/bin/xtb",
                unlimited_memory=True,
                output_dir=calc_dir / f"{name}_{cid}_unopteyxtbff",
            ).get_energy(new_mol)
            * 2625.5,
        }

        # Save optimised molecule.
        new_mol = stko.UFF(ignore_inter_interactions=False).optimize(
            mol=new_mol
        )
        new_mol = stko.XTBFF(
            xtb_path="/home/atarzia/miniconda3/envs/beves/bin/xtb",
            unlimited_memory=True,
            output_dir=calc_dir / f"{name}_{cid}_xtbff",
        ).optimize(mol=new_mol)
        new_mol.write(conf_opt_file_name)

        lig_conf_data[f"{cid}_opt"] = {
            "NcentroidN_angle": calculate_N_centroid_N_angle(new_mol),
            "NCCN_dihedral": calculate_NCCN_dihedral(new_mol),
            "NN_distance": calculate_NN_distance(new_mol),
            "NN_BCN_angles": calculate_NN_BCN_angles(new_mol),
            "UFFEnergy;kj/mol": stko.UFFEnergy(
                ignore_inter_interactions=False
            ).get_energy(new_mol)
            * 4.184,
            "xtbenergy;kj/mol": stko.XTBEnergy(
                xtb_path="/home/atarzia/miniconda3/envs/beves/bin/xtb",
                unlimited_memory=True,
                output_dir=calc_dir / f"{name}_{cid}_opteyxtbff",
            ).get_energy(new_mol)
            * 2625.5,
        }
        num_confs += 1
    logging.info(f"{num_confs} conformers generated for {name}")

    with open(conf_data_file, "w") as f:
        json.dump(lig_conf_data, f)


def stk_mol_to_pybel_mol(stk_mol, reperceive_bonds=False):
    if reperceive_bonds:
        return pb.readstring(
            "xyz", MolToXYZBlock(stk_mol.to_rdkit_mol())
        )
    else:
        return pb.readstring(
            "mol", MolToMolBlock(stk_mol.to_rdkit_mol())
        )


def pybel_mol_to_stk_mol(pybel_mol):
    rdkit_mol = MolFromMolBlock(pybel_mol.write("mol"), removeHs=False)
    rdkit.rdmolops.Kekulize(rdkit_mol)
    stk_mol = stk.BuildingBlock.init_from_rdkit_mol(rdkit_mol)
    return stk_mol


def conformer_generation_ob(
    molecule,
    name,
    conf_data_file,
    calc_dir,
):
    pybel_mol = stk_mol_to_pybel_mol(molecule)
    cs = pb.ob.OBConformerSearch()
    # Setup arguments:
    # OBMol, numConformers, numChildren, mutability, convergence.
    cs.Setup(pybel_mol.OBMol, 500, 20, 20, 20)
    cs.Search()
    cs.GetConformers(pybel_mol.OBMol)

    stk_conformers = []
    for i in range(pybel_mol.OBMol.NumConformers()):
        pybel_mol.OBMol.SetConformer(i)
        stk_conformers.append(pybel_mol_to_stk_mol(pybel_mol))

    lig_conf_data = {}
    num_confs = 0
    for cid, stkconf in enumerate(stk_conformers):
        conf_opt_file_name = calc_dir / f"{name}_c{cid}_cob.mol"
        # Update stk_mol to conformer geometry.
        new_mol = molecule.with_position_matrix(
            stkconf.get_position_matrix()
        )
        # Need to define the functional groups.
        new_mol = stk.BuildingBlock.init_from_molecule(
            molecule=new_mol,
            functional_groups=[AromaticCNCFactory()],
        )
        # Only get two FGs.
        new_mol = new_mol.with_functional_groups(
            functional_groups=get_furthest_pair_FGs(new_mol),
        )

        # Save unoptimised structure.
        lig_conf_data[f"{cid}_unopt"] = {
            "NcentroidN_angle": calculate_N_centroid_N_angle(new_mol),
            "NCCN_dihedral": calculate_NCCN_dihedral(new_mol),
            "NN_distance": calculate_NN_distance(new_mol),
            "NN_BCN_angles": calculate_NN_BCN_angles(new_mol),
            "UFFEnergy;kj/mol": stko.UFFEnergy(
                ignore_inter_interactions=False
            ).get_energy(new_mol)
            * 4.184,
        }

        # Save optimised molecule.
        new_mol = stko.UFF(ignore_inter_interactions=False).optimize(
            mol=new_mol
        )
        new_mol.write(conf_opt_file_name)

        lig_conf_data[f"{cid}_opt"] = {
            "NcentroidN_angle": calculate_N_centroid_N_angle(new_mol),
            "NCCN_dihedral": calculate_NCCN_dihedral(new_mol),
            "NN_distance": calculate_NN_distance(new_mol),
            "NN_BCN_angles": calculate_NN_BCN_angles(new_mol),
            "UFFEnergy;kj/mol": stko.UFFEnergy(
                ignore_inter_interactions=False
            ).get_energy(new_mol)
            * 4.184,
        }
        num_confs += 1
    logging.info(f"{num_confs} conformers generated for {name}")

    with open(conf_data_file, "w") as f:
        json.dump(lig_conf_data, f)


def plot_vs_energy(
    uff_results_dict,
    ob_results_dict,
    xtbffmd_results_dict,
    outname,
    yproperty,
    figure_output,
):

    options = (
        (uff_results_dict, "uff", "#c057a1", "UFFEnergy;kj/mol"),
        (ob_results_dict, "ob", "r", "UFFEnergy;kj/mol"),
        (uff_results_dict, "xtb", "gold", "xtbenergy;kj/mol"),
        (
            xtbffmd_results_dict,
            "xtbffmd",
            "skyblue",
            "xtbenergy;kj/mol",
        ),
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    for rd, label, col, elabel in options:
        all_energies = []
        all_values = []
        all_cids = []
        for cid_state in rd:
            cid, state = cid_state.split("_")
            if rd[cid_state][yproperty] is None:
                continue

            energy = rd[cid_state][elabel]

            if yproperty == "NN_BCN_angles":
                value = rd[cid_state][yproperty]["NN_BCN1"]
                all_energies.append(energy)
                all_values.append(value)
                value = rd[cid_state][yproperty]["NN_BCN2"]
                all_energies.append(energy)
                all_values.append(value)
            else:
                value = rd[cid_state][yproperty]
                all_energies.append(energy)
                all_values.append(value)
            all_cids.append(cid)

        all_energies = [i - min(all_energies) for i in all_energies]
        stable_dihedrals = []
        for i, ey in enumerate(all_energies):
            if ey < 5:
                stable_dihedrals.append(all_values[i])

        ax.scatter(
            [i for i in all_values],
            [i for i in all_energies],
            c=col,
            s=100,
            alpha=1.0,
            edgecolors="k",
            label=label,
        )

    ax.axhline(y=5, lw=2, c="k", linestyle="--")
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_xlabel("abs. NCCN dihedral [deg]", fontsize=16)
    ax.set_ylabel("rel. energy [kJmol-1]", fontsize=16)
    # ax.set_xlim(0, 180)
    # ax.set_ylim(0, None)
    ax.set_yscale("log")
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        os.path.join(figure_output, f"{outname}.pdf"),
        dpi=720,
        bbox_inches="tight",
    )
    plt.close()
    return stable_dihedrals


def ligand_smiles():
    return {
        "l1": "C1C=CC(C2C=CC(C3C=CC=NC=3)=CC=2)=CN=1",
        "l2": "C1C=CC(C2C(C)=C(C)C(C3C=NC=CC=3)=C(C)C=2C)=CN=1",
        "l3": "C1C=C(C)C(C2C(C)=C(C)C(C3C(C)=CC=NC=3)=C(C)C=2C)=CN=1",
        "l4": "N1C=C(C2=CC=C(C3C=CC=NC=3)C(C)=C2C)C=CC=1",
        "l5": "N1C(C)=C(C2=CC=C(C3C(C)=CC=NC=3)C(C)=C2C)C=CC=1",
        "l6": "N1C(C)=C(C2=CC=C(C3C(C)=CC=NC=3)C=C2)C=CC=1",
        "l7": "N1C(C)=C(C2=C(C)C(C)=C(C3C(C)=CC=NC=3)C(C)=C2C)C=CC=1",
        "l8": "N1C=C(C2=C(C)C(C)=C(C3C(C)=CC=NC=3)C(C)=C2C)C=CC=1",
        "l9": "CC1=CC(=C(C=C1C2=CN=CC=C2)C)C3=CN=CC=C3",
        "l10": "C(C1C=NC=CC=1C)1C(C)=CC(C2C(C)=CC=NC=2)=C(C)C=1",
        # Below from: 10.1002/chem.202101057 -- numbers match.
        # This one is also from: 10.1039/D3DT00248A -- L1
        "r1": "C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3",
        # This one is also from: 10.1039/D3DT00248A -- L3'
        "r2": "C1=CC(=CC(=C1)C2=CC=NC=C2)C3=CC=NC=C3",
        # This is also from 10.1021/jacs.8b12738 -- L1
        # This one is also from: 10.1039/D3DT00248A -- L1'
        "r3": "C1=CC(=CC(=C1)C#CC2=CN=CC=C2)C#CC3=CN=CC=C3",
        # This one is also from: 10.1039/D3DT00248A -- L3
        "r4": "C1=CC(=CC(=C1)C#CC2=CC=NC=C2)C#CC3=CC=NC=C3",
        "r5": "C1=CN=CC=C1C2=CC=C(S2)C3=CC=NC=C3",
        "r6": "C1=CC(=CC=C1C2=CC=NC=C2)C3=CC=NC=C3",
        "r7": "C1=CC(=CC=C1C#CC2=CC=NC=C2)C#CC3=CC=NC=C3",
        "r8": "C1C=C(OC)C(C2C=CC=C(C3C(OC)=CC=NC=3)C=2)=CN=1",
        # No longer from chem.202...
        # This one is from: 10.1039/D3DT00248A -- L2
        "r9": "C1=CC(=CN=C1)C#CC2=CN=CC=C2",
        # This one is from: 10.1039/D3DT00248A -- L2'
        # Also from 10.1021/jacs.8b12738  -- L2
        "r10": "C1=CC(=CN=C1)C2=CC=C(C=C2)C3=CN=CC=C3",
        # This from 10.1021/jacs.8b12738  -- L3
        "r11": "C1C=C(OC)C(C#CC2C=C(C#CC3C(OC)=CC=NC=3)C=CC=2)=CN=1",
        # This from 10.1021/jacs.8b12738  -- L4
        "r12": "C1C=C(OC)C(C2=CC=C(C=C2)C2C(OC)=CC=NC=2)=CN=1",
    }


def ligand_fakerama_plot(
    name,
    smiles,
    figure_output,
    calculation_output,
):
    conformer_search_dir = calculation_output / f"{name}_search"
    conformer_search_json = conformer_search_dir / "data.json"
    if not os.path.exists(conformer_search_dir):
        os.mkdir(conformer_search_dir)

    if os.path.exists(conformer_search_json):
        with open(conformer_search_json, "r") as f:
            data = json.load(f)

        initial_torsion1 = data["initial_torsion1"]
        initial_torsion2 = data["initial_torsion2"]
        dihedrals = data["dihedrals"]
        bite_angles = data["bite_angles"]
        NN_lengths = data['NN_lengths']
        energies = data["energies"]
        angle1s = data["angle1s"]
        angle2s = data["angle2s"]
        actual_angle1s = data["actual_angle1s"]
        actual_angle2s = data["actual_angle2s"]

    else:
        unopt_mol = stk.BuildingBlock(
            smiles=smiles,
            functional_groups=[AromaticCNCFactory()],
        )

        # Get two binding N atoms.
        binding_N_ids = tuple(
            fg.get_nitrogen().get_id()
            for fg in unopt_mol.get_functional_groups()
        )

        dihedrals = []
        bite_angles = []
        NN_lengths = []
        energies = []
        angle1s = []
        angle2s = []
        actual_angle1s = []
        actual_angle2s = []

        # Optimise the initial ligand structure.
        conf_mol = stko.MMFF(ignore_inter_interactions=False).optimize(
            mol=unopt_mol
        )
        conf_mol.write(conformer_search_dir / "initial.mol")

        rdk_mol = copy.deepcopy(conf_mol.to_rdkit_mol())
        rdkit.SanitizeMol(rdk_mol)
        # rdk_props = rdkit.MMFFGetMoleculeProperties(rdk_mol)

        # Iterate over smarts strings, and their possible torsions.
        target_smarts = {
            "[#7X2]@[#6X3]@[#6X3H0]-!@[#6X3H0]@[#6X3]": (
                0,
                (1, 2, 3, 4),
            ),
            "[#7X2]@[#6X3]@[#6X3]@[#6X3H0]-!@[#6X3H0]@[#6X3]": (
                0,
                (2, 3, 4, 5),
            ),
            "[#7X2]@[#6X3]@[#6X3]~[#6X2]#[#6X2]~[#6X3]@[#6X3]": (
                0,
                (1, 2, 5, 6),
            ),
            "[#7X2]@[#6X3]@[#6X3]@[#6X3]~[#6X2]#[#6X2]~[#6X3]@[#6X3]": (
                0,
                (2, 3, 6, 7),
            ),
        }

        matches_to_test = {}
        smarts_testing = []
        for smarts in target_smarts:
            query = rdkit.MolFromSmarts(smarts)
            matches = rdk_mol.GetSubstructMatches(query)
            logging.info(
                f"num matches for {name},{smarts}: {len(matches)}"
            )
            for match in matches:
                N_atom_id = match[target_smarts[smarts][0]]

                if N_atom_id not in binding_N_ids:
                    continue

                if N_atom_id in matches_to_test:
                    continue

                torsion_atom_ids = tuple(
                    match[j] for j in target_smarts[smarts][1]
                )
                matches_to_test[N_atom_id] = torsion_atom_ids
                smarts_testing.append(smarts)

        if len(matches_to_test) != 2:
            raise ValueError(
                f"{name} does not have 2 matches! ({len(matches_to_test)})"
            )
        logging.info(matches_to_test)
        logging.info(f"smarts:\n{smarts_testing}")

        # Check if the torsion atoms overlap.
        if (
            len(
                set(list(matches_to_test.values())[0]).intersection(
                    list(matches_to_test.values())[1]
                )
            )
            > 1
        ):
            logging.info(
                "there is an overlap in torsions. Only rotating one."
            )
            overlapping = True
        else:
            overlapping = False

        for N_atom_id_1, N_atom_id_2 in itertools.combinations(
            matches_to_test, 2
        ):
            torsion_atom_ids_1 = matches_to_test[N_atom_id_1]
            torsion_atom_ids_2 = matches_to_test[N_atom_id_2]

            # Get initial angles.
            initial_torsion1 = rdMolTransforms.GetDihedralDeg(
                rdk_mol.GetConformer(0),
                torsion_atom_ids_1[0],
                torsion_atom_ids_1[1],
                torsion_atom_ids_1[2],
                torsion_atom_ids_1[3],
            )
            initial_torsion2 = rdMolTransforms.GetDihedralDeg(
                rdk_mol.GetConformer(0),
                torsion_atom_ids_2[0],
                torsion_atom_ids_2[1],
                torsion_atom_ids_2[2],
                torsion_atom_ids_2[3],
            )
            logging.info(
                f"initial torsions: {initial_torsion1}, "
                f"{initial_torsion2}"
            )
            angle_range1 = [
                int(initial_torsion1 + i) for i in range(0, 362, 10)
            ]
            if overlapping:
                angle_range2 = [int(initial_torsion2)]
            else:
                angle_range2 = [
                    int(initial_torsion2 + i) for i in range(0, 362, 10)
                ]

            for angle1, angle2 in itertools.product(
                angle_range1, angle_range2
            ):
                logging.info(
                    f"{name}: doing {angle1} and {angle2} torsions in "
                    f"ranges ({min(angle_range1)} to "
                    f"{max(angle_range1)}; "
                    f"{min(angle_range2)} to {max(angle_range2)})"
                )
                cname = f"c_{angle1}_{int(angle2)}"

                # ff = rdkit.MMFFGetMoleculeForceField(
                # rdk_mol, rdk_props)
                # ff.MMFFAddTorsionConstraint(
                #     torsion_atom_ids[0],
                #     torsion_atom_ids[1],
                #     torsion_atom_ids[2],
                #     torsion_atom_ids[3],
                #     False,
                #     angle - 0.2,
                #     angle + 0.2,
                #     1000.0,
                # )
                # ff.Minimize()

                # new_conf = rdkit.Conformer(rdk_mol.GetNumAtoms())
                # for i in range(rdk_mol.GetNumAtoms()):
                #     new_conf.SetAtomPosition(
                #         i, (
                # rdk_mol.GetConformer(-1).GetAtomPosition(i))
                #     )

                # Set up directory.
                xtb_dir = conformer_search_dir / f"{cname}_xtbopt"
                out_file = xtb_dir / "output.txt"
                if not os.path.exists(out_file):
                    if os.path.exists(xtb_dir):
                        shutil.rmtree(xtb_dir)
                    os.mkdir(xtb_dir)
                    init_dir = os.getcwd()
                    os.chdir(xtb_dir)

                    # Write constraint.inp file.
                    # Xtb indexes start at 1. not 0.
                    atoms1 = (
                        f"{torsion_atom_ids_1[0]+1},"
                        f"{torsion_atom_ids_1[1]+1},"
                        f"{torsion_atom_ids_1[2]+1},"
                        f"{torsion_atom_ids_1[3]+1}"
                    )
                    atoms2 = (
                        f"{torsion_atom_ids_2[0]+1},"
                        f"{torsion_atom_ids_2[1]+1},"
                        f"{torsion_atom_ids_2[2]+1},"
                        f"{torsion_atom_ids_2[3]+1}"
                    )
                    if overlapping:
                        tcstr1 = f"dihedral: {atoms1},{float(angle1)}\n"
                        tcstr2 = ""
                    else:
                        tcstr1 = f"dihedral: {atoms1},{float(angle1)}\n"
                        tcstr2 = f"dihedral: {atoms2},{float(angle2)}\n"
                    with open(xtb_dir / "input.inp", "w") as f:
                        f.write(
                            "$constrain\n"
                            "force constant=1.0\n"
                            f"{tcstr1}"
                            f"{tcstr2}"
                            "$end\n"
                        )

                    # Run xTB opt.
                    conf_mol.write(xtb_dir / "input.xyz")
                    xtb_path = (
                        "/home/atarzia/miniconda3/envs/simple_het/bin/"
                        "xtb"
                    )
                    cmd = (
                        "ulimit -s unlimited ; "
                        f"{xtb_path} "
                        "input.xyz --gfnff "
                        "--opt normal --parallel 4 "
                        "--chrg 0 "
                        "--input input.inp"
                    )
                    with open(out_file, "w") as f:
                        # Note that sp.call will hold the program until
                        # completion of the calculation.
                        sp.call(
                            cmd,
                            stdin=sp.PIPE,
                            stdout=f,
                            stderr=sp.PIPE,
                            # Shell is required for complex arguments.
                            shell=True,
                        )
                    os.chdir(init_dir)

                # Extract structure and energy.
                # conf_mol = update_from_rdkit_conf(
                #     stk_mol=opt_mol,
                #     rdk_mol=rdk_mol,
                #     conf_id=-1,
                # )

                try:
                    conf_mol = conf_mol.with_structure_from_file(
                        xtb_dir / "xtbopt.xyz"
                    )
                    with open(out_file, "r") as f:
                        for line in f.readlines():
                            if "TOTAL ENERGY" in line:
                                XtbEnergy = float(
                                    line.strip().split()[3]
                                )

                except FileNotFoundError:
                    logging.info(
                        f"{xtb_dir} did not finish - ignoring this pair"
                    )

                    # This jumps back quite a bit in both angles, to
                    # avoid just being at a clash point again.
                    if overlapping:
                        previous_cname = f"c_{angle1-20}_{int(angle2)}"
                    else:
                        previous_cname = f"c_{angle1-20}_{angle2-20}"

                    previous_xtb_dir = (
                        conformer_search_dir
                        / f"{previous_cname}_xtbopt"
                    )
                    if not os.path.exists(
                        previous_xtb_dir / "xtbopt.xyz"
                    ):
                        if overlapping:
                            previous_cname = (
                                f"c_{angle1-40}_{int(angle2)}"
                            )
                        else:
                            previous_cname = (
                                f"c_{angle1-40}_{angle2-40}"
                            )
                        previous_xtb_dir = (
                            conformer_search_dir
                            / f"{previous_cname}_xtbopt"
                        )

                    if not os.path.exists(
                        previous_xtb_dir / "xtbopt.xyz"
                    ):
                        raise ValueError(
                            "Had to go back many times, and none found."
                            " This is broken!"
                        )
                    conf_mol = conf_mol.with_structure_from_file(
                        previous_xtb_dir / "xtbopt.xyz"
                    )
                    continue

                conf_mol.write(
                    conformer_search_dir / f"conf_{cname}.mol"
                )
                NCCN_dihedral = calculate_NCCN_dihedral(conf_mol)
                NN_BCN_angles = calculate_NN_BCN_angles(conf_mol)
                NN_length = calculate_NN_distance(conf_mol)
                # MMFFEnergy = (
                #     stko.MMFFEnergy(
                #         ignore_inter_interactions=False
                #     ).get_energy(conf_mol)
                #     * 4.184
                # )

                angle1s.append(angle1)
                angle2s.append(angle2)

                # New rdkit conformer and calculate actual torsions.
                temp_mol = copy.deepcopy(conf_mol.to_rdkit_mol())
                rdkit.SanitizeMol(temp_mol)
                actual_angle1s.append(
                    rdMolTransforms.GetDihedralDeg(
                        temp_mol.GetConformer(0),
                        torsion_atom_ids_1[0],
                        torsion_atom_ids_1[1],
                        torsion_atom_ids_1[2],
                        torsion_atom_ids_1[3],
                    )
                )
                actual_angle2s.append(
                    rdMolTransforms.GetDihedralDeg(
                        temp_mol.GetConformer(0),
                        torsion_atom_ids_2[0],
                        torsion_atom_ids_2[1],
                        torsion_atom_ids_2[2],
                        torsion_atom_ids_2[3],
                    )
                )
                bite_angle = (
                    (180 - NN_BCN_angles["NN_BCN1"])
                    - 90
                    + (180 - NN_BCN_angles["NN_BCN2"])
                    - 90
                )

                dihedrals.append(NCCN_dihedral)
                bite_angles.append(bite_angle)
                energies.append(XtbEnergy)
                NN_lengths.append(NN_length)

        data = {}
        data["initial_torsion1"] = initial_torsion1
        data["initial_torsion2"] = initial_torsion2
        data["dihedrals"] = dihedrals
        data["bite_angles"] = bite_angles
        data['NN_lengths'] = NN_lengths
        data["energies"] = energies
        data["angle1s"] = angle1s
        data["angle2s"] = angle2s
        data["actual_angle1s"] = actual_angle1s
        data["actual_angle2s"] = actual_angle2s
        with open(conformer_search_json, "w") as f:
            json.dump(data, f, indent=4)

    relative_energies = [(i - min(energies)) * 2625.5 for i in energies]

    # Plot bite-angle vs. NCCN torsion angle.
    fig, axs = plt.subplots(ncols=3, nrows=3, figsize=(16, 10))
    ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8 = axs.flatten()
    vmax = 10

    ax0.scatter(
        [i - initial_torsion1 for i in angle1s],
        [i - initial_torsion2 for i in angle2s],
        c=relative_energies,
        vmin=0,
        vmax=vmax,
        s=50,
        alpha=1.0,
        edgecolors="k",
        cmap="Blues_r",
    )
    ax0.tick_params(axis="both", which="major", labelsize=16)
    ax0.set_xlabel("delta target torsion 1 [deg]", fontsize=16)
    ax0.set_ylabel("delta target torsion 2 [deg]", fontsize=16)
    cmap = mpl.cm.Blues_r
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax0,
        orientation="vertical",
    )
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label("rel. energy [kJ/mol]", fontsize=16)

    ax1.scatter(
        [i - initial_torsion1 for i in angle1s],
        [i - initial_torsion2 for i in angle2s],
        c=dihedrals,
        vmin=-180,
        vmax=180,
        s=50,
        alpha=1.0,
        edgecolors="k",
        cmap="PiYG",
    )
    ax1.tick_params(axis="both", which="major", labelsize=16)
    ax1.set_xlabel("delta target torsion 1 [deg]", fontsize=16)
    ax1.set_ylabel("delta target torsion 2 [deg]", fontsize=16)
    cmap = mpl.cm.PiYG
    norm = mpl.colors.Normalize(vmin=-180, vmax=180)
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax1,
        orientation="vertical",
    )
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label("NCCN dihedral [deg]", fontsize=16)

    ax2.scatter(
        dihedrals,
        relative_energies,
        c="gray",
        s=50,
        alpha=1.0,
        edgecolors="k",
    )
    ax2.tick_params(axis="both", which="major", labelsize=16)
    ax2.set_xlabel("NCCN dihedral [deg]", fontsize=16)
    ax2.set_ylabel("rel. energy [kJ/mol]", fontsize=16)
    ax2.axvline(x=0, c="r", lw=2, ls="--")
    ax2.axvline(x=60, c="r", lw=2, ls="--")
    ax2.axvline(x=-60, c="r", lw=2, ls="--")
    ax2.set_xlim(-180, 180)
    ax2.set_ylim(0, 10)

    ax3.scatter(
        actual_angle1s,
        actual_angle2s,
        c=relative_energies,
        vmin=0,
        vmax=vmax,
        s=50,
        alpha=1.0,
        edgecolors="k",
        cmap="Blues_r",
    )
    ax3.tick_params(axis="both", which="major", labelsize=16)
    ax3.set_xlabel("measured torsion 1 [deg]", fontsize=16)
    ax3.set_ylabel("measured torsion 2 [deg]", fontsize=16)
    cmap = mpl.cm.Blues_r
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax3,
        orientation="vertical",
    )
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label("rel. energy [kJ/mol]", fontsize=16)

    energies_near_0 = []
    energies_near_60 = []
    for energy, dih in zip(relative_energies, dihedrals):
        if dih < 10 and dih > -10:
            energies_near_0.append(energy)
        if dih < 70 and dih > 50:
            energies_near_60.append(energy)
        if dih > -70 and dih < -50:
            energies_near_60.append(energy)

    ax4.hist(
        energies_near_0,
        bins=np.arange(0, 30.1, 1.0),
        alpha=0.9,
        density=True,
        # edgecolors="k",
        label="near 0",
    )
    ax4.hist(
        energies_near_60,
        bins=np.arange(0, 30.1, 1.0),
        alpha=0.9,
        density=True,
        # edgecolors="k",
        label="near |60|",
    )
    ax4.tick_params(axis="both", which="major", labelsize=16)
    ax4.set_ylabel("frequency", fontsize=16)
    ax4.set_xlabel("rel. energy [kJ/mol]", fontsize=16)
    ax4.legend(fontsize=16)

    ax5.scatter(
        actual_angle1s,
        dihedrals,
        c=relative_energies,
        vmin=0,
        vmax=vmax,
        s=50,
        alpha=1.0,
        edgecolors="k",
        cmap="Blues_r",
    )
    ax5.tick_params(axis="both", which="major", labelsize=16)
    ax5.set_xlabel("measured torsion 1 [deg]", fontsize=16)
    ax5.set_ylabel("NCCN dihedral [deg]", fontsize=16)
    ax5.axhline(y=0, c="r", lw=2, ls="--")
    ax5.axhline(y=60, c="r", lw=2, ls="--")
    ax5.axhline(y=-60, c="r", lw=2, ls="--")
    ax5.set_ylim(-180, 180)
    cmap = mpl.cm.Blues_r
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax5,
        orientation="vertical",
    )
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label("rel. energy [kJ/mol]", fontsize=16)

    ax6.scatter(
        actual_angle2s,
        dihedrals,
        c=relative_energies,
        vmin=0,
        vmax=vmax,
        s=50,
        alpha=1.0,
        edgecolors="k",
        cmap="Blues_r",
    )
    ax6.tick_params(axis="both", which="major", labelsize=16)
    ax6.set_xlabel("measured torsion 2 [deg]", fontsize=16)
    ax6.set_ylabel("NCCN dihedral [deg]", fontsize=16)
    ax6.axhline(y=0, c="r", lw=2, ls="--")
    ax6.axhline(y=60, c="r", lw=2, ls="--")
    ax6.axhline(y=-60, c="r", lw=2, ls="--")
    ax6.set_ylim(-180, 180)
    cmap = mpl.cm.Blues_r
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax6,
        orientation="vertical",
    )
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label("rel. energy [kJ/mol]", fontsize=16)

    ax7.scatter(
        [i for i, j in zip(bite_angles, relative_energies) if j < vmax],
        [i for i, j in zip(dihedrals, relative_energies) if j < vmax],
        c=[i for i in relative_energies if i < vmax],
        vmin=0,
        vmax=vmax,
        s=50,
        alpha=0.5,
        edgecolors="k",
        cmap="Blues_r",
    )
    ax7.tick_params(axis="both", which="major", labelsize=16)
    ax7.set_xlabel("bite angle [deg]", fontsize=16)
    ax7.set_ylabel("NCCN dihedral [deg]", fontsize=16)
    ax7.set_xlim(0, 180)
    ax7.set_ylim(-180, 180)
    cmap = mpl.cm.Blues_r
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax7,
        orientation="vertical",
    )
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label("rel. energy [kJ/mol]", fontsize=16)

    ax8.scatter(
        [i for i, j in zip(NN_lengths, relative_energies) if j < vmax],
        [i for i, j in zip(dihedrals, relative_energies) if j < vmax],
        c=[i for i in relative_energies if i < vmax],
        vmin=0,
        vmax=vmax,
        s=50,
        alpha=0.5,
        edgecolors="k",
        cmap="Blues_r",
    )
    ax8.tick_params(axis="both", which="major", labelsize=16)
    ax8.set_xlabel("NN distance [A]", fontsize=16)
    ax8.set_ylabel("NCCN dihedral [deg]", fontsize=16)
    ax8.set_xlim(5, 20)
    ax8.set_ylim(-180, 180)
    cmap = mpl.cm.Blues_r
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax8,
        orientation="vertical",
    )
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label("rel. energy [kJ/mol]", fontsize=16)

    fig.tight_layout()
    fig.savefig(
        os.path.join(figure_output, f"{name}_fakerama.pdf"),
        dpi=720,
        bbox_inches="tight",
    )
    plt.close()


def main():
    if not len(sys.argv) == 1:
        logging.info(f"Usage: {__file__}\n" "   Expected 0 arguments:")
        sys.exit()
    else:
        pass

    working_dir = pathlib.Path(
        "/home/atarzia/workingspace/cg-4p82-finder/initial_tests"
    )

    figure_output = working_dir / "figures"
    calculation_output = working_dir / "calculations"

    lsmiles = ligand_smiles()
    stable_states = {}
    for lig in lsmiles:
        logging.info(f"running {lig}")
        unopt_file = calculation_output / f"{lig}_unopt.mol"
        unopt_xyz_file = calculation_output / f"{lig}_unopt.xyz"

        ligand_fakerama_plot(
            name=lig,
            smiles=lsmiles[lig],
            figure_output=figure_output,
            calculation_output=calculation_output,
        )

        continue

        confuff_data_file = (
            calculation_output / f"{lig}_conf_uff_data.json"
        )
        confob_data_file = (
            calculation_output / f"{lig}_conf_ob_data.json"
        )
        confxtb_data_file = (
            calculation_output / f"{lig}_conf_xtbmd_data.json"
        )

        unopt_mol = stk.BuildingBlock(
            smiles=lsmiles[lig],
            functional_groups=[AromaticCNCFactory()],
        )
        unopt_mol.write(unopt_file)
        unopt_mol.write(unopt_xyz_file)

        if not os.path.exists(confxtb_data_file):
            conformer_generation_xtb(
                molecule=unopt_mol,
                name=lig,
                conf_data_file=confxtb_data_file,
                calc_dir=calculation_output,
            )
        with open(confxtb_data_file, "r") as f:
            xtbffmd_property_dict = json.load(f)

        if not os.path.exists(confob_data_file):
            conformer_generation_ob(
                molecule=unopt_mol,
                name=lig,
                conf_data_file=confob_data_file,
                calc_dir=calculation_output,
            )
        with open(confob_data_file, "r") as f:
            ob_property_dict = json.load(f)

        if not os.path.exists(confuff_data_file):
            conformer_generation_uff(
                molecule=unopt_mol,
                name=lig,
                conf_data_file=confuff_data_file,
                calc_dir=calculation_output,
            )
        with open(confuff_data_file, "r") as f:
            uff_property_dict = json.load(f)

        _ = plot_vs_energy(
            uff_results_dict=uff_property_dict,
            ob_results_dict=ob_property_dict,
            xtbffmd_results_dict=xtbffmd_property_dict,
            outname=f"{lig}_bite",
            yproperty="NN_BCN_angles",
            figure_output=figure_output,
        )
        stable_dihedrals = plot_vs_energy(
            uff_results_dict=uff_property_dict,
            ob_results_dict=ob_property_dict,
            xtbffmd_results_dict=xtbffmd_property_dict,
            outname=lig,
            yproperty="NCCN_dihedral",
            figure_output=figure_output,
        )
        stable_states[lig] = stable_dihedrals

    fig, ax = plt.subplots(figsize=(8, 5))
    for i, lig in enumerate(stable_states):
        stable_dihedrals = stable_states[lig]
        ax.scatter(
            [dihedral for dihedral in stable_dihedrals],
            [i for dihedral in stable_dihedrals],
            s=100,
            alpha=1.0,
            edgecolors="k",
        )

    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_xlabel("stable NCCN dihedrals [deg]", fontsize=16)
    ax.set_ylabel("ligand", fontsize=16)
    # ax.set_xlim(0, 180)
    # ax.set_ylim(0, None)

    fig.tight_layout()
    fig.savefig(
        os.path.join(figure_output, "stablestates.pdf"),
        dpi=720,
        bbox_inches="tight",
    )
    plt.close()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
