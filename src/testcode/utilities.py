#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for utility functions.

Author: Andrew Tarzia

"""

import os
import stk
import numpy as np
import networkx as nx
from itertools import combinations
from scipy.spatial.distance import euclidean
import pymatgen.core as pmg
from pymatgen.analysis.local_env import (
    LocalStructOrderParams,
)


class AromaticCNCFactory(stk.FunctionalGroupFactory):
    """
    A subclass of stk.SmartsFunctionalGroupFactory.

    """

    def __init__(self, bonders=(1,), deleters=()):
        """
        Initialise :class:`.AromaticCNCFactory`.

        """

        self._bonders = bonders
        self._deleters = deleters

    def get_functional_groups(self, molecule):
        generic_functional_groups = stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=self._bonders,
            deleters=self._deleters,
        ).get_functional_groups(molecule)
        for fg in generic_functional_groups:
            atom_ids = (i.get_id() for i in fg.get_atoms())
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield AromaticCNC(
                carbon1=atoms[0],
                nitrogen=atoms[1],
                carbon2=atoms[2],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )


class AromaticCNC(stk.GenericFunctionalGroup):
    """
    Represents an N atom in pyridine functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon][nitrogen][carbon]``.

    """

    def __init__(
        self,
        carbon1,
        nitrogen,
        carbon2,
        bonders,
        deleters,
    ):
        """
        Initialize a :class:`.AromaticCNC` instance.

        Parameters
        ----------
        carbon1 : :class:`.C`
            The first carbon atom.

        nitrogen : :class:`.N`
            The nitrogen atom.

        carbon2 : :class:`.C`
            The second carbon atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._carbon1 = carbon1
        self._nitrogen = nitrogen
        self._carbon2 = carbon2
        atoms = (carbon1, nitrogen, carbon2)
        super().__init__(atoms, bonders, deleters)

    def get_carbon1(self):
        return self._carbon1

    def get_carbon2(self):
        return self._carbon2

    def get_nitrogen(self):
        return self._nitrogen

    def clone(self):
        clone = super().clone()
        clone._carbon1 = self._carbon1
        clone._nitrogen = self._nitrogen
        clone._carbon2 = self._carbon2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon1 = atom_map.get(
            self._carbon1.get_id(),
            self._carbon1,
        )
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        return clone

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            f"{self._carbon1}, {self._nitrogen}, {self._carbon2}, "
            f"bonders={self._bonders})"
        )


def unit_vector(vector):
    """
    Returns the unit vector of the vector.

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249

    """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2, normal=None):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'::

        >>> angle_between((1, 0, 0), (0, 1, 0))
        1.5707963267948966
        >>> angle_between((1, 0, 0), (1, 0, 0))
        0.0
        >>> angle_between((1, 0, 0), (-1, 0, 0))
        3.141592653589793

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249

    If normal is given, the angle polarity is determined using the
    cross product of the two vectors.

    """

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    if normal is not None:
        # Get normal vector and cross product to determine sign.
        cross = np.cross(v1_u, v2_u)
        if np.dot(normal, cross) < 0:
            angle = -angle
    return angle


def convert_stk_to_pymatgen(stk_mol):
    """
    Convert stk.Molecule to pymatgen.Molecule.

    Parameters
    ----------
    stk_mol : :class:`stk.Molecule`
        Stk molecule to convert.

    Returns
    -------
    pmg_mol : :class:`pymatgen.Molecule`
        Corresponding pymatgen Molecule.

    """
    stk_mol.write("temp.xyz")
    pmg_mol = pmg.Molecule.from_file("temp.xyz")
    os.system("rm temp.xyz")

    return pmg_mol


def calculate_sites_order_values(
    molecule, site_idxs, target_species_type=None, neigh_idxs=None
):
    """
    Calculate order parameters around metal centres.

    Parameters
    ----------
    molecule : :class:`pmg.Molecule` or :class:`pmg.Structure`
        Pymatgen (pmg) molecule/structure to analyse.

    site_idxs : :class:`list` of :class:`int`
        Atom ids of sites to calculate OP of.

    target_species_type : :class:`str`
        Target neighbour element to use in OP calculation.
        Defaults to :class:`NoneType` if no target species is known.

    neigh_idxs : :class:`list` of :class:`list` of :class:`int`
        Neighbours of each atom in site_idx. Ordering is important.
        Defaults to :class:`NoneType` for when using
        :class:`pmg.Structure` - i.e. a structure with a lattice.

    Returns
    -------
    results : :class:`dict`
        Dictionary of format
        site_idx: dict of order parameters
        {
            `oct`: :class:`float`,
            `sq_plan`: :class:`float`,
            `q2`: :class:`float`,
            `q4`: :class:`float`,
            `q6`: :class:`float`
        }.

    """

    results = {}

    if target_species_type is None:
        targ_species = None
    else:
        targ_species = pmg.Species(target_species_type)

    # Define local order parameters class based on desired types.
    types = [
        "sq_plan",  # Square planar envs.
    ]
    loc_ops = LocalStructOrderParams(
        types=types,
    )
    if neigh_idxs is None:
        for site in site_idxs:
            site_results = loc_ops.get_order_parameters(
                structure=molecule, n=site, target_spec=[targ_species]
            )
            results[site] = {i: j for i, j in zip(types, site_results)}
    else:
        for site, neigh in zip(site_idxs, neigh_idxs):
            site_results = loc_ops.get_order_parameters(
                structure=molecule,
                n=site,
                indices_neighs=neigh,
                target_spec=targ_species,
            )
            results[site] = {i: j for i, j in zip(types, site_results)}

    return results


def get_order_values(mol, metal, per_site=False):
    """
    Calculate order parameters around metal centres.

    Parameters
    ----------
    mol : :class:`stk.ConstructedMolecule`
        stk molecule to analyse.

    metal : :class:`int`
        Element number of metal atom.

    per_site : :class:`bool`
        Defaults to False. True if the OPs for each site are desired.

    Returns
    -------
    results : :class:`dict`
        Dictionary of order parameter max/mins/averages if `per_site`
        is False.

    """

    pmg_mol = convert_stk_to_pymatgen(stk_mol=mol)
    # Get sites of interest and their neighbours.
    sites = []
    neighs = []
    for atom in mol.get_atoms():
        if atom.get_atomic_number() == metal:
            sites.append(atom.get_id())
            bonds = [
                i
                for i in mol.get_bonds()
                if i.get_atom1().get_id() == atom.get_id()
                or i.get_atom2().get_id() == atom.get_id()
            ]
            a_neigh = []
            for b in bonds:
                if b.get_atom1().get_id() == atom.get_id():
                    a_neigh.append(b.get_atom2().get_id())
                elif b.get_atom2().get_id() == atom.get_id():
                    a_neigh.append(b.get_atom1().get_id())
            neighs.append(a_neigh)

    order_values = calculate_sites_order_values(
        molecule=pmg_mol,
        site_idxs=sites,
        neigh_idxs=neighs,
    )

    if per_site:
        results = order_values
        return results
    else:
        # Get max, mins and averages of all OPs for the whole molecule.
        OPs = [order_values[i].keys() for i in order_values][0]
        OP_lists = {}
        for OP in OPs:
            OP_lists[OP] = [order_values[i][OP] for i in order_values]

        results = {
            # OP: (min, max, avg)
            i: {
                "min": min(OP_lists[i]),
                "max": max(OP_lists[i]),
                "avg": np.average(OP_lists[i]),
            }
            for i in OP_lists
        }

        return results


def decompose_cage(cage, metal_atom_nos):

    # Produce a graph from the cage that does not include metals.
    cage_g = nx.Graph()
    atom_ids_in_G = set()
    for atom in cage.get_atoms():
        if atom.get_atomic_number() in metal_atom_nos:
            continue
        cage_g.add_node(atom)
        atom_ids_in_G.add(atom.get_id())

    # Add edges.
    for bond in cage.get_bonds():
        a1id = bond.get_atom1().get_id()
        a2id = bond.get_atom2().get_id()
        if a1id in atom_ids_in_G and a2id in atom_ids_in_G:
            cage_g.add_edge(bond.get_atom1(), bond.get_atom2())

    # Get disconnected subgraphs as molecules.
    # Sort and sort atom ids to ensure molecules are read by RDKIT
    # correctly.
    connected_graphs = [
        sorted(subgraph, key=lambda a: a.get_id())
        for subgraph in sorted(nx.connected_components(cage_g))
    ]
    return connected_graphs


def get_organic_linkers(
    cage,
    metal_atom_nos,
    calc_dir,
    file_prefix=None,
):
    """
    Extract a list of organic linker .Molecules from a cage.

    Parameters
    ----------
    cage : :class:`stk.Molecule`
        Molecule to get the organic linkers from.

    metal_atom_nos : :class:`iterable` of :class:`int`
        The atomic number of metal atoms to remove from structure.

    file_prefix : :class:`str`, optional
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    Returns
    -------
    org_lig : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    """

    connected_graphs = decompose_cage(cage, metal_atom_nos)
    org_lig = {}
    smiles_keys = {}
    for i, cg in enumerate(connected_graphs):
        # Get atoms from nodes.
        atoms = list(cg)
        atom_ids = [i.get_id() for i in atoms]
        cage.write(
            "temporary_linker.mol",
            atom_ids=atom_ids,
        )
        temporary_linker = stk.BuildingBlock.init_from_file(
            "temporary_linker.mol"
        ).with_canonical_atom_ordering()
        os.system("rm temporary_linker.mol")
        smiles_key = stk.Smiles().get_key(temporary_linker)
        if smiles_key not in smiles_keys:
            smiles_keys[smiles_key] = len(smiles_keys.values()) + 1
        idx = smiles_keys[smiles_key]
        sgt = str(len(atoms))
        # Write to mol file.
        if file_prefix is None:
            filename_ = f"organic_linker_s{sgt}_{idx}_{i}.mol"
        else:
            filename_ = f"{file_prefix}{sgt}_{idx}_{i}.mol"

        final_path = os.path.join(calc_dir, filename_)
        if os.path.exists(final_path):
            org_lig[filename_] = stk.BuildingBlock.init_from_file(
                path=final_path,
            )
        else:
            org_lig[filename_] = temporary_linker
            # Rewrite to fix atom ids.
            org_lig[filename_].write(os.path.join(calc_dir, filename_))
            org_lig[filename_] = stk.BuildingBlock.init_from_file(
                path=final_path,
            )

    return org_lig, smiles_keys


def get_organic_linker_atoms(cage, metal_atom_nos):

    connected_graphs = decompose_cage(cage, metal_atom_nos)
    for i, cg in enumerate(connected_graphs):
        # Get atoms from nodes.
        atoms = list(cg)
        yield atoms


def read_xtb_energy(name, calc_dir):
    raise SystemExit("handle solvent")
    output_file = os.path.join(calc_dir, f"{name}_xtb.ey")
    if os.path.exists(output_file):
        with open(output_file, "r") as f:
            lines = f.readlines()
        for line in lines:
            energy = float(line.rstrip())
            break
    else:
        raise FileNotFoundError(f"{output_file} not found.")

    # In a.u.
    return energy


def get_dft_energy(name, txt_file):
    with open(txt_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        if name in line:
            number = line.strip().split()[1]
            if number == "None":
                return None
            else:
                energy = float(number)
                # kJmol-1
                return energy


def get_stoichiometry(topology_string):
    return {
        "m2": (4, 0),
        "m3": (6, 0),
        "m4": (8, 0),
        "trans": (4, 0),
        "cis": (2, 2),
    }[topology_string]


def name_parser(name):
    splits = name.split("_")
    if len(splits) == 2:
        topo, ligand1 = splits
        ligand2 = None
    elif len(splits) == 3:
        topo, ligand1, ligand2 = splits

    return topo, ligand1, ligand2


def update_from_rdkit_conf(stk_mol, rdk_mol, conf_id):
    """
    Update the structure to match `conf_id` of `mol`.

    Parameters
    ----------
    struct : :class:`stk.Molecule`
        The molecule whoce coordinates are to be updated.

    mol : :class:`rdkit.Mol`
        The :mod:`rdkit` molecule to use for the structure update.

    conf_id : :class:`int`
        The conformer ID of the `mol` to update from.

    Returns
    -------
    :class:`.Molecule`
        The molecule.

    """

    pos_mat = rdk_mol.GetConformer(id=conf_id).GetPositions()
    return stk_mol.with_position_matrix(pos_mat)


def calculate_N_centroid_N_angle(bb):
    """
    Calculate the N-centroid-N angle of a ditopic building block.

    This function will not work for cages built from FGs other than
    metals + AromaticCNC and metals + AromaticCNN.

    Parameters
    ----------
    bb : :class:`stk.BuildingBlock`
        stk molecule to analyse.

    Returns
    -------
    angle : :class:`float`
        Angle between two bonding vectors of molecule.

    """

    fg_counts = 0
    fg_positions = []
    for fg in bb.get_functional_groups():
        if isinstance(fg, AromaticCNC):
            fg_counts += 1
            # Get geometrical properties of the FG.
            # Get N position - deleter.
            (N_position,) = bb.get_atomic_positions(
                atom_ids=fg.get_nitrogen().get_id()
            )
            fg_positions.append(N_position)

    if fg_counts != 2:
        raise ValueError(
            f"{bb} does not have 2 AromaticCNC or AromaticCNN "
            "functional groups."
        )

    # Get building block centroid.
    centroid_position = bb.get_centroid()

    # Get vectors.
    fg_vectors = [i - centroid_position for i in fg_positions]

    # Calculate the angle between the two vectors.
    angle = np.degrees(angle_between(*fg_vectors))
    return angle


def calculate_NN_distance(bb):
    """
    Calculate the N-N distance of ditopic building block.

    This function will not work for cages built from FGs other than
    metals + AromaticCNC.

    Parameters
    ----------
    bb : :class:`stk.BuildingBlock`
        stk molecule to analyse.

    Returns
    -------
    NN_distance : :class:`float`
        Distance(s) between [angstrom] N atoms in functional groups.

    """

    fg_counts = 0
    N_positions = []
    for fg in bb.get_functional_groups():
        if isinstance(fg, AromaticCNC):
            fg_counts += 1
            # Get geometrical properties of the FG.
            # Get N position - deleter.
            (N_position,) = bb.get_atomic_positions(
                atom_ids=fg.get_nitrogen().get_id()
            )
            N_positions.append(N_position)

    if fg_counts != 2:
        raise ValueError(
            f"{bb} does not have 2 AromaticCNC functional groups."
        )

    NN_distance = np.linalg.norm(N_positions[0] - N_positions[1])
    return NN_distance


def calculate_NN_BCN_angles(bb):
    fg_counts = 0
    N_positions = []
    N_C_vectors = []
    for fg in bb.get_functional_groups():
        if isinstance(fg, AromaticCNC):
            fg_counts += 1
            # Get geometrical properties of the FG.
            # Get N position - deleter.
            (N_position,) = bb.get_atomic_positions(
                atom_ids=fg.get_nitrogen().get_id()
            )
            N_positions.append(N_position)
            C_atom_ids = (
                fg.get_carbon1().get_id(),
                fg.get_carbon2().get_id(),
            )
            C_centroid = bb.get_centroid(atom_ids=C_atom_ids)
            N_C_vector = N_position - C_centroid
            N_C_vectors.append(N_C_vector)

    if fg_counts != 2:
        raise ValueError(
            f"{bb} does not have 2 AromaticCNC functional groups."
        )

    NN_vector = N_positions[1] - N_positions[0]

    NN_BCN_1 = np.degrees(
        angle_between(
            N_C_vectors[0],
            -NN_vector,
        )
    )
    NN_BCN_2 = np.degrees(
        angle_between(
            N_C_vectors[1],
            NN_vector,
        )
    )

    return {"NN_BCN1": NN_BCN_1, "NN_BCN2": NN_BCN_2}


def get_furthest_pair_FGs(stk_mol):
    """
    Returns the pair of functional groups that are furthest apart.

    """

    if stk_mol.get_num_functional_groups() == 2:
        return tuple(i for i in stk_mol.get_functional_groups())
    elif stk_mol.get_num_functional_groups() < 2:
        raise ValueError(f"{stk_mol} does not have at least 2 FGs")

    fg_centroids = [
        (fg, stk_mol.get_centroid(atom_ids=fg.get_placer_ids()))
        for fg in stk_mol.get_functional_groups()
    ]

    fg_dists = sorted(
        [
            (i[0], j[0], euclidean(i[1], j[1]))
            for i, j in combinations(fg_centroids, 2)
        ],
        key=lambda x: x[2],
        reverse=True,
    )

    return (fg_dists[0][0], fg_dists[0][1])


def get_dihedral(pt1, pt2, pt3, pt4):
    """
    Calculate the dihedral between four points.

    Uses Praxeolitic formula --> 1 sqrt, 1 cross product
    Output in range (-pi to pi).

    From: https://stackoverflow.com/questions/20305272/
    dihedral-torsion-angle-from-four-points-in-cartesian-
    coordinates-in-python
    (new_dihedral(p))
    """

    p0 = np.asarray(pt1)
    p1 = np.asarray(pt2)
    p2 = np.asarray(pt3)
    p3 = np.asarray(pt4)

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def calculate_NCCN_dihedral(bb):
    fg_counts = 0
    N_positions = []
    C_centroids = []
    for fg in bb.get_functional_groups():
        if isinstance(fg, AromaticCNC):
            fg_counts += 1
            # Get geometrical properties of the FG.
            # Get N position - deleter.
            (N_position,) = bb.get_atomic_positions(
                atom_ids=fg.get_nitrogen().get_id()
            )
            N_positions.append(N_position)
            C_atom_ids = (
                fg.get_carbon1().get_id(),
                fg.get_carbon2().get_id(),
            )
            C_centroid = bb.get_centroid(atom_ids=C_atom_ids)
            C_centroids.append(C_centroid)

    if fg_counts != 2:
        raise ValueError(
            f"{bb} does not have 2 AromaticCNC functional groups."
        )

    return get_dihedral(
        pt1=N_positions[0],
        pt2=C_centroids[0],
        pt3=C_centroids[1],
        pt4=N_positions[1],
    )


def calculate_helicities(molecule, name, calc_dir):

    pos_mat = molecule.get_position_matrix()
    n_bonds = {}
    for bond in molecule.get_bonds():
        atom1 = bond.get_atom1()
        atom2 = bond.get_atom2()
        if (
            atom1.get_atomic_number() == 46
            and atom2.get_atomic_number() == 7
        ):
            if atom2.get_id() not in n_bonds:
                n_bonds[atom2.get_id()] = set()
            n_bonds[atom2.get_id()].add(atom1.get_id())

        if (
            atom2.get_atomic_number() == 46
            and atom1.get_atomic_number() == 7
        ):
            if atom1.get_id() in n_bonds:
                raise ValueError(f"only one Pd bond per N atom! {bond}")
            n_bonds[atom1.get_id()] = atom2.get_id()

    helicities = []
    for lig_atoms in get_organic_linker_atoms(
        cage=molecule,
        metal_atom_nos=(46,),
    ):
        pairs = []
        for la in lig_atoms:
            if la.get_id() in n_bonds:
                pda = n_bonds[la.get_id()]
                pairs.append((la.get_id(), pda))

        dihedral = abs(
            get_dihedral(
                pt1=pos_mat[pairs[0][0]],
                pt2=pos_mat[pairs[0][1]],
                pt3=pos_mat[pairs[1][1]],
                pt4=pos_mat[pairs[1][0]],
            )
        )

        helicities.append(dihedral)

    return helicities


def get_pore_angle(molecule, metal_atom_num):

    atom_ids = [
        i.get_id()
        for i in molecule.get_atoms()
        if i.get_atomic_number() == metal_atom_num
    ]
    if len(atom_ids) != 2:
        raise ValueError(
            f"{len(atom_ids)} metal atoms found. Expecting 2"
        )

    centroid = molecule.get_centroid()
    v1, v2 = (
        i - centroid for i in molecule.get_atomic_positions(atom_ids)
    )
    aniso_angle = np.degrees(angle_between(v1, v2))

    return aniso_angle


def get_mm_distance(molecule, metal_atom_num):

    atom_ids = [
        i.get_id()
        for i in molecule.get_atoms()
        if i.get_atomic_number() == metal_atom_num
    ]
    if len(atom_ids) != 2:
        raise ValueError(
            f"{len(atom_ids)} metal atoms found. Expecting 2"
        )

    position_matrix = molecule.get_position_matrix()

    distance = euclidean(
        u=position_matrix[atom_ids[0]], v=position_matrix[atom_ids[1]]
    )

    return float(distance)
