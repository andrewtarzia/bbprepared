Torsion scan.
=============

This is a simple script for selecting a torsion, scanning around it and saving
the energy and the torsion value. Changes will be necessary for more than one
torsion.


.. testcode:: recipe5-test

    import stk
    import bbprepared

    smiles = "C(C#CC1C=CC=C(Br)C=1)1C=CC=C(Br)C=1"

    building_block = stk.BuildingBlock(
        smiles=smiles,
        functional_groups=(stk.BromoFactory(),),
    )

    # Rotate around alkynes.
    smarts_to_rotate = "[#6X3][#6X3][#6X2H0]#!@[#6X2H0][#6X3][#6X3]"
    selected_indices = (0, 1, 4, 5)
    expected_num_atoms = 6

    selector = bbprep.selectors.BySmartsSelector(
        smarts=smarts_to_rotate,
        selected_indices=selected_indices,
    )

    # Scan with MMFF.
    generator = TorsionScanner(
        target_torsions=(
            bbprep.generators.TorsionRange(
                smarts=smarts_to_rotate,
                expected_num_atoms=expected_num_atoms,
                scanned_ids=selected_indices,
                scanned_range=range(0, 359, 20),
            ),
        ),
    )
    ensemble = generator.generate_conformers(molecule)

    energies = []
    torsions = []
    for conformer in ensemble.yield_conformers():
        molecule = conformer.molecule.with_centroid(np.array((0, 0, 0)))

        # Align to conformer 0.
        if conformer.conformer_id == 0:
            first_pos_mat = molecule.get_position_matrix()
            aligned_ = molecule.clone()

        # Else, align.
        else:
            current_pos_mat = molecule.get_position_matrix()
            alignment = rmsd.kabsch(current_pos_mat, first_pos_mat)
            aligned_pos_mat = np.dot(current_pos_mat, alignment)
            aligned_ = molecule.with_position_matrix(aligned_pos_mat)

        # Write to file.
        # aligned_.write(
        #     f"{name}_d{conformer.conformer_id}.mol"
        # )

        energy = (
            stko.MMFFEnergy(ignore_inter_interactions=False).get_energy(
                aligned_
            )
            * 4.184
        )

        atom_positions = aligned_.get_position_matrix()
        for found in selector.yield_stepwise(aligned_):
            found_str = "_".join(str(i) for i in found)
            if not any(i in atoms_to_be_constrained for i in found):
                for i in found:
                    atoms_to_be_constrained.add(i)

                torsion = float(
                    stko.calculate_dihedral(
                        pt1=atom_positions[found[0]],
                        pt2=atom_positions[found[1]],
                        pt3=atom_positions[found[2]],
                        pt4=atom_positions[found[3]],
                    )
                )

                torsions.append(torsion)

    # fig, ax = plt.subplots(figsize=(8, 5))
    # ax.scatter(
    #     torsions,
    #     energies,
    #     edgecolor="k",
    #     c="tab:blue,
    #     s=100,
    # )
    # ax.tick_params(axis="both", which="major", labelsize=16)
    # ax.set_xlabel("measured torsion [deg]", fontsize=16)
    # ax.set_ylabel("rel. MMFF energy [kJ.mol-1]", fontsize=16)
    # ax.set_xlim(None, None)
    # ax.set_ylim(0, max(energies))
    # ax.legend(fontsize=16)
    # fig.tight_layout()
    # fig.savefig(
    #     figure_dir / f"scan_{name}.png",
    #     dpi=360,
    #     bbox_inches="tight",
    # )
    # plt.close()

.. testcode:: recipe5-test
    :hide:

    import numpy as np

    assert np.isclose(
        process.calculate_score(
            conformer=min_molecule,
            conformer_id=min_molecule.conformer_id,
        ),
        0.8249489663822132,
    )
