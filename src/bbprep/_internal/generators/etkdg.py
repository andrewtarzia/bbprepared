import typing

import stk
from rdkit.Chem import AllChem as rdkit


class ETKDG:
    """
    Generate conformers as stk molecules with ETKDG__.

    __ rdkit

    """

    def __init__(self, num_confs):
        """
        Initialise ETKDG generator.

        `v3`:
            New version from DOI: 10.1021/acs.jcim.0c00025
            with improved handling of macrocycles.

        """
        self._num_confs = num_confs

    def generate_conformers(
        self,
        molecule: stk.Molecule,
    ) -> typing.Iterable[stk.Molecule]:
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit_molecule.RemoveAllConformers()

        params = rdkit.ETKDGv3()
        params.randomSeed = 1000
        cids = rdkit.EmbedMultipleConfs(
            mol=rdkit_molecule,
            numConfs=self._num_confs,
            params=params,
        )

        for cid in cids:
            pos_mat = rdkit_molecule.GetConformer(id=cid).GetPositions()
            yield molecule.with_position_matrix(pos_mat)
