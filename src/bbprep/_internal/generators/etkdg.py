import stk
from rdkit.Chem import AllChem

from bbprep._internal.ensemble.ensemble import Conformer, Ensemble

from .generator import Generator


class ETKDG(Generator):
    """Generate conformers as stk molecules with :func:`AllChem.ETKDGv3()`."""

    def __init__(self, num_confs: int) -> None:
        """Initialise ETKDG generator.

        `v3`:
            New version from DOI: 10.1021/acs.jcim.0c00025
            with improved handling of macrocycles.

        """
        self._num_confs = num_confs

    def generate_conformers(
        self,
        molecule: stk.BuildingBlock,
    ) -> Ensemble:
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit_molecule.RemoveAllConformers()

        params = AllChem.ETKDGv3()  # type: ignore[attr-defined]
        params.randomSeed = 1000
        cids = AllChem.EmbedMultipleConfs(  # type: ignore[attr-defined]
            mol=rdkit_molecule,
            numConfs=self._num_confs,
            params=params,
        )

        ensemble = Ensemble(base_molecule=molecule)
        for cid in cids:
            pos_mat = rdkit_molecule.GetConformer(id=cid).GetPositions()
            ensemble.add_conformer(
                conformer=Conformer(
                    molecule=molecule.with_position_matrix(pos_mat),
                    conformer_id=cid,
                    source="etkdg-v3",
                    permutation=None,
                ),
            )

        return ensemble
