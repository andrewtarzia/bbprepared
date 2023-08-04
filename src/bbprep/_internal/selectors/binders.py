import stk

from .selector import Selector


class BindersSelector(Selector):
    """
    Select atom ids in stk molecules by binders.

    """

    def select_atoms(self, molecule: stk.Molecule) -> tuple[int]:
        assert molecule.get_num_functional_groups() > 0

        atoms = []
        for fg in molecule.get_functional_groups():
            for id_ in fg.get_bonder_ids():
                atoms.append(id_)

        return tuple(atoms)
