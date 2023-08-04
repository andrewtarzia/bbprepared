import stk

from .selector import Selector


class DeletersSelector(Selector):
    """
    Select atom ids in stk molecules by deleters.

    """

    def __init__(self):
        """
        Initialise Selector.

        """
        pass

    def select_atoms(self, molecule: stk.Molecule) -> tuple[int]:
        assert molecule.get_num_functional_groups() > 0

        atoms = []
        for fg in molecule.get_functional_groups():
            for id_ in fg.get_deleter_ids():
                atoms.append(id_)

        return tuple(atoms)
