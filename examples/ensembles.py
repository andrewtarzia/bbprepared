import logging
from pathlib import Path

import stk
import stko

import bbprep

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> None:
    """Run the example."""
    bb1 = stk.BuildingBlock("NCCNCCN", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            orientations=(0, 0),
            num_repeating_units=1,
        )
    )

    examples_output = Path("out_ensembles")
    examples_output.mkdir(parents=True, exist_ok=True)

    calculator = bbprep.EnergyCalculator(
        name="MMFFEnergy",
        function=stko.MMFFEnergy().get_energy,
    )

    optimiser = bbprep.Optimiser(
        name="MMFF",
        function=stko.MMFF().optimize,
    )

    generator = bbprep.generators.ETKDG(num_confs=30)
    ensemble = generator.generate_conformers(
        stk.BuildingBlock.init_from_molecule(polymer)
    )
    logger.info(ensemble)

    # Get lowest energy without opt.
    lowest_energy_conformer = ensemble.get_lowest_energy_conformer(
        calculator=calculator,
    )
    logger.info(calculator.function(lowest_energy_conformer.molecule))
    lowest_energy_conformer.molecule.write(examples_output / "low_e_1.mol")

    # Optimise ensemble.
    opt_ensemble = ensemble.optimise_conformers(
        optimiser=optimiser,
    )

    # Get lowest energy conformer.
    lowest_energy_conformer = opt_ensemble.get_lowest_energy_conformer(
        calculator=calculator
    )
    logger.info(calculator.function(lowest_energy_conformer.molecule))
    lowest_energy_conformer.molecule.write(examples_output / "low_e_2.mol")

    logger.info(lowest_energy_conformer)


if __name__ == "__main__":
    main()
