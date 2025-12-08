import logging
import os
from pathlib import Path

import ensembles

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> None:
    """Run the example."""
    init_dir = Path.cwd()
    os.chdir("examples/")

    try:
        ensembles.main()
        logger.info("all examples ran, at least!")

    finally:
        os.chdir(init_dir)


if __name__ == "__main__":
    main()
