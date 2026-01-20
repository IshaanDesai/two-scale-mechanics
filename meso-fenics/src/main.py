import argparse
from .config import Config
from .simulation import Simulation


def main():
    """
    Main entry point for the FEniCSx macro-scale simulation.

    Parses command-line arguments, loads configuration, and runs
    the appropriate simulation type.

    Command-line Arguments
    ----------------------
    path : str
        Path to the configuration JSON file.
    """
    parser = argparse.ArgumentParser(
        description="Runs configurable solid mechanics macro sim using fenicsx",
        usage="%(prog)s [options]",
    )
    parser.add_argument("path", type=str, help="Path to the config file")
    args = parser.parse_args()

    config = Config()
    config.load(args.path)

    sim = Simulation.generate(config)
    sim.run()
