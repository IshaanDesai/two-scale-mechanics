import argparse
from .config import Config
from .simulation import Simulation

def main():
    parser = argparse.ArgumentParser(
        description='Runs configurable solid mechanics macro sim using fenicsx',
        usage='%(prog)s [options]'
    )
    parser.add_argument(
        'path',
        type=str,
        help='Path to the config file'
    )
    args = parser.parse_args()

    config = Config()
    config.load(args.path)

    sim = Simulation.generate(config)
    sim.run()