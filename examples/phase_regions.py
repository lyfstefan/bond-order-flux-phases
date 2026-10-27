"""Replot the recovered author phase-region maps, with no new classification."""

import argparse
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
from bond_order.phase_regions import export_regions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("docs/phase_diagrams"))
    args = parser.parse_args()
    export_regions(args.root)


if __name__ == "__main__":
    main()
