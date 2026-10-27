"""Export existing public scan values without rerunning searches or refinements."""

import argparse, json
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from bond_order import FluxModel
from bond_order.plotting import plot_phase_diagram, plot_phase_comparison, plot_bands
from bond_order.publication_style import apply_style, save


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("docs"))
    parser.add_argument("--font-path", type=Path)
    args = parser.parse_args()
    root = args.root

    def export(fig, path):
        apply_style(font_path=args.font_path)
        save(fig, path)
        fig.savefig(path.parent / (path.name + ".svg"))
        plt.close(fig)

    scans = root / "phase_diagrams"
    for name in ["unconstrained", "zero_flux"]:
        rows = json.loads((scans / name / "sweep.json").read_text())
        export(plot_phase_diagram(rows), scans / name / "phase_diagram")
        rows = json.loads((scans / "audited" / f"{name}.json").read_text())
        export(plot_phase_diagram(rows), scans / "audited" / f"{name}_order_parameters")
    a, b = [
        json.loads((scans / "audited" / f"{name}.json").read_text())
        for name in ["unconstrained", "zero_flux"]
    ]
    export(plot_phase_comparison(a, b), scans / "phase_diagram_comparison")
    settings = json.loads((root / "bands_parameters.json").read_text())
    export(
        plot_bands(FluxModel(**settings["parameters"]), settings["currents"]),
        root / "bands",
    )


if __name__ == "__main__":
    main()
