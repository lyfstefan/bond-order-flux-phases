"""Independent phase-letter anchors read from the author's original region figures."""

from pathlib import Path
import unittest
import numpy as np
from bond_order.phase_regions import load_regions

ROOT = Path(__file__).resolve().parents[1] / "docs/phase_diagrams/historical_regions"


class HistoricalPhaseRegions(unittest.TestCase):
    def test_original_figure_phase_labels(self):
        anchors = {
            "unconstrained": [
                (0.1, 1.0, "N"),
                (0.21, 4.7, "A"),
                (0.37, 4.7, "B"),
                (0.48, 4.7, "F"),
                (0.69, 4.7, "G"),
                (0.84, 4.7, "D"),
                (0.5, 1.5, "E"),
            ],
            "zero_flux": [
                (0.1, 2.3, "N"),
                (0.21, 4.7, "A"),
                (0.37, 4.7, "B"),
                (0.49, 4.7, "H"),
                (0.64, 4.7, "C"),
                (0.84, 4.7, "D"),
            ],
        }
        for name, points in anchors.items():
            data = load_regions(ROOT / f"{name}.json")
            for filling, interaction, expected in points:
                i = np.argmin(abs(np.asarray(data["interactions"]) - interaction))
                j = np.argmin(abs(np.asarray(data["fillings"]) - filling))
                self.assertEqual(
                    data["phase_labels"][i][j], expected, (name, filling, interaction)
                )


if __name__ == "__main__":
    unittest.main()
