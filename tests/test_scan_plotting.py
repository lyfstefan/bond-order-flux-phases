import unittest
import copy
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
from bond_order.plotting import scan_grid, plot_phase_diagram, flux_mod_2pi, plot_bands
from bond_order import FluxModel
from bond_order.optimization import refine_from_seed, rlex_values


def fixture():
    return [
        {
            "requested_filling": nu,
            "parameters": {
                "interaction": v,
                "mesh": 12,
                "t": 1.0,
                "tp": 0.2,
                "occupation": "canonical",
                "particle_tolerance": 3.0,
            },
            "settings": {
                "samples": 32,
                "iterations": 1,
                "symmetry": "rlex",
                "zero_flux": False,
            },
            "currents": [-nu, 0.1 * v, 0.02, 0.04],
            "energy_per_site": -v - nu,
            "solver_candidates": [{"success": True, "feasible": True}],
        }
        for v in [1.0, 2.0]
        for nu in [0.1, 0.3, 0.8]
    ]


class ScanPlotting(unittest.TestCase):
    def test_shuffled_nonuniform_coordinates_and_observables(self):
        results = fixture()
        grid = scan_grid(results[::-1])
        np.testing.assert_allclose(grid["fillings"], [0.1, 0.3, 0.8])
        np.testing.assert_allclose(grid["currents"][1, 2], [-0.8, 0.2, 0.02, 0.04])
        phi = np.arctan(2 * np.array([-0.8, 0.2, 0.02, 0.04]))
        self.assertAlmostEqual(
            grid["square_flux_over_pi"][1, 2],
            (phi[3] - phi[0] + phi[1] - phi[2]) / np.pi,
        )
        fig = plot_phase_diagram(results)
        fig.canvas.draw()
        plt.close(fig)

    def test_lower_and_upper_band_labels(self):
        fig = plot_bands(FluxModel(mesh=12), np.zeros(4))
        labels = [
            child._publication_label
            for ax in fig.axes
            for child in ax.child_axes
            if hasattr(child, "_publication_label")
        ]
        self.assertEqual(labels, [r"$E_-$", r"$E_+$"])
        plt.close(fig)

    def test_oriented_flux_is_periodic_and_zero_roundoff_is_removed(self):
        raw = np.array([-np.pi / 2, 0.0, np.pi, 2 * np.pi, 9 * np.pi / 2, -1e-12])
        wrapped = flux_mod_2pi(raw)
        np.testing.assert_allclose(
            wrapped, [3 * np.pi / 2, 0, np.pi, 0, np.pi / 2, 0], atol=1e-10
        )
        np.testing.assert_allclose(np.exp(1j * wrapped), np.exp(1j * raw), atol=1e-10)
        self.assertTrue(np.all((wrapped >= 0) & (wrapped < 2 * np.pi)))

    def test_reject_incomplete_duplicate_or_mixed_data(self):
        results = fixture()
        for bad in [results[:-1], results + [results[0]]]:
            with self.assertRaises(ValueError):
                scan_grid(bad)
        bad = copy.deepcopy(results)
        bad[0]["parameters"]["mesh"] = 16
        with self.assertRaises(ValueError):
            scan_grid(bad)

    def test_admissible_refinement_retains_seed_bound_and_constraints(self):
        model = FluxModel(mesh=12, occupation="canonical")
        seed = [-0.3, -0.3, 0.3, 0.3]
        result = refine_from_seed(model, 0.4, seed, zero_flux=True)
        self.assertLessEqual(result["energy_per_site"], model.energy(seed, 0.4) + 1e-12)
        self.assertLess(abs(result["square_flux"]), 1e-6)
        self.assertGreaterEqual(float(np.min(rlex_values(result["currents"]))), -1e-6)
        with self.assertRaises(ValueError):
            refine_from_seed(model, 0.4, [0.3, 0, 0, 0])


if __name__ == "__main__":
    unittest.main()
