import unittest
import copy
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
from bond_order.plotting import scan_grid,plot_phase_diagram
from bond_order import FluxModel
from bond_order.optimization import refine_from_seed,rlex_values


def fixture():
    return [{"requested_filling":nu,"parameters":{"interaction":v,"mesh":12,"t":1.,"tp":.2,
             "occupation":"canonical","particle_tolerance":3.},
             "settings":{"samples":32,"iterations":1,"symmetry":"rlex","zero_flux":False},
             "currents":[-nu,.1*v,.02,.04],"energy_per_site":-v-nu,
             "solver_candidates":[{"success":True,"feasible":True}]} for v in [1.,2.] for nu in [.1,.3,.8]]


class ScanPlotting(unittest.TestCase):
    def test_shuffled_nonuniform_coordinates_and_observables(self):
        results=fixture();grid=scan_grid(results[::-1])
        np.testing.assert_allclose(grid["fillings"],[.1,.3,.8])
        np.testing.assert_allclose(grid["currents"][1,2],[-.8,.2,.02,.04])
        phi=np.arctan(2*np.array([-.8,.2,.02,.04]))
        self.assertAlmostEqual(grid["square_flux_over_pi"][1,2],(phi[3]-phi[0]+phi[1]-phi[2])/np.pi)
        fig=plot_phase_diagram(results);fig.canvas.draw();plt.close(fig)

    def test_reject_incomplete_duplicate_or_mixed_data(self):
        results=fixture()
        for bad in [results[:-1],results+[results[0]]]:
            with self.assertRaises(ValueError):scan_grid(bad)
        bad=copy.deepcopy(results);bad[0]["parameters"]["mesh"]=16
        with self.assertRaises(ValueError):scan_grid(bad)

    def test_admissible_refinement_retains_seed_bound_and_constraints(self):
        model=FluxModel(mesh=12,occupation="canonical")
        seed=[-.3,-.3,.3,.3]
        result=refine_from_seed(model,.4,seed,zero_flux=True)
        self.assertLessEqual(result["energy_per_site"],model.energy(seed,.4)+1e-12)
        self.assertLess(abs(result["square_flux"]),1e-6)
        self.assertGreaterEqual(float(np.min(rlex_values(result["currents"]))),-1e-6)
        with self.assertRaises(ValueError):refine_from_seed(model,.4,[.3,0,0,0])


if __name__=="__main__":unittest.main()
