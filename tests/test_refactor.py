import unittest
from unittest.mock import patch
import numpy as np
from bond_order import FluxModel
from bond_order.occupation import occupy
from bond_order.optimization import rlex_values,lex_sectors
from reference import parallel_integral as original
from reference import zero_flux


class Refactor(unittest.TestCase):
    def test_recovered_bands_and_energy(self):
        model=FluxModel(mesh=16)
        points=[[-.32,.11,-.07,.19],[-.21,-.17,.09,.13],[-.4,.4,-.4,.4]]
        with patch.multiple(original,Nx=16,Ny=16,siteN=256,t=1.,tp=.2,error=3):
            for j in points:
                expected=np.stack([original.E1(model.kx,model.ky,*j,4),original.E2(model.kx,model.ky,*j,4)])
                np.testing.assert_allclose(model.bands(j),expected,atol=2e-13)
                for filling in [.23,.4,.5,.7]:
                    self.assertAlmostEqual(model.energy(j,filling),original.Etot(j,256*filling,4),places=12)

    def test_hermitian_matrix_dispersion(self):
        model=FluxModel(mesh=20)
        j=np.array([-.32,.11,-.07,.19])
        expanded=model.bands(j)
        np.testing.assert_allclose(expanded,model.bands(j,"complex"),atol=5e-13)
        for row,col in [(0,0),(2,8),(13,4)]:
            kx,ky=model.kx[row,col],model.ky[row,col]
            h=-sum((1+4j*j[q])*np.exp(1j*(kx*dx+ky*dy)) for q,(dx,dy) in enumerate([(1,0),(0,1),(-1,0),(0,-1)]))
            d=-.8*np.cos(kx)*np.cos(ky)
            np.testing.assert_allclose(expanded[:,row,col],np.linalg.eigvalsh([[d,h],[h.conjugate(),d]]),atol=1e-13)

    def test_occupation_modes(self):
        # Degenerate shell: the original tolerance cannot guarantee exact N.
        bands=np.array([[-2,-1,-1,-1],[1,1,1,2]])
        legacy=occupy(bands,.375,"legacy",0)
        exact=occupy(bands,.375,"canonical")
        self.assertAlmostEqual(exact.actual_filling,.375)
        self.assertEqual(exact.particle_error,0)
        self.assertFalse(legacy.bisection_converged)
        self.assertNotEqual(legacy.actual_filling,.375)

    def test_symmetry_and_canonical_counting(self):
        model=FluxModel(mesh=24,tp=0,occupation="canonical")
        j=np.array([-.32,.11,-.07,.19])
        energy=model.energy(j,.37)
        for other in [-j,np.roll(j,1)]:
            self.assertAlmostEqual(model.energy(other,.37),energy,places=12)
        self.assertAlmostEqual(model.energy(j,.63),energy,places=12)
        for filling in [0,.001,.371,.5,1]:
            result=model.evaluate(j,filling)
            self.assertEqual(result["actual_filling"],filling)

    def test_rlex_matches_paper_region(self):
        for j in np.random.default_rng(2026).uniform(-1,1,(1000,4)):
            source=(j[0]<=0 and np.all(rlex_values(j)>=0))
            paper=(j[0]<=j[1]<=j[3] and j[0]<=j[2] and j[0]<=0
                   and all(j[q]<=-j[0] for q in [1,2,3]))
            self.assertEqual(source,paper)

    def test_lex_constraints_match_recovered_source_including_slack(self):
        indices=[[0,1,2,9,10,11],[3,4,5,9,10,11],[6,7,9,10,11]]
        for j in np.random.default_rng(2026).uniform(-1,1,(20,4)):
            for index,constraints in enumerate(lex_sectors()):
                expected=[getattr(zero_flux,f"constraint{i}")(j) for i in indices[index]]
                np.testing.assert_allclose(constraints[0]["fun"](j),expected,atol=1e-15)
                if index==2:
                    self.assertEqual(constraints[1]["fun"](j),zero_flux.constraint8(j))

    def test_domain(self):
        with self.assertRaises(ValueError):FluxModel(mesh=15)
        with self.assertRaises(ValueError):FluxModel().energy([0]*4,1.1)


if __name__=="__main__":unittest.main()
