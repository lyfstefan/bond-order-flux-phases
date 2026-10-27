"""Four real bond currents, directions (+x,+y,-x,-y), energies per site.

The expanded dispersion is refactored from the identical E1/E2 expressions
in the four recovered author scripts. No script or global plotting state is
imported. The compact complex-hopping expression is available independently.
"""
from dataclasses import dataclass,asdict
import numpy as np
from .occupation import occupy


@dataclass(frozen=True)
class FluxModel:
    mesh: int=100
    t: float=1.0
    tp: float=.2
    interaction: float=4.0
    occupation: str="legacy"
    particle_tolerance: float=3.0

    def __post_init__(self):
        if not isinstance(self.mesh,int) or self.mesh<4 or self.mesh%2:
            raise ValueError("mesh must be an even integer >= 4")
        if not np.isfinite([self.t,self.tp,self.interaction,self.particle_tolerance]).all():
            raise ValueError("model parameters must be finite")
        if self.t<0 or self.interaction<0 or self.particle_tolerance<0:
            raise ValueError("t, interaction and particle tolerance must be nonnegative")
        if self.occupation not in ("legacy","canonical"):
            raise ValueError("occupation must be legacy or canonical")
        # Use the same coordinates as the original research programs.
        k=np.arange(0,2*np.pi,2*np.pi/self.mesh)-np.pi
        if len(k)!=self.mesh:raise ValueError("unexpected momentum-grid length")
        kx,ky=np.meshgrid(k,k)
        object.__setattr__(self,"kx",kx)
        object.__setattr__(self,"ky",ky)

    @staticmethod
    def currents(values):
        j=np.asarray(values,dtype=float)
        if j.shape!=(4,) or not np.isfinite(j).all():
            raise ValueError("provide four finite real currents")
        return j

    def bands(self,currents,representation="expanded"):
        j1,j2,j3,j4=self.currents(currents)
        x,y,t,v=self.kx,self.ky,self.t,self.interaction
        if representation=="expanded":
            # The author's A-F expansion is retained, including sqrt(abs(q)).
            # q is |h_AB|^2, so this is a square root of a quadratic expression.
            a=2*(t*t+j1*j3*v*v)
            b=2*(4*t*t+(j1+j3)*(j2+j4)*v*v)
            c=2*(t*t+j2*j4*v*v)
            d=2*(j3-j1)*t*v
            e=2*(j4-j2)*t*v
            f=2*(j3-j1)*(j4-j2)*v*v
            q=(4*t*t+(j1*j1+j2*j2+j3*j3+j4*j4)*v*v
               +a*np.cos(2*x)+b*np.cos(x)*np.cos(y)+c*np.cos(2*y)
               +2*e*np.cos(x)*np.sin(y)+d*np.sin(2*x)+2*d*np.cos(y)*np.sin(x)
               +f*np.sin(x)*np.sin(y)+e*np.sin(2*y))
            magnitude=np.sqrt(np.abs(q))
        elif representation=="complex":
            phases=np.stack([np.exp(1j*x),np.exp(1j*y),np.exp(-1j*x),np.exp(-1j*y)])
            magnitude=np.abs(np.sum((t+1j*v*np.array([j1,j2,j3,j4]))[:,None,None]*phases,axis=0))
        else:raise ValueError("representation must be expanded or complex")
        diagonal=-4*self.tp*np.cos(x)*np.cos(y)
        return np.stack([diagonal-magnitude,diagonal+magnitude])

    def evaluate(self,currents,filling):
        j=self.currents(currents)
        occ=occupy(self.bands(j),filling,self.occupation,self.particle_tolerance)
        penalty=float(self.interaction/4*(j@j))
        return {**asdict(occ),"requested_filling":float(filling),"current_cost_per_site":penalty,
                "energy_per_site":occ.electronic_energy_per_site+penalty}

    def energy(self,currents,filling):
        return self.evaluate(currents,filling)["energy_per_site"]

    def square_flux(self,currents):
        if self.t<=0:raise ValueError("source zero-flux convention requires t > 0")
        phi=np.arctan(self.interaction*self.currents(currents)/self.t)
        return float(phi[3]-phi[0]+phi[1]-phi[2])

    def parameters(self):
        return {k:getattr(self,k) for k in ["mesh","t","tp","interaction","occupation","particle_tolerance"]}
