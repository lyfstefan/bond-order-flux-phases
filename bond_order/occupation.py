"""Discrete zero-temperature filling, with explicit legacy/canonical modes."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Occupation:
    chemical_potential: float
    electronic_energy_per_site: float
    actual_filling: float
    particle_error: float
    bisection_converged: bool


def occupy(bands, filling, mode="legacy", particle_tolerance=3.0):
    """Every full-BZ folded band entry has weight 1/2 per physical state.

    Legacy reproduces source Mu12's fixed [-1000,1000] bracket, <= occupation,
    half-count and tolerance/2 stopping condition, including its 99 iterations.
    Canonical instead sorts and fractionally occupies the last discrete entry.
    Neither mode introduces a finite temperature or interpolates the spectrum.
    """
    levels=np.asarray(bands,dtype=float).ravel()
    if not np.isfinite(levels).all() or not np.isfinite(filling) or not 0<=filling<=1:
        raise ValueError("finite levels and filling in [0,1] required")
    sites=len(levels)/2
    if mode=="canonical":
        ordered=np.sort(levels)
        number=len(levels)*filling
        whole=int(np.floor(number))
        fraction=number-whole
        energy=float(ordered[:whole].sum())
        if whole<len(ordered): energy+=fraction*ordered[whole]
        if whole==0:mu=ordered[0]
        elif whole==len(ordered):mu=ordered[-1]
        elif fraction:mu=ordered[whole]
        else:mu=(ordered[whole-1]+ordered[whole])/2
        return Occupation(float(mu),energy/len(levels),float(filling),0.0,True)
    if mode!="legacy":
        raise ValueError("occupation must be legacy or canonical")
    if not np.isfinite(particle_tolerance) or particle_tolerance<0:
        raise ValueError("particle tolerance must be finite and nonnegative")
    if levels.min() < -1000 or levels.max() > 1000:
        raise ValueError("source legacy mu bracket is insufficient; use canonical mode")
    lo,hi=-1000.,1000.
    converged=False
    target=sites*filling
    for _ in range(99):
        mu=(lo+hi)/2
        count=np.count_nonzero(levels<=mu)/2
        if count>target+particle_tolerance/2:hi=mu
        elif count<target-particle_tolerance/2:lo=mu
        else:
            converged=True
            break
    mask=levels<=mu
    count=np.count_nonzero(mask)/2
    return Occupation(float(mu),float(levels[mask].sum()/len(levels)),
                      float(count/sites),float(count-target),converged)
