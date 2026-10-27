"""SHGO and source symmetry/zero-flux constraints, without top-level sweeps."""
import numpy as np
from scipy.optimize import shgo


def rlex_values(j):
    """Active four inequalities in parallel calculation[_integral].py.

    Together with J1 in [-1,0], these imply the paper's relaxed inequalities
    J1<=J2<=J4, J1<=J3 and J2,J3,J4<=-J1.
    """
    return np.array([-j[0]-abs(j[1]),-j[0]-abs(j[2]),-j[0]-j[3],j[3]-j[1]])


def lex_sectors():
    """Original zero-flux sectors, including their historical 1e-5 slack.

    The equality J2-J1+1e-5=0 is also retained, rather than silently
    changing it to exact equality of the two currents.
    """
    slack=1e-5
    common=lambda j:np.array([-j[0]-j[1],-j[0]-j[2],-j[0]-j[3]])+slack
    return [
        [{"type":"ineq","fun":lambda j:np.r_[np.array([j[3]-j[1],j[1]-j[2],j[2]-j[0]])+slack,common(j)]}],
        [{"type":"ineq","fun":lambda j:np.r_[np.array([j[2]-j[1],j[1]-j[0],j[3]-j[1]])+slack,common(j)]}],
        [{"type":"ineq","fun":lambda j:np.r_[np.array([j[3]-j[2],j[2]-j[1]])+slack,common(j)]},
         {"type":"eq","fun":lambda j:j[1]-j[0]+slack}]]


def solve(model,filling,zero_flux=False,samples=256,iterations=2,symmetry=None):
    """Bounded source-method optimization with every solver status reported.

    Ordinary mode defaults to source RLEX. Zero-flux mode defaults to the
    source's three LEX sectors; RLEX may be selected explicitly. This does
    not certify global optimality or convergence of phase boundaries.
    """
    if not np.isfinite(filling) or not 0<=filling<=1:raise ValueError("filling must be in [0,1]")
    if not isinstance(samples,int) or samples<16 or not isinstance(iterations,int) or iterations<1:
        raise ValueError("samples >= 16 and iterations >= 1 required")
    symmetry=symmetry or ("legacy_lex" if zero_flux else "rlex")
    if symmetry=="rlex":sectors=[[{"type":"ineq","fun":rlex_values}]]
    elif symmetry=="legacy_lex":sectors=lex_sectors()
    else:raise ValueError("symmetry must be rlex or legacy_lex")
    if zero_flux and model.t<=0:raise ValueError("zero-flux convention requires t > 0")
    bounds=[(-1,0),(-1,1),(-1,1),(-1,1)]
    candidates=[]
    best_j=np.zeros(4)
    best_energy=model.energy(best_j,filling)
    improved=False
    for index,base in enumerate(sectors):
        constraints=base+([{"type":"eq","fun":model.square_flux}] if zero_flux else [])
        result=shgo(lambda j:model.energy(j,filling),bounds,sampling_method="sobol",n=samples,
                    iters=iterations,constraints=constraints)
        j=getattr(result,"x",None)
        feasible=False
        violation=None
        if j is not None and np.isfinite(j).all():
            violations=[max(0.,float(np.max(np.maximum(-np.asarray(c["fun"](j)),0)))) if c["type"]=="ineq"
                        else float(np.max(abs(np.asarray(c["fun"](j))))) for c in constraints]
            violations += [float(max(0,-1-j[0],j[0],np.max(np.abs(j[1:]))-1))]
            violation=max(violations)
            feasible=violation<=1e-6
            energy=model.energy(j,filling)
            if feasible and energy<best_energy:
                best_j,best_energy=j.copy(),energy;improved=True
        candidates.append({"sector":index,"success":bool(result.success),"message":str(result.message),
                           "feasible":feasible,"max_constraint_violation":violation,
                           "currents":None if j is None else j.tolist(),
                           "energy":None if j is None else float(model.energy(j,filling))})
    return {"parameters":model.parameters(),"currents":best_j.tolist(),**model.evaluate(best_j,filling),
            "settings":{"method":"SHGO","samples":samples,"iterations":iterations,"sampling":"sobol",
                        "symmetry":symmetry,"legacy_lex_slack":1e-5 if symmetry=="legacy_lex" else 0.,
                        "zero_flux":zero_flux,"bounds":bounds},
            "square_flux":model.square_flux(best_j) if model.t>0 else None,
            "selected":"solver candidate" if improved else "zero-current baseline",
            "solver_candidates":candidates,"any_feasible_solver_candidate":any(c["feasible"] for c in candidates)}
