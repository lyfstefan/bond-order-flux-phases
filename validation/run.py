from pathlib import Path
from unittest.mock import patch
import json,platform
import numpy as np
import scipy
from bond_order import FluxModel,solve
from bond_order.optimization import rlex_values
from reference import parallel_integral,zero_flux


def main():
    report={"environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
            "source_status":"Confirmed project directory; final historical figure-generation version unknown",
            "source_equivalence":[],"optimization_examples":[]}
    model=FluxModel(mesh=20)
    for name,ref in [("parallel calculation_integral.py",parallel_integral),
                     ("phase diagram J1 J2 J3 J4_zero flux subspace.py",zero_flux)]:
        max_band,max_energy=0.,0.
        with patch.multiple(ref,Nx=20,Ny=20,siteN=400,t=1.,tp=.2,error=3):
            for j in [[-.32,.11,-.07,.19],[-.21,-.17,.09,.13],[-.4,.4,-.4,.4]]:
                bands=np.stack([ref.E1(model.kx,model.ky,*j,4),ref.E2(model.kx,model.ky,*j,4)])
                max_band=max(max_band,float(np.max(abs(model.bands(j)-bands))))
                for filling in [.23,.4,.5,.7]:
                    max_energy=max(max_energy,abs(model.energy(j,filling)-ref.Etot(j,400*filling,4)))
        assert max_band<1e-11 and max_energy<1e-11
        report["source_equivalence"].append({"source":name,"mesh":20,"band_max_absolute_error":max_band,
            "energy_max_absolute_error":max_energy,"currents_tested":3,"fillings_tested":[.23,.4,.5,.7]})
    trial=[-.32,.11,-.07,.19]
    report["occupation_mode_comparison"]=[{"mode":mode,**FluxModel(mesh=24,occupation=mode).evaluate(trial,.371)} for mode in ["legacy","canonical"]]
    report["canonical_fixed_trial_mesh_sequence"]=[{"mesh":n,"energy":FluxModel(mesh=n,occupation="canonical").energy(trial,.37)} for n in [24,40,80,120]]
    # Independent continuum normalization check, not an original paper number.
    free=FluxModel(mesh=120,t=.3,tp=0,interaction=0,occupation="canonical").energy([0]*4,.5)
    exact=-8*.3/np.pi**2
    assert abs(free-exact)<1e-4
    report["free_lattice_half_filling"]={"mesh":120,"computed":free,"continuum":exact,"absolute_error":abs(free-exact)}
    # Exercise both retained occupation modes and both physical subspaces.
    # These small meshes/search budgets are workflow checks, not full figures.
    for mode,zero in [("legacy",False),("canonical",False),("canonical",True)]:
        model=FluxModel(mesh=12,occupation=mode)
        result=solve(model,.4,zero_flux=zero,samples=32,iterations=1)
        assert np.isfinite(result["energy_per_site"])
        assert result["energy_per_site"]<=model.energy([0]*4,.4)+1e-12
        assert result["any_feasible_solver_candidate"]
        if zero:assert abs(result["square_flux"])<1e-6
        report["optimization_examples"].append(result)
    # Re-evaluate versioned scan data without redoing 760 global searches.
    phase_root=Path(__file__).resolve().parents[1]/"docs/phase_diagrams"
    if (phase_root/"audit.json").exists():
        loaded={name:json.loads((phase_root/"audited"/f"{name}.json").read_text()) for name in ["unconstrained","zero_flux"]}
        max_error=0.;max_zero_flux=0.;maps={}
        for name,rows in loaded.items():
            maps[name]={}
            assert len(rows)==380
            for row in rows:
                model=FluxModel(**row["parameters"])
                value=model.energy(row["currents"],row["requested_filling"])
                max_error=max(max_error,abs(value-row["energy_per_site"]))
                assert min(rlex_values(row["currents"]))>=-1e-6
                assert -1-1e-6<=row["currents"][0]<=1e-6 and max(abs(np.asarray(row["currents"])))<=1+1e-6
                if name=="zero_flux":max_zero_flux=max(max_zero_flux,abs(model.square_flux(row["currents"])))
                maps[name][(model.interaction,row["requested_filling"])]=value
        assert maps["unconstrained"].keys()==maps["zero_flux"].keys()
        gap=max(maps["unconstrained"][k]-maps["zero_flux"][k] for k in maps["unconstrained"])
        assert max_error<1e-9 and max_zero_flux<1e-6 and gap<1e-8
        report["versioned_phase_map_checks"]={"points":760,"global_searches_repeated":False,
            "max_stored_vs_reevaluated_energy_error":max_error,"max_zero_flux_residual":max_zero_flux,
            "max_variational_nesting_violation":gap}
    report["scope"]="Import-safe source refactor; source-band/energy comparisons, normalization and bounded optimization checks. No final-version identification or full paper phase-diagram reproduction. Loop-current ansatz only; no CDW competition."
    report["all_checks_passed"]=True
    Path(__file__).with_name("report.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"all_checks_passed":True,"source_equivalence":report["source_equivalence"],
                      "optimization_checks":len(report["optimization_examples"])},indent=2))


if __name__=="__main__":main()
