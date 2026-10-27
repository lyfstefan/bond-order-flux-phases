"""Explicit one-point and sweep commands; no calculation is triggered on import."""
import argparse,json
from pathlib import Path
from .model import FluxModel
from .optimization import solve
from .output import write_json,save_band_data


def main():
    parser=argparse.ArgumentParser(description="Refactored bond-order research calculations")
    parser.add_argument("command",choices=["one-point","sweep"])
    parser.add_argument("--mesh",type=int,default=24)
    parser.add_argument("--t",type=float,default=1.)
    parser.add_argument("--tp",type=float,default=.2)
    parser.add_argument("--V",type=float,default=4.)
    parser.add_argument("--filling",type=float,default=.4)
    parser.add_argument("--occupation",choices=["legacy","canonical"],default="legacy")
    parser.add_argument("--particle-tolerance",type=float,default=3.)
    parser.add_argument("--zero-flux",action="store_true")
    parser.add_argument("--symmetry",choices=["rlex","legacy_lex"])
    parser.add_argument("--samples",type=int,default=64)
    parser.add_argument("--iterations",type=int,default=1)
    parser.add_argument("--fillings",help="comma-separated fillings for sweep")
    parser.add_argument("--interactions",help="comma-separated V values for sweep")
    parser.add_argument("--output",type=Path,default=Path("examples/output"))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    def point(v,nu):
        model=FluxModel(args.mesh,args.t,args.tp,v,args.occupation,args.particle_tolerance)
        result=solve(model,nu,args.zero_flux,args.samples,args.iterations,args.symmetry)
        return model,result
    if args.command=="sweep":
        if not args.fillings or not args.interactions:parser.error("sweep requires --fillings and --interactions")
        fillings=[float(v) for v in args.fillings.split(",")]
        interactions=[float(v) for v in args.interactions.split(",")]
        results=[]
        for v in interactions:
            for nu in fillings:
                _,result=point(v,nu);results.append(result)
                # Incremental save retains completed points if a later point fails.
                write_json(args.output/"sweep.json",results)
                print(f"V={v:g}, filling={nu:g}, energy={result['energy_per_site']:.9g}",flush=True)
    else:
        model,result=point(args.V,args.filling)
        write_json(args.output/"one_point.json",result)
        save_band_data(args.output/"bands.npz",model,result["currents"])
        import matplotlib
        matplotlib.use("Agg")
        from .plotting import plot_bands
        from matplotlib import pyplot as plt
        fig=plot_bands(model,result["currents"]);fig.savefig(args.output/"bands.png",dpi=160);plt.close(fig)
        print(json.dumps({k:result[k] for k in ["currents","energy_per_site","requested_filling","actual_filling","particle_error","selected"]},indent=2))


if __name__=="__main__":main()
