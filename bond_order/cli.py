"""Explicit one-point and sweep commands; no calculation is triggered on import."""

import argparse, json
from pathlib import Path
from .model import FluxModel
from .optimization import solve
from .output import write_json, save_band_data


def main():
    parser = argparse.ArgumentParser(
        description="Refactored bond-order research calculations"
    )
    parser.add_argument("command", choices=["one-point", "sweep", "plot-sweep"])
    parser.add_argument("--input", type=Path, help="saved sweep JSON for plot-sweep")
    parser.add_argument(
        "--resume", action="store_true", help="resume a compatible saved sweep"
    )
    parser.add_argument("--mesh", type=int, default=24)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--tp", type=float, default=0.2)
    parser.add_argument("--V", type=float, default=4.0)
    parser.add_argument("--filling", type=float, default=0.4)
    parser.add_argument(
        "--occupation", choices=["legacy", "canonical"], default="legacy"
    )
    parser.add_argument("--particle-tolerance", type=float, default=3.0)
    parser.add_argument("--zero-flux", action="store_true")
    parser.add_argument("--symmetry", choices=["rlex", "legacy_lex"])
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--iterations", type=int, default=1)
    parser.add_argument("--fillings", help="comma-separated fillings for sweep")
    parser.add_argument("--interactions", help="comma-separated V values for sweep")
    parser.add_argument("--output", type=Path, default=Path("examples/output"))
    parser.add_argument(
        "--font-path", type=Path, help="legacy option; LaTeX selects Computer Modern"
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    def plot_scan(results):
        import matplotlib

        matplotlib.use("Agg")
        from .plotting import plot_phase_diagram
        from matplotlib import pyplot as plt

        fig = plot_phase_diagram(results)
        from .publication_style import apply_style, save

        apply_style(font_path=args.font_path)
        save(fig, args.output / "phase_diagram")
        plt.close(fig)

    def point(v, nu):
        model = FluxModel(
            args.mesh, args.t, args.tp, v, args.occupation, args.particle_tolerance
        )
        result = solve(
            model, nu, args.zero_flux, args.samples, args.iterations, args.symmetry
        )
        return model, result

    if args.command == "plot-sweep":
        if not args.input:
            parser.error("plot-sweep requires --input")
        plot_scan(json.loads(args.input.read_text()))
    elif args.command == "sweep":
        if not args.fillings or not args.interactions:
            parser.error("sweep requires --fillings and --interactions")
        fillings = [float(v) for v in args.fillings.split(",")]
        interactions = [float(v) for v in args.interactions.split(",")]
        if len(set(fillings)) != len(fillings) or len(set(interactions)) != len(
            interactions
        ):
            parser.error("sweep coordinates must be unique")
        path = args.output / "sweep.json"
        results = json.loads(path.read_text()) if args.resume and path.exists() else []
        wanted = {(v, nu) for v in interactions for nu in fillings}
        existing = set()
        expected = FluxModel(
            args.mesh, args.t, args.tp, args.V, args.occupation, args.particle_tolerance
        ).parameters()
        symmetry = args.symmetry or ("legacy_lex" if args.zero_flux else "rlex")
        for r in results:
            key = (r["parameters"]["interaction"], r["requested_filling"])
            if key not in wanted or key in existing:
                parser.error("checkpoint coordinates differ from this sweep")
            if any(
                r["parameters"][k] != expected[k]
                for k in expected
                if k != "interaction"
            ):
                parser.error("checkpoint model settings differ")
            if any(
                r["settings"][k] != v
                for k, v in {
                    "samples": args.samples,
                    "iterations": args.iterations,
                    "zero_flux": args.zero_flux,
                    "symmetry": symmetry,
                }.items()
            ):
                parser.error("checkpoint solver settings differ")
            existing.add(key)
        for v in interactions:
            for nu in fillings:
                if (v, nu) in existing:
                    continue
                _, result = point(v, nu)
                results.append(result)
                # Incremental save retains completed points if a later point fails.
                write_json(args.output / "sweep.json", results)
                print(
                    f"V={v:g}, filling={nu:g}, energy={result['energy_per_site']:.9g}",
                    flush=True,
                )
        if len(fillings) > 1 and len(interactions) > 1:
            plot_scan(results)
    else:
        model, result = point(args.V, args.filling)
        write_json(args.output / "one_point.json", result)
        save_band_data(args.output / "bands.npz", model, result["currents"])
        import matplotlib

        matplotlib.use("Agg")
        from .plotting import plot_bands
        from matplotlib import pyplot as plt
        from .publication_style import apply_style, save

        fig = plot_bands(model, result["currents"])
        apply_style(font_path=args.font_path)
        save(fig, args.output / "bands")
        plt.close(fig)
        print(
            json.dumps(
                {
                    k: result[k]
                    for k in [
                        "currents",
                        "energy_per_site",
                        "requested_filling",
                        "actual_filling",
                        "particle_error",
                        "selected",
                    ]
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
