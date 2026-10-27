"""Calculate, audit and plot paired phase-space scans. No paper labels guessed."""

import argparse, copy, json, subprocess, sys, platform
from pathlib import Path
import numpy as np
import scipy
from bond_order import FluxModel, solve
from bond_order.optimization import refine_from_seed
from bond_order.output import write_json
from bond_order.publication_style import apply_style, save


def key(row):
    return row["parameters"]["interaction"], row["requested_filling"]


def successful(row):
    return any(c["success"] and c["feasible"] for c in row["solver_candidates"])


def update_candidate(row, seed, source, zero_flux):
    model = FluxModel(**row["parameters"])
    result = refine_from_seed(
        model, row["requested_filling"], seed["currents"], zero_flux
    )
    if result["energy_per_site"] >= row["energy_per_site"] - 1e-10 and successful(row):
        return None
    if result["energy_per_site"] > row["energy_per_site"] + 1e-12:
        return None
    change = {
        "interaction": key(row)[0],
        "filling": key(row)[1],
        "subspace": "zero_flux" if zero_flux else "unconstrained",
        "previous_energy": row["energy_per_site"],
        "new_energy": result["energy_per_site"],
        "source": source,
        "local_success": result["local_success"],
        "local_message": result["local_message"],
    }
    row.update(
        {
            k: result[k]
            for k in [
                "currents",
                "energy_per_site",
                "chemical_potential",
                "electronic_energy_per_site",
                "current_cost_per_site",
                "actual_filling",
                "particle_error",
                "bisection_converged",
                "square_flux",
            ]
        }
    )
    row["selected"] = "admissible seed and local refinement"
    row["solver_candidates"].append(
        {
            "sector": "paired_scan_seed",
            "method": "SLSQP plus retained feasible seed",
            "success": result["local_success"] or successful(seed),
            "feasible": True,
            "currents": result["currents"],
            "energy": result["energy_per_site"],
            "seed_source": source,
            "local_success": result["local_success"],
        }
    )
    row["any_feasible_solver_candidate"] = True
    return change


def assemble(output, refinement_checks=True, font_path=None):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from bond_order.plotting import scan_grid, plot_phase_diagram, plot_phase_comparison

    raw = {
        name: json.loads((output / name / "sweep.json").read_text())
        for name in ["unconstrained", "zero_flux"]
    }
    # Validate coordinates/settings before comparing any energies.
    a, b = [scan_grid(raw[n]) for n in ["unconstrained", "zero_flux"]]
    assert np.array_equal(a["fillings"], b["fillings"]) and np.array_equal(
        a["interactions"], b["interactions"]
    )
    scans = copy.deepcopy(raw)
    ordinary = {key(r): r for r in scans["unconstrained"]}
    zero = {key(r): r for r in scans["zero_flux"]}
    before = [
        ordinary[k]["energy_per_site"] - zero[k]["energy_per_site"] for k in ordinary
    ]
    changes = []
    for coordinate in ordinary:
        general, restricted = ordinary[coordinate], zero[coordinate]
        # An ordinary solution with zero flux is also admissible in the
        # restricted problem. No nonzero-flux seed is ever projected by hand.
        if (
            abs(general["square_flux"]) <= 1e-7
            and general["energy_per_site"] < restricted["energy_per_site"] - 1e-10
        ):
            change = update_candidate(
                restricted, general, "unconstrained scan, already zero flux", True
            )
            if change:
                changes.append(change)
        # Every restricted solution is an admissible unrestricted trial.
        if general["energy_per_site"] > restricted[
            "energy_per_site"
        ] + 1e-10 or not successful(general):
            change = update_candidate(general, restricted, "zero-flux SHGO scan", False)
            if change:
                changes.append(change)
    after = [
        ordinary[k]["energy_per_site"] - zero[k]["energy_per_site"] for k in ordinary
    ]
    assert max(after) < 1e-8
    derived = output / "audited"
    derived.mkdir(parents=True, exist_ok=True)
    audit = {
        "method": "Independent SHGO scans; admissible cross-subspace seeds; selected 60x60 seeds reevaluated at 40x40; forward/backward phase-preserving neighbor seeds and SLSQP. Raw scans preserved unchanged.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "points_per_subspace": len(ordinary),
        "initial_unrestricted_energy_above_restricted_gt_1e_minus7": sum(
            v > 1e-7 for v in before
        ),
        "max_initial_variational_violation": max(before),
        "max_final_variational_violation": max(after),
        "candidate_updates": changes,
        "selected_point_refinement_checks": [],
    }
    # Use independently refined, finer-grid solutions as additional *trials*
    # at the original grid, never as interpolated output values.
    refinement_points = [
        (0.3, 2.0),
        (0.5, 1.5),
        (0.4, 4.0),
        (0.5, 4.0),
        (0.7, 4.5),
        (0.9, 4.75),
    ]
    if refinement_checks:
        for filling, v in refinement_points:
            for name, constrained, rows in [
                ("unconstrained", False, ordinary),
                ("zero_flux", True, zero),
            ]:
                coarse = rows[v, filling]
                fine_model = FluxModel(**{**coarse["parameters"], "mesh": 60})
                fine = solve(
                    fine_model,
                    filling,
                    zero_flux=constrained,
                    symmetry="rlex",
                    samples=512,
                    iterations=3,
                )
                seeded = refine_from_seed(
                    fine_model, filling, coarse["currents"], constrained
                )
                seed = (
                    fine
                    if fine["energy_per_site"] <= seeded["energy_per_site"]
                    else {**coarse, **seeded}
                )
                change = update_candidate(
                    coarse,
                    seed,
                    "independent 60x60 refinement seed, reevaluated at 40x40",
                    constrained,
                )
                if change:
                    changes.append(change)
    # Continue neighboring minima in both directions. At a different V,
    # rescale J to preserve V*J and hence all link phases; reject transfers
    # outside the original box. Every retained value is re-minimized at its
    # actual (nu,V), and every objective change is recorded.
    for name, constrained, rows in [
        ("unconstrained", False, ordinary),
        ("zero_flux", True, zero),
    ]:
        for reverse in [False, True]:
            fs = list(a["fillings"])
            vs = list(a["interactions"])
            if reverse:
                fs.reverse()
                vs.reverse()
            count = 0
            for vi, v in enumerate(vs):
                for fi, filling in enumerate(fs):
                    neighbors = []
                    if fi:
                        neighbors.append(rows[v, fs[fi - 1]])
                    if vi:
                        neighbors.append(rows[vs[vi - 1], filling])
                    for neighbor in neighbors:
                        seed = copy.deepcopy(neighbor)
                        seed["currents"] = (
                            np.asarray(seed["currents"])
                            * neighbor["parameters"]["interaction"]
                            / v
                        ).tolist()
                        if max(abs(np.asarray(seed["currents"]))) > 1:
                            continue
                        change = update_candidate(
                            rows[v, filling],
                            seed,
                            "neighbor continuation at target nu,V",
                            constrained,
                        )
                        if change:
                            changes.append(change)
                            count += 1
            print(
                f"{name} {'backward' if reverse else 'forward'} continuation: {count} energy improvements",
                flush=True,
            )
    # Restore the variational nesting relation after independent continuation.
    for coordinate in ordinary:
        general, restricted = ordinary[coordinate], zero[coordinate]
        if (
            abs(general["square_flux"]) <= 1e-7
            and general["energy_per_site"] < restricted["energy_per_site"] - 1e-10
        ):
            change = update_candidate(
                restricted,
                general,
                "unconstrained continuation, already zero flux",
                True,
            )
            if change:
                changes.append(change)
        if general["energy_per_site"] > restricted["energy_per_site"] + 1e-10:
            change = update_candidate(
                general, restricted, "zero-flux continuation scan", False
            )
            if change:
                changes.append(change)
    after = [
        ordinary[k]["energy_per_site"] - zero[k]["energy_per_site"] for k in ordinary
    ]
    assert max(after) < 1e-8
    audit["max_final_variational_violation"] = max(after)
    for name, rows in scans.items():
        grid = scan_grid(rows)
        write_json(derived / f"{name}.json", rows)
        fig = plot_phase_diagram(rows)
        apply_style(font_path=font_path)
        save(fig, derived / f"{name}_order_parameters")
        fig.savefig(derived / f"{name}_order_parameters.svg")
        plt.close(fig)
        audit[name] = {
            "failed_solver_cells": int(np.count_nonzero(~grid["solver_success"])),
            "max_absolute_square_flux": float(
                np.max(abs(grid["square_flux_over_pi"])) * np.pi
            ),
            "parameters": grid["parameters"],
            "settings": grid["settings"],
        }
        np.savez_compressed(
            derived / f"{name}.npz",
            **{k: v for k, v in grid.items() if isinstance(v, np.ndarray)},
        )
    fig = plot_phase_comparison(scans["unconstrained"], scans["zero_flux"])
    apply_style(font_path=font_path)
    save(fig, output / "phase_diagram_comparison")
    fig.savefig(output / "phase_diagram_comparison.svg")
    plt.close(fig)
    if refinement_checks:
        for filling, v in refinement_points:
            for name, constrained in [("unconstrained", False), ("zero_flux", True)]:
                coarse = (ordinary if not constrained else zero)[v, filling]
                parameters = {**coarse["parameters"], "mesh": 60}
                fine_model = FluxModel(**parameters)
                fixed = fine_model.energy(coarse["currents"], filling)
                fine = solve(
                    fine_model,
                    filling,
                    zero_flux=constrained,
                    symmetry="rlex",
                    samples=512,
                    iterations=3,
                )
                seeded = refine_from_seed(
                    fine_model, filling, coarse["currents"], constrained
                )
                best = min(fine["energy_per_site"], seeded["energy_per_site"])
                audit["selected_point_refinement_checks"].append(
                    {
                        "filling": filling,
                        "interaction": v,
                        "subspace": name,
                        "coarse_mesh": 40,
                        "fine_mesh": 60,
                        "fine_SHGO_samples": 512,
                        "fine_SHGO_iterations": 3,
                        "coarse_energy": coarse["energy_per_site"],
                        "fine_fixed_current_energy": fixed,
                        "fine_best_energy": best,
                        "fixed_current_mesh_delta": fixed - coarse["energy_per_site"],
                        "search_gain_on_fine_mesh": fixed - best,
                        "fine_vs_coarse_energy_delta": best - coarse["energy_per_site"],
                        "fine_SHGO_solver_success": successful(fine),
                        "seeded_local_success": seeded["local_success"],
                    }
                )
                write_json(output / "audit.json", audit)
    audit["scope"] = (
        "Finite-grid order-parameter/flux maps, not certified A-H phase boundaries; selected refinement checks do not establish full-grid convergence. Loop-current ansatz only."
    )
    write_json(output / "audit.json", audit)
    print(
        json.dumps(
            {
                "points_per_subspace": len(ordinary),
                "candidate_updates": len(changes),
                "max_final_variational_violation": max(after),
                "refinement_checks": len(audit["selected_point_refinement_checks"]),
                "comparison_figure": str(output / "phase_diagram_comparison.png"),
            },
            indent=2,
        )
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=Path("examples/output/phase_diagrams")
    )
    parser.add_argument(
        "--plot-only", action="store_true", help="reuse existing independent raw scans"
    )
    parser.add_argument("--skip-refinement-checks", action="store_true")
    parser.add_argument(
        "--font-path", type=Path, help="optional Roman font file to pin at export"
    )
    args = parser.parse_args()
    output = args.output.resolve()
    if not args.plot_only:
        fillings = ",".join(f"{x:.2f}" for x in np.arange(1, 20) * 0.05)
        interactions = ",".join(f"{x:.2f}" for x in np.arange(1, 21) * 0.25)
        for name, zero in [("unconstrained", False), ("zero_flux", True)]:
            cmd = [
                sys.executable,
                "-m",
                "bond_order.cli",
                "sweep",
                "--fillings",
                fillings,
                "--interactions",
                interactions,
                "--occupation",
                "canonical",
                "--mesh",
                "40",
                "--samples",
                "256",
                "--iterations",
                "2",
                "--symmetry",
                "rlex",
                "--resume",
                "--output",
                str(output / name),
            ]
            if zero:
                cmd.append("--zero-flux")
            subprocess.run(cmd, check=True)
    assemble(output, not args.skip_refinement_checks, font_path=args.font_path)


if __name__ == "__main__":
    main()
