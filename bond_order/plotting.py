"""Compact LaTeX figures with explicit numerical solver-status markers."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from .publication_style import apply_style, width, panel_colorbar


def flux_mod_2pi(flux):
    """Oriented flux in [0, 2 pi), with 1e-10 rad endpoint roundoff removed."""
    value = np.remainder(flux, 2 * np.pi)
    endpoint = np.isclose(value, 0, atol=1e-10, rtol=0) | np.isclose(
        value, 2 * np.pi, atol=1e-10, rtol=0
    )
    return np.where(endpoint, 0.0, value)


def scan_grid(results):
    """Index a complete rectangular scan by its actual coordinates.

    Input order need not match a reshape convention. Mixed settings,
    duplicates and missing cells are errors, not silently interpolated data.
    """
    if not results:
        raise ValueError("scan contains no points")
    fillings = np.array(sorted({r["requested_filling"] for r in results}))
    interactions = np.array(sorted({r["parameters"]["interaction"] for r in results}))
    if len(fillings) < 2 or len(interactions) < 2:
        raise ValueError(
            "a phase map requires at least two fillings and two interactions"
        )
    parameter_keys = ["mesh", "t", "tp", "occupation", "particle_tolerance"]
    setting_keys = ["samples", "iterations", "symmetry", "zero_flux"]
    signature = lambda r: tuple(r["parameters"][k] for k in parameter_keys) + tuple(
        r["settings"][k] for k in setting_keys
    )
    expected = signature(results[0])
    j = np.full((len(interactions), len(fillings), 4), np.nan)
    energy = np.full(j.shape[:2], np.nan)
    trusted = np.zeros(j.shape[:2], dtype=bool)
    for r in results:
        if signature(r) != expected:
            raise ValueError("scan mixes incompatible calculation settings")
        row = int(np.searchsorted(interactions, r["parameters"]["interaction"]))
        col = int(np.searchsorted(fillings, r["requested_filling"]))
        if np.isfinite(energy[row, col]):
            raise ValueError("duplicate scan coordinate")
        values = np.asarray(r["currents"], dtype=float)
        if (
            values.shape != (4,)
            or not np.isfinite(values).all()
            or not np.isfinite(r["energy_per_site"])
        ):
            raise ValueError("invalid scan result")
        j[row, col] = values
        energy[row, col] = r["energy_per_site"]
        trusted[row, col] = any(
            c.get("success") and c.get("feasible") for c in r["solver_candidates"]
        )
    if not np.isfinite(energy).all():
        raise ValueError("scan has missing cells")
    t = results[0]["parameters"]["t"]
    if t <= 0:
        raise ValueError("phase-map flux convention requires t > 0")
    phi = np.arctan(interactions[:, None, None] * j / t)
    flux = (phi[:, :, 3] - phi[:, :, 0] + phi[:, :, 1] - phi[:, :, 2]) / np.pi
    hopping = np.sqrt(t * t + (interactions[:, None, None] * j) ** 2)
    return {
        "fillings": fillings,
        "interactions": interactions,
        "currents": j,
        "energy": energy,
        "solver_success": trusted,
        "square_flux_over_pi": flux,
        "square_flux_mod_2pi_over_pi": flux_mod_2pi(flux * np.pi) / np.pi,
        "current_rms": np.sqrt(np.mean(j * j, axis=2)),
        "bond_amplitude_spread": np.ptp(hopping, axis=2),
        "additional_seed_refinement": any(
            "seed_source" in c for r in results for c in r["solver_candidates"]
        ),
        "parameters": results[0]["parameters"],
        "settings": results[0]["settings"],
    }


def _edges(values, lower=None, upper=None):
    edges = np.r_[
        values[0] - (values[1] - values[0]) / 2,
        (values[:-1] + values[1:]) / 2,
        values[-1] + (values[-1] - values[-2]) / 2,
    ]
    if lower is not None:
        edges[0] = max(edges[0], lower)
    if upper is not None:
        edges[-1] = min(edges[-1], upper)
    return edges


def _map(ax, grid, values, **kwargs):
    im = ax.pcolormesh(
        _edges(grid["fillings"], 0, 1),
        _edges(grid["interactions"], 0),
        values,
        shading="flat",
        rasterized=True,
        **kwargs,
    )
    ax.set(xlabel=r"$\nu$", ylabel=r"$V/t$")
    ax.set_xlim(grid["fillings"][[0, -1]])
    ax.set_ylim(grid["interactions"][[0, -1]])
    # Convert energy-coordinate labels only when t differs from one.
    if grid["parameters"]["t"] != 1:
        from matplotlib.ticker import FuncFormatter

        t = grid["parameters"]["t"]
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value/t:g}"))
    failed = ~grid["solver_success"]
    if failed.any():
        xx, yy = np.meshgrid(grid["fillings"], grid["interactions"])
        ax.scatter(xx[failed], yy[failed], marker="x", c="black", s=12, linewidths=0.7)
    return im


def _caption(grid):
    p, s = grid["parameters"], grid["settings"]
    caption = (
        f"k mesh {p['mesh']} x {p['mesh']}; t' / t={p['tp']/p['t']:g}; "
        f"{len(grid['fillings'])} x {len(grid['interactions'])} parameter points; "
        f"{p['occupation']} occupation\nSHGO {s['samples']} samples x {s['iterations']} iterations"
    )
    if grid["additional_seed_refinement"]:
        caption += "; admissible-seed refinements (see audit)"
    return caption


def plot_phase_diagram(results):
    """Four order-parameter maps, without guessing published phase labels."""
    apply_style()
    grid = scan_grid(results)
    fig, axes = plt.subplots(
        2, 2, figsize=(width(2), width(2)), constrained_layout=True
    )
    for index, ax in enumerate(axes.flat):
        im = _map(
            ax,
            grid,
            grid["currents"][:, :, index],
            cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
        )
        ax.set_title(f"({chr(97+index)})")
        panel_colorbar(fig, ax, im, rf"$J_{index+1}$")
    subspace = r"$\Phi=0$" if grid["settings"]["zero_flux"] else "Unconstrained"
    fig.suptitle(subspace)
    if not grid["solver_success"].all():
        from matplotlib.lines import Line2D

        fig.legend(
            [Line2D([], [], marker="x", color="black", linestyle="none", markersize=4)],
            ["Solver not confirmed"],
            loc="outside lower center",
            frameon=False,
            handlelength=1,
        )
    return fig


def plot_phase_comparison(ordinary, zero_flux):
    """Compare currents, bond-amplitude anisotropy and plaquette flux."""
    grids = [scan_grid(ordinary), scan_grid(zero_flux)]
    a, b = grids
    if a["settings"]["zero_flux"] or not b["settings"]["zero_flux"]:
        raise ValueError("provide unconstrained and zero-flux scans in that order")
    for key in ["fillings", "interactions"]:
        if not np.array_equal(a[key], b[key]):
            raise ValueError("comparison grids differ")
    for key in ["mesh", "t", "tp", "occupation", "particle_tolerance"]:
        if a["parameters"][key] != b["parameters"][key]:
            raise ValueError("comparison model settings differ")
    apply_style()
    fig, axes = plt.subplots(
        3, 2, figsize=(width(2), width(2) * 1.4), constrained_layout=True
    )
    fields = ["current_rms", "bond_amplitude_spread", "square_flux_mod_2pi_over_pi"]
    labels = [
        r"$J_{\mathrm{rms}}$",
        r"$\Delta |t_{\mathrm{eff}}|$",
        r"$\Phi/\pi$",
    ]
    for row, (field, label) in enumerate(zip(fields, labels)):
        high = max(max(float(np.max(abs(g[field]))) for g in grids), 1e-8)
        for col, g in enumerate(grids):
            kwargs = (
                dict(cmap="twilight", vmin=0, vmax=2)
                if row == 2
                else dict(cmap="viridis", vmin=0, vmax=high)
            )
            im = _map(axes[row, col], g, g[field], **kwargs)
            axes[row, col].set_title(
                f"({chr(97+row*2+col)}) "
                + (("Unconstrained" if col == 0 else r"$\Phi=0$") if row == 0 else "")
            )
            panel_colorbar(fig, axes[row, col], im, label)

    return fig


def plot_bands(model, currents):
    apply_style()
    bands = model.bands(currents)
    fig, axes = plt.subplots(
        1, 2, figsize=(width(2), width(2) * 0.58), constrained_layout=True
    )
    for ax, band, title, symbol in zip(
        axes, bands, ["(a)", "(b)"], [r"$E_-$", r"$E_+$"]
    ):
        im = ax.pcolormesh(
            model.kx / np.pi, model.ky / np.pi, band, shading="auto", cmap="viridis"
        )
        ax.set(xlabel=r"$k_x/\pi$", ylabel=r"$k_y/\pi$", title=title)
        panel_colorbar(fig, ax, im, symbol)
    fig.suptitle(rf"$V/t={model.interaction/model.t:g}$")
    return fig
