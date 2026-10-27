"""Restyle recovered historical author phase regions without inventing labels."""

from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from .publication_style import apply_style, width, save, C

LABELS = ["N", "A", "B", "C", "D", "E", "F", "G", "H"]
COLORS = [
    C["light"],
    C["blue"],
    C["sky"],
    C["orange"],
    C["green"],
    C["yellow"],
    C["purple"],
    C["red"],
    C["gray"],
]


def load_regions(path):
    data = json.loads(Path(path).read_text())
    phase = np.asarray(data["phase_labels"])
    if phase.shape != (len(data["interactions"]), len(data["fillings"])):
        raise ValueError("phase labels must match the actual saved coordinate grid")
    if not set(phase.flat) <= set(LABELS):
        raise ValueError("unknown historical phase label")
    return data


def draw(ax, data, title):
    phase = np.asarray(data["phase_labels"])
    indices = np.vectorize(LABELS.index)(phase)
    ax.pcolormesh(
        data["fillings"],
        data["interactions"],
        indices,
        shading="nearest",
        cmap=ListedColormap(COLORS),
        norm=BoundaryNorm(np.arange(len(LABELS) + 1) - 0.5, len(LABELS)),
        rasterized=True,
    )
    ax.set(
        xlabel=r"$\nu$",
        ylabel=r"$V/t$",
        title=title,
        xlim=(0, 1),
        ylim=(min(data["interactions"]), max(data["interactions"])),
    )
    for label, point in data["original_label_seeds"].items():
        color = "white" if label in ["A", "F", "G", "H"] else C["black"]
        ax.text(
            point["filling"],
            point["interaction"],
            rf"${label}$",
            color=color,
            ha="center",
            va="center",
        )


def export_regions(root):
    root = Path(root)
    apply_style()
    data = [
        load_regions(root / "historical_regions" / f"{name}.json")
        for name in ["unconstrained", "zero_flux"]
    ]
    for item, title in zip(data, ["Unconstrained", r"$\Phi=0$"]):
        fig, ax = plt.subplots(
            figsize=(width(1), width(1) + 0.3), constrained_layout=True
        )
        draw(ax, item, title)
        stem = root / f"phase_regions_{item['subspace']}"
        save(fig, stem)
        plt.close(fig)
    fig, axes = plt.subplots(
        1, 2, figsize=(width(2), width(2) * 0.60), constrained_layout=True
    )
    for ax, item, title in zip(axes, data, ["(a) Unconstrained", r"(b) $\Phi=0$"]):
        draw(ax, item, title)
    fig.legend(
        [Patch(facecolor=c) for c in COLORS],
        [rf"${label}$" for label in LABELS],
        loc="outside lower center",
        ncol=9,
        frameon=False,
        handlelength=0.8,
        handletextpad=0.35,
        columnspacing=0.8,
    )
    fig.suptitle(r"$t=1$, $t^{\prime}=0.2$")
    stem = root / "phase_regions"
    save(fig, stem)
    plt.close(fig)
