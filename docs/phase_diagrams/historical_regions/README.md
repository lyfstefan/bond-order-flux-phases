# Categorical phase-region data

The JSON files store the N/A–H phase regions in the submitted figures for
the four-current model at $t=1$ and $t'=0.2$.

| File | Square-flux constraint | Filling | $V/t$ | Grid |
| --- | --- | --- | --- | --- |
| `unconstrained.json` | None | 0–1 | 0–5 | $20\times20$ |
| `zero_flux.json` | Zero square flux | 0–1 | 2–5 | $20\times20$ |

Each record contains the coordinate grid, phase labels, phase counts and
the SHA-256 checksum of its source figure. The 800 labels were read directly
from the visible vector cells in `Phase_diagram.pdf` and
`Phase_diagram_theta2.pdf`. The color-to-phase correspondence was checked
at the letter annotations in both figures.

These files preserve the plotted phase classification, including the
boundaries chosen in those figures. They do not contain the underlying
optimized currents. The numerical current and flux scans are stored
separately in [../audited/](../audited/).

## Plot the region maps

From the repository root:

```sh
python -m examples.phase_regions
```

This exports the combined phase-region figure and the two individual maps
under `docs/phase_diagrams/`.

To extract the labels again from locally available source figures:

```sh
python -m pip install pymupdf
python -m examples.recover_historical_regions --source /path/to/original/figures
```

The extraction preserves the cell labels and does not perform minimization,
reclassification or boundary smoothing.
