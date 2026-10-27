# Recovered original phase-region maps

The author-submitted figures `Phase_diagram.pdf` and
`Phase_diagram_theta2.pdf` each contain a visible 20x20 grid of vector
rectangles. All 800 cells were recovered directly from those PDF objects,
without raster sampling, interpolation, new minimization or new phase
classification. Palette-to-letter mappings were independently verified
at every original on-map A-H/N annotation.

The JSON files record source SHA-256 fingerprints, exact phase labels,
cell counts and the recovered coordinate grid. The unconstrained figure
uses filling 0..1 and V/t 0..5; the zero-flux figure uses filling 0..1 and
V/t 2..5. The original zero-flux script independently confirms this range.
Parameters are t=1 and t'=0.2, as recorded by the accompanying manuscript.

These are historical *plotted classifications*. The old research script
notes that boundaries needed hand adjustments after inspecting currents.
Recovering its figure cells does not recover or certify the historical
optimizer outputs. No cells were adjusted or smoothed during this recovery.
They are deliberately separate from the 2026 finite-mesh optimization
results in `../audited/`, which have no automatically assigned A-H labels.

`phase_regions.png/pdf/svg` compares the two subspaces. Individual compact
square maps are provided as `phase_regions_unconstrained.*` and
`phase_regions_zero_flux.*`. Labels and colors have the same correspondence
in all three exports; all text is rendered by LaTeX.

To repeat the recovery with the original author figures available locally,
install the optional `pymupdf` dependency and run:

```sh
python -m examples.recover_historical_regions --source /path/to/original/figures
```

No original absolute filesystem paths or original figure PDFs are bundled.
