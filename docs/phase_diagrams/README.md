# Phase diagrams and order-parameter maps

This directory contains categorical phase-region maps and numerical scans
of the four-current mean-field ansatz. All fillings are particles per
physical site; energies are measured in units of $t$.

## Phase regions

![Categorical phase regions](phase_regions.png)

The N/A–H labels are the classifications shown in the submitted figures.
They are stored as two $20\times20$ grids in
[historical_regions/](historical_regions/README.md). The grids cover
$0\leq\nu\leq1$, with $0\leq V/t\leq5$ for unconstrained square flux and
$2\leq V/t\leq5$ for zero square flux. Both use $t'/t=0.2$.

## Current, bond-amplitude and flux maps

![Numerical comparison of the two subspaces](phase_diagram_comparison.png)

The numerical scans use $t=1$, $t'=0.2$, fillings
$\nu=0.05,0.10,\ldots,0.95$ and interactions
$V=0.25,0.50,\ldots,5$. There are 380 points per flux subspace.
Each point uses a $40\times40$ momentum mesh, canonical occupation, RLEX
symmetry sectors, and SHGO with 256 Sobol samples and two iterations.

The comparison shows the root-mean-square current, effective-bond-amplitude
spread and square flux. The individual maps show $J_1,\ldots,J_4$. Their
cells are calculated parameter points and are not interpolated.

| Files | Contents |
| --- | --- |
| `unconstrained/sweep.json` | Independent unconstrained searches |
| `zero_flux/sweep.json` | Independent zero-flux searches |
| `audited/unconstrained.json`, `.npz` | Unconstrained values used in the comparison |
| `audited/zero_flux.json`, `.npz` | Zero-flux values used in the comparison |
| [Unconstrained currents](audited/unconstrained_order_parameters.png) | Four refined current-component maps |
| [Zero-flux currents](audited/zero_flux_order_parameters.png) | Four refined current-component maps |
| [audit.json](audit.json) | Candidate updates, variational comparison and mesh checks |

## Minimization and convergence

Since the zero-flux variational space is contained in the unconstrained
space, its minimum energy cannot be lower at the same parameters. The
independent scans had 34 violations larger than $10^{-7}$, with a maximum
of $1.0588\times10^{-3}$ per site. Feasible cross-subspace seeds,
selected finer-grid seeds and neighbor continuation were therefore used
for additional local searches at the actual parameter points.

The audit records 163 accepted updates. After refinement, the largest
variational-nesting residual is $1.05\times10^{-12}$. Twelve selected
checks use a $60\times60$ mesh, 512 SHGO samples and three iterations.
Among those checks, the largest fixed-current mesh-energy change is
$1.08821\times10^{-4}$ per site, and the largest additional search gain
is $1.40115\times10^{-4}$ per site.

These checks assess the stated finite-grid calculations within the
loop-current ansatz. Phase-boundary convergence and competition with other
orders require separate calculations. A cross in a raw map indicates that
no saved candidate reported both solver success and feasibility. The two
such points in the raw unconstrained scan, at $\nu=0.1$ and
$V/t=3.0,3.25$, are resolved in the refined data.

## Reproduction

From the repository root, redraw the saved results without minimization:

```sh
python -m examples.replot_saved
```

To calculate the paired scans and perform the refinement checks:

```sh
OPENBLAS_NUM_THREADS=1 python -m examples.phase_diagrams
```

To repeat the refinement procedure from the independent scans in this
directory, without repeating their global searches:

```sh
OPENBLAS_NUM_THREADS=1 python -m examples.phase_diagrams --plot-only --output docs/phase_diagrams
```

PDF, PNG and SVG exports share the same coordinate grids and normalization.
