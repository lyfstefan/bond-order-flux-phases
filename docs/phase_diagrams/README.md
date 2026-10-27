# Bond-order phase diagrams

![Historical phase regions](phase_regions.png)

The previously missing categorical overview is now provided from the
original author figure cells. [Recovery provenance](historical_regions/README.md)
distinguishes those historical region labels from the new calculations below.

# New calculated phase-space maps

![Comparison](phase_diagram_comparison.png)

The maps compare unconstrained and zero-square-flux mean-field calculations.
They show current strength, effective-bond-amplitude ordering and square flux.
Four individual current-component maps are also provided below. No A-H phase
letters, interpolated cells or manually edited phase-boundary curves are used.

## Calculation

- t=1, t'=0.2; filling .05,.10,...,.95; V=.25,.50,...,5.
- 19x20=380 points per subspace, 760 independent SHGO calculations in total.
- Momentum mesh 40x40; canonical fixed filling; RLEX; SHGO 256 Sobol samples,
  two iterations. The phase-map occupation choice is explicit; it does not
  silently alter the package's legacy default.
- Every point records currents, requested/actual filling, energy components,
  solver settings, statuses and feasibility. Original Physics sources are unchanged.

## Files

| Artifact | Purpose |
|---|---|
| [Comparison PNG](phase_diagram_comparison.png), [SVG](phase_diagram_comparison.svg) | Current-strength, bond-ordering and flux maps |
| [Unconstrained J1-J4](audited/unconstrained_order_parameters.png) | Four current-component maps after recorded candidate refinement |
| [Zero-flux J1-J4](audited/zero_flux_order_parameters.png) | Corresponding restricted-subspace maps |
| [Raw ordinary scan](unconstrained/sweep.json) | Independent unmodified SHGO results |
| [Raw zero-flux scan](zero_flux/sweep.json) | Independent unmodified SHGO results |
| [Audited ordinary data](audited/unconstrained.json), [NPZ](audited/unconstrained.npz) | Values actually plotted in the audited maps |
| [Audited zero-flux data](audited/zero_flux.json), [NPZ](audited/zero_flux.npz) | Corresponding restricted-subspace values |
| [Numerical audit](audit.json) | Every refinement, variational comparison and selected finer-grid check |

## Candidate refinement and scope

The independent scans initially had 34 points where the ordinary
solution energy exceeded the restricted energy by more than 1e-7, with maximum
violation 0.0010587777. That indicates missed minimizers,
since each restricted solution is also an admissible ordinary trial.

The plotted data therefore include explicit admissible cross-subspace seeds,
selected finer-grid seeds evaluated again at the coarse grid, and forward/backward
neighbor continuation. For transfers between V values, J is rescaled to keep
V*J (and all link phases) fixed; transfers outside the original bounds are rejected.
Every accepted candidate is re-evaluated and locally optimized at its actual
(nu,V). Neither energies nor currents are interpolated into the map. Raw independent
scans remain unchanged. The audit lists 163 accepted candidate updates.
The final maximum variational-nesting residual is 1.05e-12.

Twelve selected checks use a 60x60 mesh, 512 SHGO samples and three iterations,
also retaining refined feasible seeds. Among those checks, the largest fixed-current
mesh energy change is 0.000108821, and the largest extra search gain on the fine mesh
is 0.000140115, in the common energy unit per site. These finite checks do not establish
that every phase boundary is converged. The figures are exploratory finite-grid
maps within the loop-current ansatz; there is no CDW competition or certified
reproduction of the paper's exact A-H phase diagram.

Regenerate the calculations and audited figures from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -m examples.phase_diagrams --output examples/output/phase_diagrams
```

Regenerate from these saved independent scans:

```sh
OPENBLAS_NUM_THREADS=1 python -m examples.phase_diagrams --plot-only --output docs/phase_diagrams
```

The latter repeats the recorded seed/refinement procedure and selected checks,
but does not repeat the 760 independent global searches. `--skip-refinement-checks`
skips the finer-grid checks and their seeds, so it is not identical to this release.
