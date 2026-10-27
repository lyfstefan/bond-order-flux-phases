# Source audit and refactor boundary (October 2026)

The author confirmed that the project is in
`Research_QC&HighTc/Loop current phase/Numerics/Python Code`, and stated that
the exact final version is not remembered. That uncertainty is preserved.
The separate Affleck project is not the source of this refactor.

| Candidate | t | t' | Mesh | Active method / interpretation |
|---|---:|---:|---|---|
| `phase diagram J1 J2 J3 J4.py` | 1 | 0 | 10x10 | Early small-grid three-LEX-sector SHGO; top-level sweep |
| `phase diagram J1 J2 J3 J4_zero flux subspace.py` | 1 | .2 | 100x100 | Three LEX sectors plus zero-square-flux equality; top-level sweep |
| `parallel calculation.py` | 1 | .2 | 50x50 | Active RLEX, SHGO Sobol n=256, iters=2 |
| `parallel calculation_integral.py` | 1 | .2 | 100x100 | Active RLEX, SHGO Sobol n=256, iters=2; despite filename, energy uses a discrete mesh |

AST comparison finds identical E1/E2 function bodies in all four. Their
expanded radicand is the quadratic quantity |h_AB|^2, and taking its square
root is correct. The square-root-of-a-modulus issue previously found in a
separate project does not describe these bond-order functions.

`reference/parallel_integral.py` extracts E1, E2, Etot and Order from the
100x100 parallel candidate. `reference/zero_flux.py` extracts E1, E2,
particle number, chemical potential, energy and the constraints/Order from
the zero-flux candidate. Function text is retained with only trailing
whitespace removed; imports and required default globals are supplied.
No original script was imported or run at top level. The manifests record
input SHA-256, function AST fingerprints and selected functions. All original
Physics files remain unchanged.

The new `bond_order` interface refactors the shared expressions into explicit
model instances. Legacy occupation reproduces the original bisection and
particle-count tolerance, while computing bands once per energy call.
Canonical occupation is an opt-in numerical change. The new solver preserves
the active RLEX inequalities and zero-flux LEX constraints, but always reports
the zero-flux source's 1e-5 LEX slack (including its shifted current equality),
its sampling budget, feasibility and solver messages. Its explicit Sobol
settings are not a recovery of every implicit default in older versions.

No full phase diagram was rerun. Fixed-trial source comparison and mesh
checks cannot identify which historical optimizer run generated a figure.
The clean package, example curves and numerical reports are a 2026 refactor
release, not the final paper's original output archive. The publication's
loop-current-only scope remains in force.
