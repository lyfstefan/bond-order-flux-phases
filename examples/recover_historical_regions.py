"""Recover the existing author's categorical cells from supplied vector PDFs.

Optional recovery dependency: PyMuPDF. Replotting saved JSON does not need it.
"""

from pathlib import Path
import argparse, hashlib, json
from collections import Counter
import numpy as np
from matplotlib import colormaps


def recover(source_directory, output_directory):
    import fitz

    SOURCE = Path(source_directory)
    OUT = Path(output_directory)
    OUT.mkdir(parents=True, exist_ok=True)
    plans = [
        (
            "unconstrained",
            "Phase_diagram.pdf",
            0,
            5,
            {0: "N", 1: "A", 2: "B", 2.5: "F", 3.5: "E", 4.5: "G", 6: "D"},
        ),
        (
            "zero_flux",
            "Phase_diagram_theta2.pdf",
            2,
            5,
            {0: "N", 1: "A", 2: "B", 4: "H", 5: "C", 6: "D"},
        ),
    ]
    for key, filename, vmin, vmax, codes in plans:
        path = SOURCE / filename
        with fitz.open(path) as doc:
            page = doc[0]
            cells = [
                v
                for v in page.get_drawings()
                if v["fill"] is not None
                and 5 < v["rect"].width < 25
                and 5 < v["rect"].height < 25
                and v["rect"].x0 > 70
                and v["rect"].x1 < 365
                and v["rect"].y0 > 0
                and v["rect"].y1 < 330
                and len(v["items"]) == 1
                and v["items"][0][0] == "re"
            ]
            assert len(cells) == 400
            xs = sorted({round((v["rect"].x0 + v["rect"].x1) / 2, 3) for v in cells})
            ys = sorted(
                {round((v["rect"].y0 + v["rect"].y1) / 2, 3) for v in cells},
                reverse=True,
            )
            assert len(xs) == len(ys) == 20
            palette = {
                value: np.array(colormaps["viridis"](value / 6)[:3]) for value in codes
            }
            labels = np.full((20, 20), "", dtype="<U1")
            occupied = set()
            largest = 0
            for cell in cells:
                x = round((cell["rect"].x0 + cell["rect"].x1) / 2, 3)
                y = round((cell["rect"].y0 + cell["rect"].y1) / 2, 3)
                i, j = ys.index(y), xs.index(x)
                assert (i, j) not in occupied
                occupied.add((i, j))
                value = min(
                    codes,
                    key=lambda value: np.linalg.norm(palette[value] - cell["fill"]),
                )
                error = float(np.max(abs(palette[value] - cell["fill"])))
                largest = max(largest, error)
                assert error < 1e-4, (value, cell["fill"])
                labels[i, j] = codes[value]
            # Verify phase-name/color correspondence at the original on-map letter centers.
            seeds = {}
            for b in page.get_text("rawdict")["blocks"]:
                for line in b.get("lines", []):
                    for span in line["spans"]:
                        if span["size"] < 30:
                            continue
                        for ch in span.get("chars", []):
                            if ch["c"] not in set(codes.values()):
                                continue
                            box = fitz.Rect(ch["bbox"])
                            center = (box.tl + box.br) / 2
                            if not (
                                min(xs) < center.x < max(xs)
                                and min(ys) < center.y < max(ys)
                            ):
                                continue
                            i = min(range(20), key=lambda i: abs(ys[i] - center.y))
                            j = min(range(20), key=lambda j: abs(xs[j] - center.x))
                            assert labels[i, j] == ch["c"], (key, ch["c"], labels[i, j])
                            seeds[ch["c"]] = {
                                "filling": float(np.linspace(0, 1, 20)[j]),
                                "interaction": float(np.linspace(vmin, vmax, 20)[i]),
                            }
            assert set(seeds) == set(codes.values()), (key, seeds)
            report = {
                "source_figure": filename,
                "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "provenance": "Recovered from existing author-submitted vector figure cells; historical region labels, not newly classified 2026 optimization results",
                "extraction": "400 visible filled rectangles, unique 20x20 coordinates; categorical color mapping verified at every original phase-letter annotation",
                "scope": "Preserves historical plotted phase regions only. The original script notes hand-adjusted boundaries; underlying historical optimizer outputs are not recovered here. Not a certification of phase boundaries.",
                "parameters": {"t": 1.0, "tp": 0.2},
                "subspace": key,
                "fillings": np.linspace(0, 1, 20).tolist(),
                "interactions": np.linspace(vmin, vmax, 20).tolist(),
                "phase_labels": labels.tolist(),
                "phase_counts": dict(Counter(labels.flat)),
                "original_label_seeds": seeds,
                "max_source_palette_component_error": largest,
                "visible_cells": 400,
            }
            (OUT / f"{key}.json").write_text(json.dumps(report, indent=2) + "\n")
            print(key, report["phase_counts"], "all original labels verified")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        type=Path,
        required=True,
        help="directory containing the two original author figure PDFs",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("docs/phase_diagrams/historical_regions")
    )
    args = parser.parse_args()
    recover(args.source, args.output)


if __name__ == "__main__":
    main()
