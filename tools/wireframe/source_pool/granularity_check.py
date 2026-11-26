"""Supplementary check outside the frozen rule: how much M2 rewards short curves.

The artist-guided greedy can pick short pieces that lie on artist lines. This
cuts every FDoG and M-LSD curve into consecutive runs of at most k raster
pixels and reruns M2 on the baseline pool alone. The per-curve exports are
sorted, deduplicated column-major indices, not arc-length order, so a run can
join pixels that are not adjacent along the curve. k = 2 and k = 13 are the median per-curve ink of the
structured-edge and FDoG-512 sources.

    py tools/wireframe/source_pool/granularity_check.py \
        .agent-work/contour-study/draw-toronto .agent-work/source-pool
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from screen import ROSTER, greedy, png_mask  # noqa: E402

SIZES = (2, 13)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("draw", type=Path)
    parser.add_argument("pools", type=Path)
    parser.add_argument(
        "output", type=Path, nargs="?", help="optional CSV of pooled results"
    )
    args = parser.parse_args()
    with ROSTER.open() as f:
        roster = list(csv.DictReader(f))
    groups = ["all", *sorted({r["category"] for r in roster})]
    totals = {(g, k): 0 for g in groups for k in SIZES}
    artists = {g: 0 for g in groups}
    curves = {(g, k): 0 for g in groups for k in SIZES}
    for r in roster:
        i = r["id"]
        artist = png_mask(args.draw / f"{i}.artist.drawing.png")
        h, w = artist.shape
        n = h * w
        af = artist.ravel(order="F")
        for g in ("all", r["category"]):
            artists[g] += int(af.sum())
        cap = int(png_mask(args.draw / f"{i}.fdog.drawing.png").sum())
        pool = json.loads((args.pools / "pools" / f"{i}.json").read_text())
        offsets = [
            (dy, dx)
            for dy in range(-2, 3)
            for dx in range(-2, 3)
            if dx * dx + dy * dy <= 4
        ]

        def dilated(path):
            y, x = path % h, path // h
            pieces = []
            for dy, dx in offsets:
                yy, xx = y + dy, x + dx
                ok = (yy >= 0) & (yy < h) & (xx >= 0) & (xx < w)
                pieces.append(yy[ok] + h * xx[ok])
            reach = np.unique(np.concatenate(pieces))
            return reach[af[reach]]

        for k in SIZES:
            paths = []
            for m in ("fdog", "mlsd"):
                for p in pool[m]:
                    p = np.asarray(p, dtype=int)
                    _, first = np.unique(p, return_index=True)
                    p = p[np.sort(first)]
                    paths.extend(np.unique(p[s : s + k]) for s in range(0, len(p), k))
            seen, covered = greedy(paths, [dilated(p) for p in paths], n, cap)
            assert seen.sum() <= cap
            for g in ("all", r["category"]):
                totals[g, k] += int(covered.sum())
                curves[g, k] += len(paths)
    rows = [
        dict(
            group=g,
            pool="fdog_mlsd",
            max_piece_pixels=k,
            pieces=curves[g, k],
            artist_pixels=artists[g],
            covered=totals[g, k],
            m2_coverage_pct=round(100 * totals[g, k] / artists[g], 2),
        )
        for g in groups
        for k in SIZES
    ]
    for row in rows:
        if row["group"] == "all":
            print(row)
    if args.output:
        with args.output.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
