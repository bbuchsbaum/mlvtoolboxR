"""Registration check: score each source drawing under small integer shifts.

If the working-frame mapping is right, pooled F1 against the artist ink
should peak at zero shift for the new sources, as for the existing ones.

    py tools/wireframe/source_pool/alignment_check.py \
        .agent-work/contour-study/draw-toronto .agent-work/source-pool
"""

import argparse
import csv
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt

ROSTER = (
    Path(__file__).resolve().parents[3]
    / "docs/evidence/wireframe/curve-rank/development.csv"
)


def png_mask(path):
    image = np.asarray(Image.open(path))
    return (image[..., 0] if image.ndim == 3 else image) > 0


def shifted(v, dx, dy):
    out = np.zeros_like(v)
    h, w = v.shape
    ys, yd = (
        (slice(0, h - dy), slice(dy, h))
        if dy >= 0
        else (slice(-dy, h), slice(0, h + dy))
    )
    xs, xd = (
        (slice(0, w - dx), slice(dx, w))
        if dx >= 0
        else (slice(-dx, w), slice(0, w + dx))
    )
    out[yd, xd] = v[ys, xs]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("draw", type=Path)
    parser.add_argument("pools", type=Path)
    parser.add_argument("output", type=Path, nargs="?", help="optional CSV of pooled results")
    args = parser.parse_args()
    with ROSTER.open() as f:
        ids = [r["id"] for r in csv.DictReader(f)]
    shifts = [(dx, dy) for dy in range(-2, 3) for dx in range(-2, 3)]
    files = {
        "fdog": lambda i: args.draw / f"{i}.fdog.drawing.png",
        "mlsd": lambda i: args.draw / f"{i}.mlsd.drawing.png",
        "fdog512": lambda i: args.pools / "draw" / f"{i}.fdog512.drawing.png",
        "sed": lambda i: args.pools / "draw" / f"{i}.sed.drawing.png",
    }
    totals = {(m, s): np.zeros(3) for m in files for s in shifts}
    artist_total = 0
    for i in ids:
        artist = png_mask(args.draw / f"{i}.artist.drawing.png")
        artist_total += artist.sum()
        near = distance_transform_edt(~artist) <= 2
        for m, path in files.items():
            v0 = png_mask(path(i))
            for s in shifts:
                v = shifted(v0, *s)
                if not v.any():
                    continue
                hit = (distance_transform_edt(~v) <= 2) & artist
                totals[m, s] += (hit.sum(), (near & v).sum(), v.sum())
    rows = []
    for m in files:
        f1 = {}
        for s in shifts:
            c, k, n = totals[m, s]
            r, p = c / artist_total, k / n
            f1[s] = 2 * p * r / (p + r)
            rows.append(
                dict(
                    source=m,
                    dx=s[0],
                    dy=s[1],
                    covered=int(c),
                    correct_ink=int(k),
                    ink=int(n),
                    coverage_pct=round(100 * r, 2),
                    precision_pct=round(100 * p, 2),
                    f1=round(f1[s], 4),
                )
            )
        best = max(f1, key=f1.get)
        print(
            m,
            "best shift",
            best,
            "F1",
            round(f1[best], 4),
            "F1 at (0,0)",
            round(f1[0, 0], 4),
        )
    if args.output:
        with args.output.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
