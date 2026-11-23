"""Build masked, shuffled photo cards for junction and near-miss review."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from score import cross, read_csv, references, write_csv


X_REVIEW_CUTOFF = 0.186


def rank(text):
    return hashlib.sha256(text.encode()).digest()


def gap_site(folder, row):
    i = row["id"]
    size = np.array([int(row["width"]), int(row["height"])])
    lines = np.loadtxt(folder / (i + ".annot.csv"), delimiter=",", skiprows=1).reshape(
        -1, 2, 2
    )
    sites = references(lines, size)
    targets = lines[:, 1] - lines[:, 0]
    lengths = np.linalg.norm(targets, axis=1)
    best = None
    for a, stem in enumerate(lines):
        if lengths[a] < 10:
            continue
        for reverse in (0, 1):
            s, e = stem[:: 1 if reverse == 0 else -1]
            axis = (e - s) / lengths[a]
            den = axis[0] * targets[:, 1] - axis[1] * targets[:, 0]
            valid = np.flatnonzero((lengths >= 20) & (abs(den) >= 0.65 * lengths))
            for b in valid:
                diff = lines[b, 0] - e
                gap = cross(diff, targets[b]) / den[b]
                along = cross(diff, axis) / den[b] * lengths[b]
                if not (4 <= gap <= 10 and 8 <= along <= lengths[b] - 8):
                    continue
                q = e + gap * axis
                if (
                    np.any(q < 14)
                    or np.any(q > size - 13)
                    or any(np.linalg.norm(q - p) < 8 for _, p in sites)
                ):
                    continue
                key = rank(f"issue3-control/{i}/{a}/{b}/{reverse}")
                if best is None or key < best[0]:
                    best = (
                        key,
                        dict(id=i, x=q[0], y=q[1], stem=a, target=b, reverse=reverse),
                    )
    return best[1] if best else None


def card(folder, row, gap):
    i = row["id"]
    xy = np.array([float(row["x"]), float(row["y"])]) - 1
    with Image.open(folder / (i + ".small.png")) as source:
        photo = source.convert("RGB")
    panels = []
    for half, zoom in [(26, 3), (18, 4)]:
        left, top = np.rint(xy).astype(int) - half
        length = 2 * half + 1
        tile = photo.crop((left, top, left + length, top + length)).resize(
            (length * zoom, length * zoom), Image.Resampling.NEAREST
        )

        def screen(pt):
            return tuple((np.asarray(pt) - [left, top] + 0.5) * zoom - 0.5)

        draw = ImageDraw.Draw(tile)
        if gap and len(panels):
            lines = (
                np.loadtxt(
                    folder / (i + ".annot.csv"), delimiter=",", skiprows=1
                ).reshape(-1, 2, 2)
                - 1
            )
            for idx, color in [
                (int(row["stem"]), "cyan"),
                (int(row["target"]), "yellow"),
            ]:
                draw.line([screen(p) for p in lines[idx]], fill=color, width=2)
            draw.line(
                [screen(lines[int(row["stem"]), 1 - int(row["reverse"])]), screen(xy)],
                fill="red",
                width=2,
            )
        elif not gap and len(panels):
            x, y = screen(xy)
            r = 4 * zoom
            for a, b in [(-1, -1), (1, -1), (-1, 1), (1, 1)]:
                sx, sy = x + a * r, y + b * r
                draw.line(
                    [(sx - a * 5, sy), (sx, sy), (sx, sy - b * 5)], fill="cyan", width=2
                )
        panels.append(tile)
    card = Image.new("RGB", (318, 183), "white")
    for j, panel in enumerate(panels):
        card.paste(panel, (160 * j, 21))
    ImageDraw.Draw(card).text((5, 1), row["token"], fill="black", font_size=16)
    return card


def render(folder, rows, out, gap=False):
    out.mkdir(parents=True, exist_ok=True)
    for first in range(0, len(rows), 9):
        page = Image.new("RGB", (954, 549), "white")
        for j, row in enumerate(rows[first : first + 9]):
            page.paste(card(folder, row, gap), ((j % 3) * 318, (j // 3) * 183))
        page.save(out / f"sheet-{first//9+1:02}.png")


def gap_pack(folder, out, key):
    rows = []
    for photo in sorted(
        read_csv(folder / "panel.csv"), key=lambda r: rank("blind-gap-v1/" + r["id"])
    ):
        site = gap_site(folder, photo)
        if site is not None:
            rows.append(site)
        if len(rows) == 45:
            break
    for n, r in enumerate(rows, 1):
        r["token"] = f"G{n:02}"
    write_csv(key, rows, list(rows[0]))
    render(folder, rows, out, gap=True)
    print(len(rows), "gap candidates")


def prediction_pack(folder, old, new, refs, out, key):
    rows = []
    for name, events in [("previous", old), ("candidate", new)]:
        for k in ("T", "X"):
            pool = [r for r in events if r["type"] == k]
            rng = np.random.default_rng(8214 + (name == "candidate") * 100 + (k == "X"))
            choice = rng.choice(len(pool), size=24, replace=False)
            for j in choice:
                r = pool[j]
                rows.append(dict(id=r["id"], x=r["x"], y=r["y"], type=k, source=name))
    for k in ("T", "X"):
        pool = [r for r in refs if r["kind"] == k]
        for j in np.random.default_rng(632 + (k == "X")).choice(
            len(pool), size=6, replace=False
        ):
            r = pool[j]
            rows.append(
                dict(id=r["id"], x=r["x"], y=r["y"], type=k, source="reference")
            )
    rng = np.random.default_rng(71309)
    rng.shuffle(rows)
    for n, r in enumerate(rows, 1):
        r["token"] = f"P{n:03}"
    write_csv(key, rows, list(rows[0]))
    render(folder, rows, out)


def score_band_pack(folder, out, key):
    preds = read_csv(folder / "topology-remainder-predictions.csv")
    refs = read_csv(folder / "topology-remainder-sites.csv")
    xs = [r for r in preds if r["method"] == "topology" and r["type"] == "X"]
    rng = np.random.default_rng(536704)
    rows = []
    for label, pool, n in [
        ("upper_X", [r for r in xs if float(r["confidence"]) >= X_REVIEW_CUTOFF], 36),
        ("lower_X", [r for r in xs if float(r["confidence"]) < X_REVIEW_CUTOFF], 36),
        ("reference_T", [r for r in refs if r["kind"] == "T"], 9),
        ("reference_X", [r for r in refs if r["kind"] == "X"], 9),
    ]:
        for j in rng.choice(len(pool), size=n, replace=False):
            r = pool[j]
            rows.append(
                dict(
                    id=r["id"],
                    x=r["x"],
                    y=r["y"],
                    source=label.split("_")[0],
                    type=label.split("_")[1],
                    score=r.get("confidence", ""),
                )
            )
    rng.shuffle(rows)
    for i, r in enumerate(rows, 1):
        r["token"] = f"S{i:02}"
    write_csv(key, rows, list(rows[0]))
    render(folder, rows, out)


def reference_pack(folder, out, key):
    sites = read_csv(folder / "topology-prospective-sites.csv")
    rng = np.random.default_rng(47823)
    rows = []
    for kind, status, size in [("X", "0", 36), ("X", "1", 12), ("T", "1", 12)]:
        pool = [r for r in sites if r["kind"] == kind and r["topology"] == status]
        for index in rng.choice(len(pool), size=size, replace=False):
            r = pool[index]
            rows.append(dict(id=r["id"], x=r["x"], y=r["y"], kind=kind, matched=status))
    rng.shuffle(rows)
    for i, r in enumerate(rows, 1):
        r["token"] = f"V{i:02}"
    write_csv(key, rows, list(rows[0]))
    render(folder, rows, out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("sheets", type=Path)
    parser.add_argument("key", type=Path)
    parser.add_argument(
        "--predictions",
        type=Path,
        help="A prospective score run with previous and topology methods",
    )
    parser.add_argument(
        "--references",
        action="store_true",
        help="Mask a sample of matched and missed reference sites",
    )
    parser.add_argument(
        "--score-bands",
        action="store_true",
        help="Mask the secondary score-band evaluation",
    )
    args = parser.parse_args()
    if args.score_bands:
        score_band_pack(args.folder, args.sheets, args.key)
    elif args.references:
        reference_pack(args.folder, args.sheets, args.key)
    elif args.predictions:
        all_rows = read_csv(args.predictions)
        prediction_pack(
            args.folder,
            [r for r in all_rows if r["method"] == "supported"],
            [r for r in all_rows if r["method"] == "topology"],
            read_csv(args.folder / "topology-prospective-sites.csv"),
            args.sheets,
            args.key,
        )
    else:
        gap_pack(args.folder, args.sheets, args.key)
