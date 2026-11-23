"""Generate blinded native-pixel cards and independent regions for junction review."""
import argparse
import json
import shutil
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from review_topology import gap_site, rank
from score import read_csv, write_csv
from arm_experiment import sites


def build(folder, predictions, packet, keys):
    packet.mkdir(parents=True, exist_ok=True)
    keys.mkdir(parents=True, exist_ok=True)
    roster = read_csv(folder / "panel.csv")
    photos = {r["id"]: r for r in roster}
    pool = read_csv(predictions)
    rng = np.random.default_rng(517190)
    cards = []
    for old, new in (("X", "X"), ("T", "X"), ("X", "T"), ("T", "T")):
        group = [r for r in pool if r["context"] == old and r["native"] == new]
        n = min(20, len(group))
        for j in rng.choice(len(group), n, replace=False):
            cards.append(
                dict(group[j], stratum=old + new, population=len(group), sample=n)
            )
    controls = {"T": [], "X": [], "gap": []}
    for p in roster:
        i = p["id"]
        for k, pts in zip("TX", sites(folder, p)):
            controls[k].extend(dict(id=i, x=x, y=y) for x, y in pts)
        g = gap_site(folder, p)
        if g:
            controls["gap"].append({k: g[k] for k in ("id", "x", "y")})
    for k, n in (("T", 8), ("X", 8), ("gap", 12)):
        group = controls[k]
        count = min(n, len(group))
        for j in rng.choice(len(group), count, replace=False):
            cards.append(
                dict(group[j], stratum="ref_" + k, population=len(group), sample=count)
            )
    rng.shuffle(cards)
    pages = []

    def from_working(q, small, native):
        return (np.asarray(q) - 0.5) * native / small - 0.5

    for j, r in enumerate(cards):
        i = r["id"]
        p = photos[i]
        token = f"C{j+1:03}"
        r["token"] = token
        small = np.array([int(p["width"]), int(p["height"])])
        with Image.open(folder / (i + ".png")) as im:
            size = np.array(im.size)
            x, y = from_working([float(r["x"]), float(r["y"])], small, size)
            card = Image.new("RGB", (342, 191), "white")
            for left, half, zoom, mark in [(1, 83, 1, False), (170, 41, 2, True)]:
                a, b = np.rint([x, y]).astype(int) - half
                patch = im.crop((a, b, a + 2 * half + 1, b + 2 * half + 1))
                patch = patch.resize(
                    (patch.width * zoom, patch.height * zoom), Image.Resampling.NEAREST
                )
                if mark:
                    draw = ImageDraw.Draw(patch)
                    cx = (x - a + 0.5) * zoom - 0.5
                    cy = (y - b + 0.5) * zoom - 0.5
                    rx, ry = 4 * size / small * zoom
                    for u in (-1, 1):
                        for v in (-1, 1):
                            sx, sy = cx + u * rx, cy + v * ry
                            draw.line(
                                [(sx - u * 6, sy), (sx, sy), (sx, sy - v * 6)],
                                fill="cyan",
                                width=2,
                            )
                card.paste(patch, (left, 22))
        ImageDraw.Draw(card).text((5, 2), token, fill="black", font_size=16)
        if j % 9 == 0:
            pages.append(Image.new("RGB", (1026, 573), "white"))
        pages[-1].paste(card, ((j % 3) * 342, (j % 9 // 3) * 191))
    for j, pg in enumerate(pages):
        pg.save(packet / f"cards-{j+1:02}.png")
    fields = [
        "token",
        "id",
        "x",
        "y",
        "stratum",
        "population",
        "sample",
        "context",
        "working",
        "scores",
        "native",
    ]
    write_csv(keys / "cards.csv", cards, fields)
    regions = []
    for j, p in enumerate(
        sorted(roster, key=lambda r: rank("native-arm-region/" + r["id"]))[:24]
    ):
        i = p["id"]
        token = f"R{j+1:02}"
        with Image.open(folder / (i + ".png")) as im:
            w, h = im.size
            x0 = int(rng.integers(0, w - 240 + 1))
            y0 = int(rng.integers(0, h - 200 + 1))
            tile = im.crop((x0, y0, x0 + 240, y0 + 200))
            tile.save(packet / (token + ".png"))
        regions.append(
            dict(
                token=token,
                id=i,
                left=x0,
                top=y0,
                width=w,
                height=h,
                small_width=int(p["width"]),
                small_height=int(p["height"]),
            )
        )
    write_csv(keys / "regions.csv", regions, list(regions[0]))
    config = (
        "window.reviewConfig="
        + json.dumps(dict(cards=len(cards), regions=len(regions)))
        + ";\n"
    )
    (packet / "config.js").write_text(config)
    shutil.copyfile(Path(__file__).with_name("arm_review.html"), packet / "index.html")
    print(len(cards), "cards,", len(regions), "independent regions", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("folder", "predictions", "packet", "keys"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    build(args.folder, args.predictions, args.packet, args.keys)
