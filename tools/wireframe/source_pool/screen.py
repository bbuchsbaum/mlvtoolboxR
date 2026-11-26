"""Frozen development screen of additional candidate curve sources (M1-M4).

Generalizes .agent-work/contour-proposal/headroom.py to any set of sources,
reproduces its FDoG, union and artist-guided baselines exactly, adds the
FDoG-512 and structured-edge sources, and applies the screening rule of
docs/evidence/wireframe/source-pool/PLAN.md. Development photos only.

    py tools/wireframe/source_pool/screen.py \
        .agent-work/contour-study/draw-toronto .agent-work/contour-proposal \
        .agent-work/source-pool .agent-work/source-pool/arms/traced-pool-development.csv \
        docs/evidence/wireframe/source-pool
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[3]
ROSTER = ROOT / "docs/evidence/wireframe/curve-rank/development.csv"
NEW = ("fdog512", "sed")
CONFIGS = {
    "baseline": ("fdog", "mlsd"),
    "with_fdog512": ("fdog", "mlsd", "fdog512"),
    "with_sed": ("fdog", "mlsd", "sed"),
    "with_both": ("fdog", "mlsd", "fdog512", "sed"),
}
ARM_COLUMN = {
    "baseline": "combined",
    "with_fdog512": "with_fdog512",
    "with_sed": "with_sed",
    "with_both": "with_both",
}
M2_POINTS = 2.0
M4_CROSSINGS = 5
M4_PHOTOS = 3
# Earlier ledger, pooled over the 60 photos: ink, coverage (%).
LEDGER = {
    "fdog": (134708, 42.65),
    "both_sources": (223878, 58.92),
    "reference_selection": (133548, 49.94),
}


def incidence(lists, n):
    lengths = list(map(len, lists))
    data = np.ones(sum(lengths), dtype=np.int32)
    indices = np.concatenate(lists) if lists else np.zeros(0, int)
    return csr_matrix(
        (data, indices, np.r_[0, np.cumsum(lengths)]), shape=(len(lists), n)
    ).tocsc()


def png_mask(path):
    image = np.asarray(Image.open(path))
    if image.ndim == 3:
        image = image[..., 0]
    return image > 0


def greedy(paths, coverage, n, cap):
    """The earlier headroom rule: best artist gain per new ink pixel, complete curves."""
    ip, cp = incidence(paths, n), incidence(coverage, n)
    cost = np.array(list(map(len, paths)))
    gain = np.array(list(map(len, coverage)))
    seen = np.zeros(n, bool)
    covered = np.zeros(n, bool)
    remaining = cap
    while remaining:
        valid = (cost > 0) & (cost <= remaining) & (gain > 0)
        if not valid.any():
            break
        ratio = np.divide(gain, cost, out=np.zeros(len(paths)), where=valid)
        top = np.flatnonzero(valid & (ratio == ratio.max()))
        j = top[np.argmax(gain[top])]
        ink = paths[j][~seen[paths[j]]]
        refs = coverage[j][~covered[coverage[j]]]
        remaining -= len(ink)
        seen[ink] = True
        covered[refs] = True
        cost -= np.bincount(ip[:, ink].indices, minlength=len(paths))
        gain -= np.bincount(cp[:, refs].indices, minlength=len(paths))
    assert seen.sum() <= cap
    return seen, covered


def ordered(paths, key, n, cap):
    """Photo-only drawing: complete curves by decreasing key until the next exceeds cap."""
    seen = np.zeros(n, bool)
    used = 0
    taken = 0
    for j in np.argsort(-np.asarray(key, float), kind="stable"):
        new = paths[j][~seen[paths[j]]]
        if used + len(new) > cap:
            break
        seen[new] = True
        used += len(new)
        taken += 1
    return seen, taken


def evaluate(r, args):
    i = r["id"]
    draw = args.draw
    artist = png_mask(draw / f"{i}.artist.drawing.png")
    h, w = artist.shape
    n = h * w
    af = artist.ravel(order="F")

    def dilated(path):
        y, x = path % h, path // h
        pieces = []
        for dy in range(-2, 3):
            for dx in range(-2, 3):
                if dx * dx + dy * dy > 4:
                    continue
                yy, xx = y + dy, x + dx
                ok = (yy >= 0) & (yy < h) & (xx >= 0) & (xx < w)
                pieces.append(yy[ok] + h * xx[ok])
        reach = np.unique(np.concatenate(pieces)) if pieces else np.zeros(0, int)
        return reach[af[reach]]

    pool = json.loads((args.pools / "pools" / f"{i}.json").read_text())
    earlier = json.loads((args.proposal / f"{i}.json").read_text())
    for m in ("fdog", "mlsd"):
        assert pool[m] == earlier[m], f"{i}: {m} pool differs from the earlier export"
    masks = {
        "fdog": png_mask(draw / f"{i}.fdog.drawing.png"),
        "mlsd": png_mask(draw / f"{i}.mlsd.drawing.png"),
        "fdog512": png_mask(args.pools / "draw" / f"{i}.fdog512.drawing.png"),
        "sed": png_mask(args.pools / "draw" / f"{i}.sed.drawing.png"),
    }
    paths, cover = {}, {}
    for m in CONFIGS["with_both"]:
        entries = [np.unique(np.asarray(p, dtype=int)) for p in pool[m]]
        rebuilt = np.zeros(n, bool)
        if entries:
            rebuilt[np.concatenate(entries)] = True
        assert np.array_equal(rebuilt, masks[m].ravel(order="F")), (i, m)
        paths[m] = entries
        cover[m] = [dilated(p) for p in entries]
    near = distance_transform_edt(~artist) <= 2

    def scored(v):
        hit = (
            (distance_transform_edt(~v) <= 2) & artist
            if v.any()
            else np.zeros_like(artist)
        )
        return int(hit.sum()), int((near & v).sum()), int(v.sum())

    cap = int(masks["fdog"].sum())
    row = dict(id=i, category=r["category"], artist=int(af.sum()), fdog_ink_budget=cap)
    for m in CONFIGS["with_both"]:
        row[m + "_curves"] = len(paths[m])
    row.update(
        zip(
            ("m3_fdog_covered", "m3_fdog_correct", "m3_fdog_ink"), scored(masks["fdog"])
        )
    )
    for label, sources in CONFIGS.items():
        union = np.zeros_like(artist)
        for m in sources:
            union |= masks[m]
        row.update(
            zip(
                (f"m1_{label}_" + k for k in ("covered", "correct", "ink")),
                scored(union),
            )
        )
        ps = [p for m in sources for p in paths[m]]
        cs = [c for m in sources for c in cover[m]]
        seen, covered = greedy(ps, cs, n, cap)
        sel = seen.reshape((h, w), order="F")
        got = scored(sel)
        assert got[0] == covered.sum()
        row.update(
            zip((f"m2_{label}_" + k for k in ("covered", "correct", "ink")), got)
        )
        # Covered artist pixels that no FDoG or M-LSD curve reaches at any budget.
        if label == "baseline":
            base_reach = ((distance_transform_edt(~union) <= 2) & artist).ravel(
                order="F"
            )
        row[f"m2_{label}_beyond_baseline_union"] = int((covered & ~base_reach).sum())
    keys = {"fdog512": pool["fdog512_length"], "sed": pool["sed_strength"]}
    for m in NEW:
        assert len(keys[m]) == len(paths[m])
        seen, taken = ordered(paths[m], keys[m], n, cap)
        got = scored(seen.reshape((h, w), order="F"))
        row.update(zip((f"m3_{m}_" + k for k in ("covered", "correct", "ink")), got))
        row[f"m3_{m}_curves"] = taken
    return row


def rates(covered, correct, ink, artist):
    recall = covered / artist
    precision = correct / ink if ink else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return recall, precision, f1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("draw", type=Path)
    parser.add_argument("proposal", type=Path)
    parser.add_argument("pools", type=Path)
    parser.add_argument("arms", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with ROSTER.open() as f:
        roster = list(csv.DictReader(f))
    assert len(roster) == 60 and {r["panel"] for r in roster} == {"development"}
    rows = [evaluate(r, args) for r in roster]

    # Baseline reproduction against the earlier ledger, photo by photo.
    with (args.proposal / "headroom.csv").open() as f:
        ledger = {r["id"]: r for r in csv.DictReader(f)}
    assert set(ledger) == {r["id"] for r in rows}
    pairs = {
        "fdog": "m3_fdog",
        "both_sources": "m1_baseline",
        "reference_selection": "m2_baseline",
    }
    reproduction = {}
    for old, new in pairs.items():
        same = all(
            int(ledger[r["id"]][f"{old}_{k}"]) == r[f"{new}_{k}"]
            for r in rows
            for k in ("covered", "correct", "ink")
        )
        ink = sum(r[new + "_ink"] for r in rows)
        cov = (
            100
            * sum(r[new + "_covered"] for r in rows)
            / sum(r["artist"] for r in rows)
        )
        reproduction[old] = dict(
            per_photo_identical=same,
            ink=ink,
            coverage=round(cov, 2),
            ledger_ink=LEDGER[old][0],
            ledger_coverage=LEDGER[old][1],
        )
        assert same and ink == LEDGER[old][0] and round(cov, 2) == LEDGER[old][1], (
            old,
            ink,
            cov,
        )

    # Arm opportunity, per photo.
    with args.arms.open() as f:
        arms = [r for r in csv.DictReader(f) if r["usable"] == "1"]
    assert len(arms) == 84
    for r in rows:
        mine = [a for a in arms if a["id"] == r["id"]]
        r["m4_measurable"] = len(mine)
        for label, column in ARM_COLUMN.items():
            r[f"m4_{label}_four_arm"] = sum(a[column] == "4" for a in mine)
        for m in NEW:
            r[f"m4_{m}_alone_four_arm"] = sum(a[m] == "4" for a in mine)

    with (args.output / "source-pool-photo-scores.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    groups = ["all", *sorted({r["category"] for r in rows})]
    summary = []
    measures = [("M1", label) for label in CONFIGS] + [
        ("M2", label) for label in CONFIGS
    ]
    measures += [("M3", m) for m in ("fdog", *NEW)]
    for g in groups:
        sub = [r for r in rows if g == "all" or r["category"] == g]
        artist = sum(r["artist"] for r in sub)
        for measure, label in measures:
            prefix = f"{measure.lower()}_{label}_"
            covered = sum(r[prefix + "covered"] for r in sub)
            correct = sum(r[prefix + "correct"] for r in sub)
            ink = sum(r[prefix + "ink"] for r in sub)
            recall, precision, f1 = rates(covered, correct, ink, artist)
            summary.append(
                dict(
                    group=g,
                    photos=len(sub),
                    measure=measure,
                    pool=label,
                    artist_pixels=artist,
                    covered=covered,
                    correct_ink=correct,
                    ink=ink,
                    fdog_ink_budget=sum(r["fdog_ink_budget"] for r in sub),
                    coverage_pct=round(100 * recall, 2),
                    precision_pct=round(100 * precision, 2),
                    f1=round(f1, 4),
                    measurable_crossings="",
                    four_arm_crossings="",
                    four_arm_photos="",
                )
            )
        for label in CONFIGS:
            hit = [r for r in sub if r[f"m4_{label}_four_arm"]]
            summary.append(
                dict(
                    group=g,
                    photos=len(sub),
                    measure="M4",
                    pool=label,
                    measurable_crossings=sum(r["m4_measurable"] for r in sub),
                    four_arm_crossings=sum(r[f"m4_{label}_four_arm"] for r in sub),
                    four_arm_photos=len(hit),
                )
            )
    fields = list(summary[0])
    with (args.output / "source-pool-summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(summary)

    # Screening rule, pooled over the 60 photos.
    total = sum(r["artist"] for r in rows)

    def m2(label):
        return 100 * sum(r[f"m2_{label}_covered"] for r in rows) / total

    base_ids = {(a["id"], a["reference"]) for a in arms if a["combined"] == "4"}
    decision = dict(
        protocol="docs/evidence/wireframe/source-pool/PLAN.md at a2ffd64",
        structured_edges=json.loads((args.pools / "prepare-meta.json").read_text()),
        baseline_reproduction=reproduction,
        rule=dict(m2_points=M2_POINTS, m4_crossings=M4_CROSSINGS, m4_photos=M4_PHOTOS),
        baseline=dict(
            m2_coverage_pct=round(m2("baseline"), 4),
            m4_four_arm=len(base_ids),
            m4_photos=len({a["id"] for a in arms if a["combined"] == "4"}),
        ),
        sources={},
    )
    for label in ("with_fdog512", "with_sed", "with_both"):
        column = ARM_COLUMN[label]
        hit = [a for a in arms if a[column] == "4"]
        gained = [a for a in hit if (a["id"], a["reference"]) not in base_ids]
        lost = base_ids - {(a["id"], a["reference"]) for a in hit}
        assert not lost
        delta = m2(label) - m2("baseline")
        m4_pass = (
            len(gained) >= M4_CROSSINGS and len({a["id"] for a in gained}) >= M4_PHOTOS
        )
        entry = dict(
            m2_coverage_pct=round(m2(label), 4),
            m2_gain_points=round(delta, 4),
            m2_beyond_baseline_union_points=round(
                100 * sum(r[f"m2_{label}_beyond_baseline_union"] for r in rows) / total,
                4,
            ),
            m2_pass=delta >= M2_POINTS,
            m4_four_arm=len(hit),
            m4_gained=len(gained),
            m4_gained_photos=len({a["id"] for a in gained}),
            m4_pass=m4_pass,
        )
        if label != "with_both":
            entry["carry_forward"] = entry["m2_pass"] or m4_pass
        else:
            entry[
                "note"
            ] = "reported for information; the rule applies to each source added alone"
        decision["sources"][label] = entry
    (args.output / "screening-decision.json").write_text(
        json.dumps(decision, indent=2) + "\n"
    )

    print(json.dumps(decision, indent=2))
    for s in summary:
        if s["group"] == "all":
            print(s)


if __name__ == "__main__":
    main()
