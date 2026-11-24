"""Compare atomic crossing selection against the frozen whole-contour selector."""

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates
from threadpoolctl import threadpool_limits

from curve_rank import (
    EVIDENCE,
    SOURCES,
    Selection,
    best,
    emit,
    fit,
    fill_singles,
    label_junctions,
    mask,
    predict_paths,
    read_pool,
    roster,
    training_rows,
)
from score import read_csv, write_csv


@dataclass
class Crossing:
    source: str
    index: int  # Within this source's X sites in export order.
    members: np.ndarray
    pixels: np.ndarray
    tube: np.ndarray
    support: float


def arm_support(rgb, site):
    """Weakest two-sided RGB contrast, discounted for nearly coincident rays."""
    angles = -np.deg2rad(np.array(site["arms"], dtype=float))
    if angles.shape != (4,) or not np.all(np.isfinite(angles)):
        return 0.0
    ordered = np.sort(angles % (2 * np.pi))
    gap = np.diff(ordered, append=ordered[0] + 2 * np.pi).min()
    if gap <= 1e-12:
        return 0.0
    directions = np.column_stack((np.cos(angles), np.sin(angles)))
    normals = np.column_stack((-directions[:, 1], directions[:, 0]))
    center = np.array([site["x"], site["y"]]) - 1
    samples = center + directions[:, None, :] * np.array([3, 5, 7, 9])[None, :, None]
    sides = []
    for delta in (-1.5, 1.5):
        xy = samples + delta * normals[:, None, :]
        yx = xy[:, :, ::-1].transpose(2, 0, 1)
        sides.append(
            np.stack(
                [
                    map_coordinates(
                        channel, yx, output=float, order=1, mode="constant", cval=np.nan
                    )
                    for channel in rgb.transpose(2, 0, 1)
                ],
                axis=-1,
            )
        )
    contrasts = np.nan_to_num(np.abs(sides[1] - sides[0]).mean(axis=-1))
    weakest = np.median(contrasts, axis=1).min()
    return float(np.sin(gap) * weakest / (weakest + 12))


def read_crossings(staged, base, name, pool):
    with Image.open(staged / (name + ".small.png")) as image:
        rgb = np.array(image.convert("RGB"))
    crossings = []
    offset = 0
    for source in SOURCES:
        drawing = json.loads((base / f"{name}.{source}.contours.json").read_text())
        contours = drawing["contours"]
        ids = {c["id"]: i + offset for i, c in enumerate(contours)}
        for index, site in enumerate(s for s in drawing["sites"] if s["type"] == "X"):
            members = np.array(
                sorted({ids[m["contour_id"]] for m in site["members"]}), dtype=int
            )
            if not len(members):
                continue
            pixels = np.unique(np.concatenate([pool.pixels[j] for j in members]))
            area = np.unique(np.concatenate([pool.tubes[j] for j in members]))
            support = arm_support(rgb, site) if len(members) >= 2 else 0.0
            crossings.append(Crossing(source, index, members, pixels, area, support))
        offset += len(contours)
    if offset != len(pool.paths):
        raise ValueError("Crossing and contour pools differ")
    return crossings


def take_group(state, members, selected):
    for j in members:
        if not selected[j]:
            state.take(j, allow_covered=True)
            selected[j] = True


def select_groups(pool, model, crossings):
    state = Selection(pool)
    selected = np.zeros(len(pool.paths), dtype=bool)
    supported = [g for g in crossings if g.support > 0]
    actions = 0
    while True:
        scale = state.scale()
        affordable = state.eligible()
        # Every member of an affordable group is individually affordable.
        density = np.zeros(len(pool.paths))
        if affordable.any():
            density[affordable] = np.maximum(
                0, model.predict(state.features()[affordable])
            )
            rewards = density * scale
            j = best(state, rewards, affordable)
            chosen = np.array([j])
            key = (rewards[j] / state.ink[j], rewards[j])
            is_group = False
        else:
            chosen = None
            key = (-np.inf, -np.inf)
            is_group = False
        for group in supported:
            missing = group.members[~selected[group.members]]
            if not len(missing):
                continue
            cost = int(np.count_nonzero(~state.used[group.pixels]))
            if cost > state.remaining:
                continue
            area = np.count_nonzero(~state.covered[group.tube])
            weights = scale[missing]
            utility = 28 * group.support
            if weights.sum():
                utility += max(cost, area) * np.average(
                    density[missing], weights=weights
                )
            priority = (utility / cost if cost else np.inf, utility)
            if priority > key:
                chosen, key, is_group = missing, priority, True
        if chosen is None:
            break
        take_group(state, chosen, selected)
        actions += int(is_group)
    return state, actions


def force_then_fill(pool, model, groups):
    state = Selection(pool)
    selected = np.zeros(len(pool.paths), dtype=bool)
    forced = 0
    for group in groups:
        if np.count_nonzero(~state.used[group.pixels]) <= state.remaining:
            take_group(state, group.members, selected)
            forced += 1
    return fill_singles(state, model), forced


def summarize_selection(pool, state, row, output, crossings, actions):
    summary = emit(pool, state, row, output)
    selected = set(state.order)
    summary.update(
        source_crossings=len(crossings),
        positive_support_groups=sum(g.support > 0 for g in crossings),
        joint_actions=actions,
        complete_source_crossings=sum(
            selected.issuperset(g.members) for g in crossings
        ),
    )
    return summary


def experiment(args):
    rows = roster("development")
    pinned = {r["id"]: r for r in read_csv(EVIDENCE / "development/selection.csv")}
    pools, crossings, training = {}, {}, {}
    for row in rows:
        name = row["id"]
        p = read_pool(args.staged, args.base, args.raster, row)
        pools[name] = p
        crossings[name] = read_crossings(args.staged, args.base, name, p)
        artist = mask(args.base / (name + ".artist.drawing.png")).ravel(order="F")
        refs = json.loads((args.staged / (name + ".json")).read_text())
        training[name] = training_rows(p, artist, label_junctions(p, refs), name)
    folders = {k: args.output / k for k in ("baseline", "grouped", "diagnostic")}
    results = {k: [] for k in folders}
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=True)
    models = {}
    for f in range(1, 6):
        train = [r for r in rows if int(r["fold"]) != f]
        held = [r for r in rows if int(r["fold"]) == f]
        if len(train) != 48 or len(held) != 12:
            raise ValueError("Unexpected fold sizes")
        models[f] = fit([training[r["id"]] for r in train])
        for row in held:
            name = row["id"]
            p = pools[name]
            groups = crossings[name]
            old = predict_paths(p, models[f])
            result = summarize_selection(p, old, row, folders["baseline"], groups, 0)
            if result["segments_sha256"] != pinned[name]["segments_sha256"]:
                raise ValueError(f"Frozen single-path replay differs for {name}")
            results["baseline"].append(result)
            new, actions = select_groups(p, models[f], groups)
            results["grouped"].append(
                summarize_selection(p, new, row, folders["grouped"], groups, actions)
            )
        print(
            f"Fold {f}: 12 baseline replays and grouped predictions complete",
            flush=True,
        )
    # The diagnostic uses held-photo labels only after every primary prediction exists.
    from score_crossing_groups import four_arm_pairs

    for row in rows:
        name = row["id"]
        p = pools[name]
        groups = crossings[name]
        refs = json.loads((args.staged / (name + ".json")).read_text())
        drawing = json.loads((args.base / (name + ".mlsd.contours.json")).read_text())
        pairs = four_arm_pairs(refs, drawing)
        matched = {j for _, j in pairs}
        requested = [g for g in groups if g.source == "mlsd" and g.index in matched]
        if len(requested) != len(matched):
            raise ValueError("Missing diagnostic members")
        state, forced = force_then_fill(p, models[int(row["fold"])], requested)
        result = summarize_selection(
            p, state, row, folders["diagnostic"], groups, forced
        )
        result["reference_groups_requested"] = len(requested)
        results["diagnostic"].append(result)
    ledgers = [dict(method=key, **result) for key in folders for result in results[key]]
    write_csv(args.output / "selections.csv", ledgers, list(ledgers[-1]))
    for key, folder in folders.items():
        write_csv(folder / "selection.csv", results[key], list(results[key][0]))
        write_csv(
            folder / "r-plan.csv",
            [dict(id=r["id"], method="curve_rank") for r in rows],
            ["id", "method"],
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("staged", "base", "raster", "output"):
        parser.add_argument(argument, type=Path)
    with threadpool_limits(limits=2):
        experiment(parser.parse_args())
