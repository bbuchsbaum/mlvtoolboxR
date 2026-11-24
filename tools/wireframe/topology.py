"""Measure local T and X arm support from a predicted vector graph and image edges."""

import cv2
import numpy as np
from scipy.ndimage import maximum_filter1d
from scipy.spatial import cKDTree

BINS = 36


def image_lines(rgb):
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    gray = cv2.createCLAHE(clipLimit=2, tileGridSize=(8, 8)).apply(gray)
    zoom = cv2.resize(gray, None, fx=2, fy=2, interpolation=cv2.INTER_LINEAR)
    lines = cv2.createLineSegmentDetector(cv2.LSD_REFINE_STD).detect(zoom)[0]
    if lines is None:
        return np.empty((0, 5))
    xy = (lines.reshape(-1, 4).astype(float) + 0.5) / 2 + 0.5
    return np.column_stack([xy, np.ones(len(xy))])


def geometry(lines):
    start = lines[:, :2]
    delta = lines[:, 2:4] - start
    length = np.linalg.norm(delta, axis=1)
    keep = (length >= 3) & (lines[:, 4] >= 0.015)
    return start[keep], (delta[keep] / length[keep, None]), length[keep], lines[keep, 4]


def cross(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def intersections(geom):
    start, unit, length, confidence = geom
    a, b = np.triu_indices(len(length), 1)
    den = cross(unit[a], unit[b])
    keep = np.abs(den) >= 0.42
    a, b, den = a[keep], b[keep], den[keep]
    diff = start[b] - start[a]
    u, v = cross(diff, unit[b]) / den, cross(diff, unit[a]) / den
    keep = (u >= -3) & (u <= length[a] + 3) & (v >= -3) & (v <= length[b] + 3)
    return start[a[keep]] + u[keep, None] * unit[a[keep]]


def axis_features(hist):
    # A 10-degree angular tolerance on either side allows slightly bent arms.
    hist = maximum_filter1d(hist, size=3, axis=1, mode="wrap")
    pos, neg = np.split(hist, 2, axis=1)
    weak, strong = np.minimum(pos, neg), np.maximum(pos, neg)
    # Two axes must be separated by at least 30 degrees modulo 180.
    opposite = np.stack([np.roll(weak, shift, axis=1) for shift in range(3, 16)]).max(
        axis=0
    )
    x = np.minimum(weak, opposite)
    t = np.minimum(strong, opposite)
    single = np.minimum(strong * (1 - weak), opposite)
    return np.column_stack(
        [
            strong.max(1),
            weak.max(1),
            t.max(1),
            x.max(1),
            single.max(1),
            (hist > 0.15).mean(1),
            (hist > 0.50).mean(1),
        ]
    )


def evidence(points, geom):
    start, unit, length, confidence = geom
    if len(start) == 0:
        return np.zeros((len(points), 37), dtype=np.float32)
    orient = (
        np.rint(np.arctan2(unit[:, 1], unit[:, 0]) * BINS / (2 * np.pi)).astype(int)
        % BINS
    )
    result = []
    for chunk in np.array_split(points, max(1, (len(points) + 95) // 96)):
        diff = chunk[:, None, :] - start
        along = (diff * unit).sum(axis=2)
        perpendicular = np.abs(cross(diff, unit))
        outside = np.maximum(0, np.maximum(-along, along - length))
        dist = np.hypot(perpendicular, outside)
        features = [
            dist.min(1),
            np.hypot(perpendicular, np.minimum(abs(along), abs(length - along))).min(1),
        ]
        for sigma, arm, power in [
            (1, 5, 1),
            (2.5, 6, 1),
            (2, 12, 0.35),
            (3, 16, 1),
            (1.3, 5, 0.25),
        ]:
            strength = np.exp(-0.5 * (dist / sigma) ** 2) * confidence**power
            forward = strength * np.clip((length - along) / arm, 0, 1)
            back = strength * np.clip(along / arm, 0, 1)
            hist = np.zeros((len(chunk), BINS))
            ix = np.arange(len(chunk))[:, None]
            np.maximum.at(hist, (ix, orient), forward)
            np.maximum.at(hist, (ix, (orient + BINS // 2) % BINS), back)
            features.extend(axis_features(hist).T)
        result.append(np.column_stack(features))
    return np.vstack(result).astype(np.float32)


def features(rgb, wire, alternate=None):
    h, w = rgb.shape[:2]
    learned = geometry(wire)
    edges = geometry(image_lines(rgb))
    # Exact learned vertices precede inferred intersections in each 2-pixel bin.
    vertices = np.unique(wire[wire[:, 4] > 0.02, :4].reshape(-1, 2), axis=0)
    sources = [learned, edges]
    candidates = [vertices, intersections(learned), intersections(edges)]
    if alternate is not None:
        detail = geometry(alternate)
        sources.append(detail)
        candidates += [
            np.unique(alternate[:, :4].reshape(-1, 2), axis=0),
            intersections(detail),
        ]
    proposals = np.vstack(candidates)
    proposals = proposals[
        np.all((proposals >= 3) & (proposals <= [w - 2, h - 2]), axis=1)
    ]
    _, first = np.unique(np.rint(proposals / 2).astype(int), axis=0, return_index=True)
    xy = proposals[np.sort(first)]
    if not len(xy):
        return xy, np.empty((0, 37 * len(sources) + 1), dtype=np.float32)
    nearest = (
        cKDTree(vertices).query(xy)[0] if len(vertices) else np.full(len(xy), 256.0)
    )
    f = np.column_stack([*(evidence(xy, source) for source in sources), nearest])
    return xy, f.astype(np.float32)


def suppress(xy, prob, cutoffs, radius=4):
    """Choose one type per point and suppress competing sites in confidence order."""
    if not len(xy):
        return []
    # Thresholds reflect different positive base rates. Rank by calibrated odds
    # relative to the selected per-class threshold.
    odds = prob / np.maximum(1 - prob, 1e-8)
    quality = odds * ((1 - np.asarray(cutoffs)) / np.asarray(cutoffs))
    kind = quality.argmax(axis=1)
    value = quality.max(axis=1)
    tree = cKDTree(xy)
    blocked = np.zeros(len(xy), dtype=bool)
    kept = []
    for i in np.argsort(-value, kind="stable"):
        if value[i] < 1:
            break
        if blocked[i]:
            continue
        blocked[tree.query_ball_point(xy[i], radius)] = True
        kept.append(
            dict(
                x=float(xy[i, 0]),
                y=float(xy[i, 1]),
                type=("T", "X")[kind[i]],
                confidence=float(prob[i, kind[i]]),
            )
        )
    return kept


def tree_response(features, model):
    """Sum the fitted bias and leaf values of a numeric tree ensemble."""
    total = np.full(len(features), model["bias"], dtype=float)
    for tree in model["trees"]:
        nodes = np.zeros(len(features), dtype=int)
        leaves = np.asarray(tree["leaf"], dtype=bool)
        f = np.asarray(tree["feature"])
        t = np.asarray(tree["threshold"])
        left = np.asarray(tree["left"])
        right = np.asarray(tree["right"])
        active = np.flatnonzero(~leaves[nodes])
        while len(active):
            n = nodes[active]
            nodes[active] = np.where(
                features[active, f[n]] <= t[n], left[n], right[n]
            )
            active = active[~leaves[nodes[active]]]
        total += np.asarray(tree["value"])[nodes]
    return total


def probabilities(features, models):
    from scipy.special import expit

    return np.column_stack([expit(tree_response(features, model)) for model in models])


def anchor(xy, prob, prior):
    """Transfer nearby vector evidence to measured contour sites of either type."""
    base = np.array([[float(r["x"]), float(r["y"])] for r in prior]).reshape(-1, 2)
    if len(base) == 0 or len(xy) == 0:
        return xy, prob, np.zeros(len(xy), dtype=bool)
    tree = cKDTree(xy)
    borrowed = np.zeros((len(base), 2))
    for j, idx in enumerate(tree.query_ball_point(base, 5)):
        if len(idx):
            delta = xy[idx] - base[j]
            weight = np.exp(-np.sum(delta * delta, axis=1) / 18)
            borrowed[j] = (prob[idx] * weight[:, None]).max(axis=0)
    # A same-type vote from the earlier tracer supports retention of a baseline
    # site. It is never by itself a reason to introduce a new location.
    kinds = np.array([r["type"] == "X" for r in prior], int)
    borrowed[np.arange(len(base)), kinds] += 0.055
    borrowed = np.minimum(0.99, borrowed)
    novel = cKDTree(base).query(xy)[0] > 4
    return (
        np.vstack([base, xy[novel]]),
        np.vstack([borrowed, prob[novel]]),
        np.r_[np.ones(len(base), bool), np.zeros(novel.sum(), bool)],
    )


def photographic_arms(rgb, points):
    """Measure gradients normal to 36 candidate arms near and far from the center."""
    if not len(points):
        return np.empty((0, 28), np.float32)
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY).astype(np.float32) / 255
    gray = cv2.GaussianBlur(gray, (3, 3), 0.65)
    dx = cv2.Scharr(gray, cv2.CV_32F, 1, 0) / 16
    dy = cv2.Scharr(gray, cv2.CV_32F, 0, 1) / 16
    angles = np.arange(BINS) * 2 * np.pi / BINS
    directions = np.column_stack([np.cos(angles), np.sin(angles)])
    normals = np.column_stack([-directions[:, 1], directions[:, 0]])
    radii = np.array([1.5, 3, 4.5, 6, 8, 10, 13, 17])
    offset = (
        radii[None, :, None, None] * directions[:, None, None, :]
        + np.arange(-2, 3)[None, None, :, None] * normals[:, None, None, :]
    )
    all_features = []
    for chunk in np.array_split(points, max(1, (len(points) + 47) // 48)):
        grid = (chunk[:, None, None, None, :] - 1 + offset).astype(np.float32)
        shape = grid.shape[:-1]
        grid = grid.reshape(len(chunk), -1, 2)
        gx = cv2.remap(
            dx, grid, None, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE
        ).reshape(shape)
        gy = cv2.remap(
            dy, grid, None, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE
        ).reshape(shape)
        signed = gx * normals[:, 0, None, None] + gy * normals[:, 1, None, None]
        # On each straight ray the same transverse offset must explain all radii.
        # Absolute polarity is immaterial for a light or dark edge.
        groups = []
        for indices in (slice(0, 3), slice(2, 6), slice(4, 8), slice(0, 8)):
            observed = np.abs(signed[:, :, indices, :].mean(axis=2)).max(axis=2)
            groups.append(axis_features(np.tanh(5 * observed)))
        all_features.append(np.column_stack(groups))
    return np.vstack(all_features).astype(np.float32)


def resolve(xy, prob, base, cutoffs):
    """Use different score cutoffs at traced sites and at new locations."""
    if len(cutoffs) == 2:
        return suppress(xy, prob, cutoffs)
    xy, prob, anchored = anchor(xy, prob, base)
    thresholds = np.where(anchored[:, None], [cutoffs[:2]], [cutoffs[2:]])
    return suppress(xy, prob, thresholds)


def model_features(rgb, wire, variant, detail=None):
    if (variant == "multiscale") != (detail is not None):
        raise ValueError("Multiscale scoring requires cropped wire predictions")
    # The benchmark uses decimal round trips of the FP32 line output. Promote
    # before intersections so direct inference uses identical geometry.
    wire = np.asarray(wire, dtype=np.float64)
    if detail is not None:
        detail = np.asarray(detail, dtype=np.float64)
    xy, measured = features(rgb, wire, detail)
    if variant != "arms":
        measured = np.column_stack([measured, photographic_arms(rgb, xy)])
    return xy, measured
