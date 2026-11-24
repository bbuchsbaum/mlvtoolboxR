"""Score junction geometry and raster agreement on paired artist drawings."""

import cv2
import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.spatial import cKDTree

from score import positions


def draw_paths(shape, paths):
    """Rasterize one-based pixel-center paths at subpixel precision."""
    mask = np.zeros(shape, np.uint8)
    lines = [
        np.rint((np.asarray(p) - 1) * 256).astype(np.int32)
        for p in paths
        if len(p) >= 2
    ]
    if lines:
        cv2.polylines(mask, lines, False, 1, 1, cv2.LINE_8, shift=8)
    return mask.astype(bool)


def covered(source, target):
    if not target.any():
        return 0
    return int(np.count_nonzero(source & (distance_transform_edt(~target) <= 2)))


def reference_sites(data, kind):
    rows = [j for j in data["junctions"] if j["type"] == kind]
    xy = np.array([j["xy"] for j in rows], dtype=float).reshape(-1, 2)
    directions = []
    for row in rows:
        unit = np.asarray(row["incident_directions"]).reshape(-1, 2)
        directions.append(np.arctan2(unit[:, 1], unit[:, 0]))
    return xy, directions


def junction_counts(truth, predicted, reference_angles=None, predicted_angles=None):
    """Maximum one-to-one matches. Optionally require every arm within 25 degrees."""
    if not len(truth) or not len(predicted):
        return 0
    groups = cKDTree(predicted).query_ball_point(truth, 4)
    rows, columns = [], []
    for i, group in enumerate(groups):
        for j in group:
            if reference_angles is not None:
                a, b = reference_angles[i], predicted_angles[j]
                if len(a) != len(b) or len(a) not in (3, 4):
                    continue
                # Sorted circular orders differ only by a cyclic shift. Use a
                # bijection so two observed arms cannot explain the same ray.
                ordered_a, ordered_b = np.sort(a), np.sort(b)
                errors = [
                    np.abs(
                        (ordered_a - np.roll(ordered_b, shift) + np.pi) % (2 * np.pi)
                        - np.pi
                    ).max()
                    for shift in range(len(a))
                ]
                if min(errors) > np.deg2rad(25):
                    continue
            rows.append(i)
            columns.append(j)
    graph = csr_matrix(
        (np.ones(len(rows), dtype=bool), (rows, columns)),
        shape=(len(truth), len(predicted)),
    )
    return int(
        np.count_nonzero(maximum_bipartite_matching(graph, perm_type="column") >= 0)
    )


def point_scores(events, data, baseline):
    row = {}
    for kind in ("T", "X"):
        truth, _ = reference_sites(data, kind)
        ranked = sorted(
            (r for r in events if r["type"] == kind), key=lambda r: -r["confidence"]
        )
        budget = sum(r["type"] == kind for r in baseline)
        row.update(
            {
                f"{kind}_ref": len(truth),
                f"{kind}_pred": len(ranked),
                f"{kind}_hit": junction_counts(truth, positions(ranked)),
                f"{kind}_at_context_count": junction_counts(
                    truth, positions(ranked[:budget])
                ),
                f"{kind}_capped_count": min(budget, len(ranked)),
            }
        )
    return row
