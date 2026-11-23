"""Add pretrained image activations and neighboring line evidence to a candidate."""
import numpy as np
from scipy.spatial import cKDTree

from predict_hawp import sample_context


# Fixed subset of the 140 multiscale descriptors from topology.model_features.
# Includes arm support from the full, edge and cropped line proposals, and pixel gradients.
SCENE_COLUMNS = [3, 4, 5, 11, 18, 25, 32, 40, 41, 48, 77, 78, 85, 114, 115, 121, 122]


def neighborhood(xy, geometry, size):
    values = geometry[:, SCENE_COLUMNS]
    if not len(xy):
        return np.empty((0, 5 * len(SCENE_COLUMNS) + 4), dtype=np.float16)
    tree = cKDTree(xy)
    summary = []
    for radius in (12, 48):
        groups = tree.query_ball_point(xy, radius)
        summary.extend(
            [
                np.stack([values[g].mean(0) for g in groups]),
                np.stack([values[g].max(0) for g in groups]),
                np.log1p([len(g) for g in groups])[:, None],
            ]
        )
    summary.extend(
        [np.broadcast_to(values.mean(0), values.shape), xy / np.asarray(size)]
    )
    return np.column_stack(summary).astype(np.float32).astype(np.float16)


def descriptors(geometry, maps, xy, size):
    # Fitting used float32 geometry and float16 auxiliary measurements.
    return np.column_stack(
        [geometry, sample_context(maps, xy, size), neighborhood(xy, geometry, size)]
    ).astype(np.float32)
