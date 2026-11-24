"""Fixed, photo-only observations for already extracted contour curves.

The module deliberately knows nothing about annotations, fitted scorers, or
selection.  Pixel indices are zero-based NumPy/Ravel-F-compatible indices
(``row + height * column``); segment coordinates remain the one-based pixel
centres used by the contour JSON files.
"""

from __future__ import annotations

import numpy as np


COLUMN_NAMES = (
    "source_mlsd",
    "pixel_count",
    "arc_length",
    "straightness",
    "edge_contrast_mean",
    "edge_contrast_p90",
    "normal_gradient_mean",
    "tangent_gradient_mean",
    "other_source_near_fraction",
    "endpoint_near_other_curve",
)
# More descriptive public spelling for callers that assemble model matrices.
OBSERVATION_COLUMNS = COLUMN_NAMES


def _gray(rgb):
    rgb = np.asarray(rgb)
    if rgb.ndim != 3 or rgb.shape[2] < 3:
        raise ValueError("rgb must have shape (height, width, at least 3)")
    return np.asarray(rgb[..., :3], dtype=float) @ np.array([0.299, 0.587, 0.114])


def _as_segments(path):
    path = np.asarray(path, dtype=float)
    if path.size == 0:
        return np.empty((0, 4), dtype=float)
    path = path.reshape(-1, 4)
    if not np.isfinite(path).all():
        raise ValueError("paths must contain finite x1, y1, x2, y2 values")
    return path


def _curve_direction(segments):
    if not len(segments):
        return np.array([1.0, 0.0])
    points = np.vstack((segments[0, :2], segments[:, 2:]))
    delta = points[-1] - points[0]
    norm = np.linalg.norm(delta)
    if norm > 1e-9:
        return delta / norm
    delta = segments[:, 2:] - segments[:, :2]
    lengths = np.linalg.norm(delta, axis=1)
    if lengths.sum() > 1e-9:
        direction = (delta * lengths[:, None]).sum(axis=0)
        norm = np.linalg.norm(direction)
        if norm > 1e-9:
            return direction / norm
    return np.array([1.0, 0.0])


def _bilinear(image, x, y):
    """Sample zero-based x/y locations, clamping only at the photo border."""
    height, width = image.shape
    x = np.clip(np.asarray(x, float), 0, width - 1)
    y = np.clip(np.asarray(y, float), 0, height - 1)
    x0, y0 = np.floor(x).astype(int), np.floor(y).astype(int)
    x1, y1 = np.minimum(x0 + 1, width - 1), np.minimum(y0 + 1, height - 1)
    wx, wy = x - x0, y - y0
    return (
        (1 - wx) * (1 - wy) * image[y0, x0]
        + wx * (1 - wy) * image[y0, x1]
        + (1 - wx) * wy * image[y1, x0]
        + wx * wy * image[y1, x1]
    )


def _min_squared_distance(a, b, chunk=2048):
    if not len(a) or not len(b):
        return np.full(len(a), np.inf)
    result = np.full(len(a), np.inf)
    for start in range(0, len(a), chunk):
        part = a[start : start + chunk]
        result[start : start + len(part)] = (
            ((part[:, None] - b[None, :]) ** 2).sum(2).min(1)
        )
    return result


def observations(rgb, pixels, paths, sources):
    """Return one fixed feature row per existing FDoG/M-LSD contour.

    ``pixels[i]`` supplies that contour's zero-based Fortran-order pixel
    indices. ``paths[i]`` supplies its one-based ``(x1, y1, x2, y2)`` segments
    and ``sources[i]`` is ``"fdog"`` or ``"mlsd"``.  No fitted state or labels
    are accepted, so all columns are available before any split or selection.
    """
    gray = _gray(rgb)
    height, width = gray.shape
    n = len(paths)
    if len(pixels) != n or len(sources) != n:
        raise ValueError("pixels, paths, and sources must have the same length")
    if any(source not in {"fdog", "mlsd"} for source in sources):
        raise ValueError("sources must contain only 'fdog' or 'mlsd'")

    pixel_xy = []
    clean_paths = []
    for indices, path in zip(pixels, paths):
        indices = np.asarray(indices, dtype=np.int64).reshape(-1)
        if len(indices) and ((indices < 0).any() or (indices >= height * width).any()):
            raise ValueError("pixel indices are outside the image")
        # R matrices and numpy.ravel(order='F') both use row + height * column.
        pixel_xy.append(
            np.column_stack((indices // height, indices % height)).astype(float)
        )
        clean_paths.append(_as_segments(path))

    gy, gx = np.gradient(gray)
    answer = np.empty((n, len(COLUMN_NAMES)), dtype=float)
    for i, (xy, segments, source) in enumerate(zip(pixel_xy, clean_paths, sources)):
        direction = _curve_direction(segments)
        normal = np.array([-direction[1], direction[0]])
        if len(segments):
            deltas = segments[:, 2:] - segments[:, :2]
            arc = float(np.linalg.norm(deltas, axis=1).sum())
            chord = float(np.linalg.norm(segments[-1, 2:] - segments[0, :2]))
            straightness = chord / arc if arc > 1e-9 else 0.0
        else:
            arc, straightness = 0.0, 0.0
        if len(xy):
            x, y = xy[:, 0], xy[:, 1]
            contrast = np.abs(
                _bilinear(gray, x + normal[0], y + normal[1])
                - _bilinear(gray, x - normal[0], y - normal[1])
            )
            grad_x, grad_y = _bilinear(gx, x, y), _bilinear(gy, x, y)
            normal_grad = np.abs(grad_x * normal[0] + grad_y * normal[1])
            tangent_grad = np.abs(grad_x * direction[0] + grad_y * direction[1])
        else:
            contrast = normal_grad = tangent_grad = np.zeros(1)
        other_source = (
            np.vstack(
                [p for j, p in enumerate(pixel_xy) if j != i and sources[j] != source]
            )
            if any(
                j != i and sources[j] != source and len(p)
                for j, p in enumerate(pixel_xy)
            )
            else np.empty((0, 2))
        )
        other_curves = (
            np.vstack([p for j, p in enumerate(pixel_xy) if j != i and len(p)])
            if any(j != i and len(p) for j, p in enumerate(pixel_xy))
            else np.empty((0, 2))
        )
        if len(segments):
            endpoints = np.vstack((segments[0, :2], segments[-1, 2:]))
            endpoints -= 1.0  # path coordinates are one-based; pixel_xy is zero-based.
        elif len(xy):
            endpoints = xy[[0, -1]]
        else:
            endpoints = np.empty((0, 2))
        answer[i] = (
            float(source == "mlsd"),
            len(xy),
            arc,
            straightness,
            float(np.mean(contrast)),
            float(np.percentile(contrast, 90)),
            float(np.mean(normal_grad)),
            float(np.mean(tangent_grad)),
            (
                float(np.mean(_min_squared_distance(xy, other_source) <= 4.0))
                if len(xy)
                else 0.0
            ),
            (
                float(np.mean(_min_squared_distance(endpoints, other_curves) <= 9.0))
                if len(endpoints)
                else 0.0
            ),
        )
    return answer
