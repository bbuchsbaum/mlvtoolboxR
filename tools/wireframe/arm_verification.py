"""Measure continuous arms on a fixed set of sites, in working-photo coordinates."""
import cv2
import numpy as np


ANGLES = np.arange(36) * (2 * np.pi / 36)
DIRECTION = np.column_stack((np.cos(ANGLES), np.sin(ANGLES)))
NORMAL = np.column_stack((-np.sin(ANGLES), np.cos(ANGLES)))
RADII = np.arange(0.75, 18.5, 0.9)
OFFSETS = np.array([-1.5, -0.75, 0.0, 0.75, 1.5])


def arm_profiles(gray, points, size):
    """For every ray choose one parallel edge track and estimate its continuous extent."""
    # Shape after projection: site, angle, offset, radius, coordinate.
    base = (
        points[:, None, None, None, :]
        - 0.5
        + DIRECTION[None, :, None, None, :] * RADII[None, None, None, :, None]
    )
    center = base + NORMAL[None, :, None, None, :] * OFFSETS[None, None, :, None, None]
    scale = np.asarray(gray.shape[::-1]) / size

    def sample(offset):
        xy = ((center + NORMAL[None, :, None, None, :] * offset) * scale - 0.5).astype(
            np.float32
        )
        # Each row in the sampling grid is one ray at one transverse offset.
        grid = xy.reshape(-1, len(RADII), 2)
        return cv2.remap(
            gray, grid, None, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE
        ).reshape(xy.shape[:-1])

    contrasts = np.concatenate([(sample(w) - sample(-w)) for w in (0.6, 1.25)], axis=2)
    floor = np.maximum(3.0, 0.7 * np.median(np.abs(contrasts), axis=(1, 2, 3)))
    # Use one polarity over the arm. A single weak sample is tolerated.
    sign = np.where(contrasts[..., :10].mean(axis=-1) >= 0, 1.0, -1.0)
    signed = contrasts * sign[..., None]
    weak = signed < floor[:, None, None, None]
    broken = weak[..., 1:] & weak[..., :-1]
    length = np.where(broken.any(axis=-1), broken.argmax(axis=-1), len(RADII))
    near, middle, far = [
        np.maximum(0, signed[..., part].mean(axis=-1))
        for part in (slice(0, 4), slice(4, 10), slice(10, None))
    ]
    strength = np.sqrt((near + 1) * (middle + 1)) * (0.35 + length / len(RADII))
    best = strength.argmax(axis=2)
    a, b = np.indices(best.shape)
    # The last two quantities distinguish an early break from a distant, unrelated edge.
    return np.stack(
        [
            v[a, b, best]
            for v in (
                strength,
                near,
                middle,
                far,
                length * 0.9,
                (~weak[..., :4]).mean(axis=-1),
            )
        ],
        axis=-1,
    ).astype(np.float32)


def measurements(rgb, xy, size):
    """Describe the best joint three- and four-arm configurations at each site."""
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY).astype(np.float32)
    result = []
    axis = np.arange(18)
    separation = np.abs(axis[:, None] - axis)
    admissible = (separation >= 3) & (separation <= 15)
    for points in (xy[k : k + 16] for k in range(0, len(xy), 16)):
        rays = arm_profiles(gray, points, np.asarray(size))
        energy = rays[:, :, 0]
        through = np.minimum(energy[:, :18], energy[:, 18:])
        stem = np.maximum(energy[:, :18], energy[:, 18:])
        rows = np.arange(len(points))
        modes = []
        for second in (stem, through):
            support = np.where(admissible, through[:, :, None] * second[:, None, :], -1)
            pair = support.reshape(len(points), -1).argmax(1)
            first, next_axis = np.divmod(pair, 18)
            arms = []
            for a in (first, next_axis):
                stronger = np.where(energy[rows, a] >= energy[rows, a + 18], a, a + 18)
                arms.extend([rays[rows, stronger], rays[rows, (stronger + 18) % 36]])
            angle = np.minimum(abs(first - next_axis), 18 - abs(first - next_axis)) * 10
            modes.append(np.column_stack([*arms, angle]))
        result.append(np.column_stack(modes))
    return (
        np.vstack(result).astype(np.float32)
        if result
        else np.empty((0, 50), np.float32)
    )


def baseline_features(events):
    is_x = np.asarray([r["type"] == "X" for r in events])
    log_score = np.log(np.clip([float(r["score"]) for r in events], 1e-8, 1))
    return np.column_stack([is_x, log_score, np.where(is_x, log_score, 0)])


def exchange(is_x, advantage):
    """Trade types only in opposite pairs, up to a fifth of the photo's X count."""
    x = np.flatnonzero(is_x)
    t = np.flatnonzero(~is_x)
    count = min(len(x) // 5, len(t))
    drop = x[np.argsort(advantage[x], kind="stable")[:count]]
    add = t[np.argsort(-advantage[t], kind="stable")[:count]]
    beneficial = advantage[add] > advantage[drop]
    result = is_x.copy()
    result[drop[beneficial]] = False
    result[add[beneficial]] = True
    return result
