"""Propose additional X junctions from four locally supported arms in a photograph."""
import argparse
from pathlib import Path

import cv2
import numpy as np
from PIL import Image

from score import read_csv, write_csv


def support_intervals(segments, lengths, directions, seed):
    start = segments[seed, 0]
    axis = directions[seed]
    normal = np.array([-axis[1], axis[0]])
    delta = segments - start
    aligned = np.abs(directions @ axis) >= np.cos(np.deg2rad(5))
    close = np.abs(delta @ normal).max(axis=1) <= 2
    intervals = np.sort(delta[aligned & close] @ axis, axis=1)
    # Keep the consecutive fragments around the seed. Intervals retain the
    # unsupported gaps so they cannot vote as observed ink near a crossing.
    merged = []
    for lo, hi in sorted(map(tuple, intervals)):
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    ranges = np.array(merged)
    breaks = np.flatnonzero(ranges[1:, 0] - ranges[:-1, 1] > 10) + 1
    for group in np.split(ranges, breaks):
        if group[0, 0] <= lengths[seed] / 2 <= group[-1, 1]:
            return group
    raise ValueError("The seed line must support its own midpoint")


def arm_strength(intervals, at):
    def coverage(lo, hi):
        return np.maximum(
            0, np.minimum(intervals[:, 1], hi) - np.maximum(intervals[:, 0], lo)
        ).sum()

    left = coverage(at - 10, at - 1)
    right = coverage(at + 1, at + 10)
    if min(left, right) < 5.5:
        return 0.0
    return (
        min(left, right)
        + min(coverage(at - 24, at - 1), coverage(at + 1, at + 24)) / 23
    )


def crossings(rgb):
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    gray = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8)).apply(gray)
    zoom = cv2.resize(gray, None, fx=2, fy=2, interpolation=cv2.INTER_LINEAR)
    lines = cv2.createLineSegmentDetector(cv2.LSD_REFINE_STD).detect(zoom)[0]
    if lines is None:
        return []
    # Invert the pixel-center resize mapping and deliver one-based coordinates.
    segments = (lines.reshape(-1, 2, 2).astype(float) + 0.5) / 2 + 0.5
    vectors = segments[:, 1] - segments[:, 0]
    lengths = np.linalg.norm(vectors, axis=1)
    keep = lengths >= 4
    segments, vectors, lengths = segments[keep], vectors[keep], lengths[keep]
    dirs = vectors / lengths[:, None]
    seeds = np.flatnonzero(lengths >= 8)
    intervals = [support_intervals(segments, lengths, dirs, s) for s in seeds]
    axes, starts = dirs[seeds], segments[seeds, 0]
    scores = []
    for i, axis in enumerate(axes):
        others = np.arange(i + 1, len(axes))
        denom = axis[0] * axes[others, 1] - axis[1] * axes[others, 0]
        eligible = np.abs(denom) >= 0.5  # At least 30 degrees between the lines.
        others, denom = others[eligible], denom[eligible]
        delta = starts[others] - starts[i]
        u = (delta[:, 0] * axes[others, 1] - delta[:, 1] * axes[others, 0]) / denom
        v = (delta[:, 0] * axis[1] - delta[:, 1] * axis[0]) / denom
        group = intervals[i]
        local = (u >= group[0, 0] + 6) & (u <= group[-1, 1] - 6)
        for j, a, b in zip(others[local], u[local], v[local]):
            other = intervals[j]
            if b < other[0, 0] + 6 or b > other[-1, 1] - 6:
                continue
            p = starts[i] + a * axis
            if np.any(p < 5) or np.any(p > [gray.shape[1] - 4, gray.shape[0] - 4]):
                continue
            strength = min(arm_strength(group, a), arm_strength(other, b))
            if strength:
                scores.append((strength, p[0], p[1]))
    return sorted(scores, reverse=True)


def predict(folder):
    rows = []
    for photo in read_csv(folder / "panel.csv"):
        image_id = photo["id"]
        with Image.open(folder / f"{image_id}.small.png") as image:
            points = crossings(np.array(image.convert("RGB")))
        rows.extend(dict(id=image_id, strength=s, x=x, y=y) for s, x, y in points)
    write_csv(folder / "fragment-candidates.csv", rows, ["id", "strength", "x", "y"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    predict(parser.parse_args().directory)
