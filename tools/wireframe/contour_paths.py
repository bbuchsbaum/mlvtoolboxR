"""Experimental, deterministic contour checks for photographic junctions.

The module deliberately separates *ranking* from the published context scorer.
``resolve_contour_guided(..., mode="baseline")`` delegates exactly to
``topology.resolve``.  ``mode="both"`` uses the same calibrated ordering,
but removes a T only when its proposed stem has a measured break and lets the
next ranked candidate fill the vacant location.  It also treats close X sites
as duplicates only if their two observed line tracks agree.

Coordinates use the wireframe convention: one-based pixel centres.  Returned
records are JSON-ready and contain enough line/ray geometry to redraw the
decision, including the sample-level gap mask used for a T veto.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import atan2, degrees

import cv2
import numpy as np

from topology import anchor, resolve


@dataclass(frozen=True)
class _Track:
    number: int
    members: tuple[int, ...]
    direction: np.ndarray
    offset: float
    extent: tuple[float, float]


def _as_gray(image):
    image = np.asarray(image)
    if image.ndim == 3:
        image = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
    if image.ndim != 2:
        raise ValueError("image must be a grayscale or RGB raster")
    return image.astype(np.float32)


def _line_array(lines):
    lines = np.asarray(lines, dtype=float)
    if lines.size == 0:
        return np.empty((0, 4), float)
    lines = lines.reshape(-1, lines.shape[-1])
    if lines.shape[1] < 4:
        raise ValueError("lines must have x1, y1, x2, y2 columns")
    lines = lines[:, :4]
    delta = lines[:, 2:] - lines[:, :2]
    return lines[np.linalg.norm(delta, axis=1) > 1e-6]


def _union_find(n):
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def join(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[max(a, b)] = min(a, b)

    return find, join


def build_tracks(lines, *, angle_degrees=10.0, offset_pixels=1.5, join_gap=4.0):
    """Group overlapping, nearly collinear observed segments into stable tracks.

    Track IDs are deterministic and refer to the input row numbers (``segment_id``).
    The modest tolerances describe duplicate detector fragments, not junction
    proximity; intersecting tracks are never merged.
    """
    lines = _line_array(lines)
    n = len(lines)
    if not n:
        return [], lines
    delta = lines[:, 2:] - lines[:, :2]
    length = np.linalg.norm(delta, axis=1)
    unit = delta / length[:, None]
    # Give a collinear line one canonical direction, making parallel fragments
    # compare identically even if their detector orientation differs.
    flip = (unit[:, 0] < 0) | ((np.abs(unit[:, 0]) < 1e-12) & (unit[:, 1] < 0))
    unit[flip] *= -1
    normal = np.column_stack((-unit[:, 1], unit[:, 0]))
    # A global origin is easier: projection on the unit direction is translation
    # invariant, so segment extents can be compared after the same direction test.
    find, join = _union_find(n)
    cosine = np.cos(np.deg2rad(angle_degrees))
    for i in range(n):
        for j in range(i):
            # Compare the actual perpendicular separation of the two observed
            # fragments.  This is invariant to translating the whole drawing;
            # it does not compare arbitrary intercepts from different origins.
            separation = abs(
                ((lines[i, :2] + lines[i, 2:] - lines[j, :2] - lines[j, 2:]) / 2)
                @ normal[i]
            )
            if abs(unit[i] @ unit[j]) < cosine or separation > offset_pixels:
                continue
            # Project both fragments onto i's frame before comparing their
            # extents.  Separate local directions are not a common coordinate
            # system, even when the two fragments are almost parallel.
            i_extent = np.sort([lines[i, :2] @ unit[i], lines[i, 2:] @ unit[i]])
            j_extent = np.sort([lines[j, :2] @ unit[i], lines[j, 2:] @ unit[i]])
            if (
                max(i_extent[0], j_extent[0])
                <= min(i_extent[1], j_extent[1]) + join_gap
            ):
                join(i, j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    tracks = []
    for number, members in enumerate(sorted(groups.values(), key=lambda x: min(x))):
        ref = members[0]
        direction = unit[ref].copy()
        # Orient all fragments consistently with the representative.
        projections = []
        for i in members:
            projections.extend([(lines[i, :2] @ direction), (lines[i, 2:] @ direction)])
        tracks.append(
            _Track(
                number,
                tuple(members),
                direction,
                float(lines[ref, :2] @ np.array([-direction[1], direction[0]])),
                (float(min(projections)), float(max(projections))),
            )
        )
    return tracks, lines


def _track_geometry(track, lines):
    segments = [
        {"segment_id": int(i), "polyline": lines[i].reshape(2, 2).tolist()}
        for i in track.members
    ]
    return {
        "track_id": int(track.number),
        "segment_ids": [int(i) for i in track.members],
        "segments": segments,
        "direction": track.direction.tolist(),
        "angle_degrees": float(degrees(atan2(track.direction[1], track.direction[0]))),
    }


def incident_tracks(point, tracks, lines, *, distance=2.5, endpoint_gap=2.5):
    """Return tracks with an observed local arm at ``point`` and their rays."""
    point = np.asarray(point, float)
    result = []
    for track in tracks:
        normal = np.array([-track.direction[1], track.direction[0]])
        perpendicular = abs(point @ normal - track.offset)
        coordinate = point @ track.direction
        if perpendicular > distance or not (
            track.extent[0] - endpoint_gap
            <= coordinate
            <= track.extent[1] + endpoint_gap
        ):
            continue
        rays = []
        for sign, end in ((-1.0, track.extent[0]), (1.0, track.extent[1])):
            reach = sign * (end - coordinate)
            if reach > distance:
                endpoint = point + sign * reach * track.direction
                rays.append(
                    {
                        "direction": (sign * track.direction).tolist(),
                        "endpoint": endpoint.tolist(),
                        "length": float(reach),
                    }
                )
        if rays:
            item = _track_geometry(track, lines)
            item["perpendicular_distance"] = float(perpendicular)
            item["rays"] = rays
            result.append(item)
    return result


def _image_context(image):
    gray = _as_gray(image)
    gx = cv2.Scharr(gray, cv2.CV_32F, 1, 0) / 16
    gy = cv2.Scharr(gray, cv2.CV_32F, 0, 1) / 16
    scale = np.percentile(np.hypot(gx, gy), 95) if gray.size else 0
    return gray, float(max(4.0, 0.10 * scale))


def _sample_stem(
    gray,
    point,
    ray,
    track,
    *,
    contrast_threshold,
    max_gap_samples=1,
    near_bar_pixels=14,
):
    """Trace one ray using transverse contrast, returning drawable veto evidence."""
    length = float(ray["length"])
    direction = np.asarray(ray["direction"], float)
    normal = np.array([-direction[1], direction[0]])
    # Skip the ambiguous centre but include the endpoint; samples retain the
    # one-based convention expected by the benchmark and drawing renderer.
    distances = np.arange(1.0, max(1.0, length) + 0.001, 1.0)
    xy = np.asarray(point, float)[None, :] + distances[:, None] * direction

    def remap(offset):
        # OpenCV image coordinates are zero-based pixel centres.
        return cv2.remap(
            gray,
            (xy + offset * normal - 1).astype(np.float32).reshape(-1, 1, 2),
            None,
            cv2.INTER_LINEAR,
            borderMode=cv2.BORDER_REPLICATE,
        ).reshape(-1)

    centre = remap(0)
    contrasts = []
    for width in (1.2, 2.5, 4.0):
        plus, minus = remap(width), remap(-width)
        # A HAWP line can lie on an edge (bilateral contrast) or down the
        # centre of a thick stroke (centre-to-side contrast).
        contrasts.append(
            np.maximum.reduce(
                [np.abs(plus - minus), np.abs(centre - plus), np.abs(centre - minus)]
            )
        )
    contrast = np.maximum.reduce(contrasts)
    # A contrast-derived threshold is image-local and has no fitted/photo-panel
    # parameter.  The floor handles uniformly lit synthetic images.
    threshold = contrast_threshold
    gap_mask = contrast < threshold
    longest = current = 0
    for missing in gap_mask:
        current = current + 1 if missing else 0
        longest = max(longest, current)
    local = gap_mask & (distances > 2) & (distances <= near_bar_pixels)
    # Reject only a witnessed break next to the bar when the same ray resumes
    # farther out.  A dark/noisy/short ray is unknown evidence, not a veto.
    veto = False
    veto_mask = np.zeros(len(gap_mask), dtype=bool)
    for start in np.flatnonzero(local):
        # The first two pixels lie in the bar/stem overlap and cannot identify
        # a missing join.  They are retained as evidence but never form a veto.
        if start and local[start - 1]:
            continue
        end = start
        while end + 1 < len(gap_mask) and gap_mask[end + 1]:
            end += 1
        if end - start + 1 > max_gap_samples and (~gap_mask[end + 1 :]).sum() >= 3:
            veto = True
            veto_mask[start : end + 1] = True
            break
    return {
        "track_id": int(track["track_id"]),
        "segment_id": int(track["segment_ids"][0]),
        "ray_direction": direction.tolist(),
        "polyline": [np.asarray(point, float).tolist(), ray["endpoint"]],
        "sample_xy": xy.tolist(),
        "contrast": contrast.tolist(),
        "contrast_threshold": threshold,
        "gap_mask": gap_mask.astype(bool).tolist(),
        "veto_mask": veto_mask.astype(bool).tolist(),
        "max_gap_samples": int(longest),
        "near_bar_pixels": int(near_bar_pixels),
        "veto": bool(veto),
        "accepted": bool(not veto),
    }


def t_stem_evidence(image, point, tracks, lines, *, max_gap_samples=1, context=None):
    """Test whether a locally supported stem is visibly interrupted near a through-bar."""
    arms = incident_tracks(point, tracks, lines)
    # A T has a through-bar (two observed rays) and a distinct one-sided stem.
    bars = [a for a in arms if len(a["rays"]) >= 2]
    stems = [
        a
        for a in arms
        if len(a["rays"]) == 1
        and any(
            abs(
                a["direction"][0] * bar["direction"][1]
                - a["direction"][1] * bar["direction"][0]
            )
            >= 0.5
            for bar in bars
        )
    ]
    # A genuinely open stem may end outside the contact tolerance. Consider
    # its continuation only when the bar is observed and the stem is unambiguous.
    if bars and not stems:
        nearby = incident_tracks(point, tracks, lines, endpoint_gap=10)
        possible = [
            t
            for t in nearby
            if len(t["rays"]) == 1
            and any(
                abs(
                    t["direction"][0] * b["direction"][1]
                    - t["direction"][1] * b["direction"][0]
                )
                >= 0.5
                for b in bars
            )
        ]
        if len(possible) == 1:
            stems = possible
            arms = arms + possible
    if not bars or not stems:
        return {
            "accepted": True,
            "veto": False,
            "reason": "open_geometry",
            "incident_tracks": arms,
        }
    stem = max(stems, key=lambda a: a["rays"][0]["length"])
    gray, threshold = _image_context(image) if context is None else context
    result = _sample_stem(
        gray,
        point,
        stem["rays"][0],
        stem,
        contrast_threshold=threshold,
        max_gap_samples=max_gap_samples,
    )
    result["reason"] = "contour_gap" if result["veto"] else "continuous"
    result["incident_tracks"] = arms
    return result


def _x_signature(incident):
    """Two closest nonparallel through-tracks, or an empty unproven identity."""
    through = sorted(
        (a for a in incident if len(a["rays"]) >= 2),
        key=lambda a: (a["perpendicular_distance"], a["track_id"]),
    )
    for i, first in enumerate(through):
        u = np.asarray(first["direction"])
        for second in through[i + 1 :]:
            v = np.asarray(second["direction"])
            if abs(u[0] * v[1] - u[1] * v[0]) >= 0.5:
                return tuple(sorted((first["track_id"], second["track_id"])))
    return ()


def _ranked_from_resolve_inputs(xy, prob, base, cutoffs):
    """Mirror topology.resolve's candidate type/ordering before suppression."""
    xy, prob = np.asarray(xy, float), np.asarray(prob, float)
    thresholds = np.asarray(cutoffs, float)
    if len(thresholds) == 2:
        threshold = np.broadcast_to(thresholds, prob.shape)
    else:
        xy, prob, anchored = anchor(xy, prob, base)
        threshold = np.where(anchored[:, None], thresholds[:2], thresholds[2:])
    odds = prob / np.maximum(1 - prob, 1e-8)
    quality = odds * ((1 - threshold) / threshold)
    kind = quality.argmax(1)
    return [
        dict(
            x=float(p[0]),
            y=float(p[1]),
            type=("T", "X")[kind[i]],
            confidence=float(prob[i, kind[i]]),
            quality=float(quality[i, kind[i]]),
        )
        for i, p in enumerate(xy)
        if quality[i, kind[i]] >= 1
    ]


def _prepared_candidates(xy, prob, base, cutoffs, image, lines, max_gap_samples):
    """Compute line groups, image scale, and candidate evidence once per photo."""
    tracks, clean_lines = build_tracks(lines)
    context = _image_context(image)
    candidates = sorted(
        _ranked_from_resolve_inputs(xy, prob, base, cutoffs),
        key=lambda c: -c["quality"],
    )
    for candidate in candidates:
        point = np.array([candidate["x"], candidate["y"]])
        candidate["incident_tracks"] = incident_tracks(point, tracks, clean_lines)
        if candidate["type"] == "T":
            stem = t_stem_evidence(
                image,
                point,
                tracks,
                clean_lines,
                max_gap_samples=max_gap_samples,
                context=context,
            )
            candidate["stem_path"] = stem
            candidate["contour_score"] = 1.0 - float(not stem["accepted"])
        else:
            candidate["track_signature"] = _x_signature(candidate["incident_tracks"])
    return candidates, tracks, clean_lines, context


def resolve_contour_guided(
    xy,
    prob,
    base,
    cutoffs,
    image,
    lines,
    *,
    mode="both",
    radius=4,
    max_gap_samples=1,
    return_audit=False,
    _prepared=None,
):
    """Select one factorial contour mode from existing scorer inputs.

    The baseline route is intentionally a literal call to ``topology.resolve``.
    In contour mode, only a broken T is vetoed.  Processing continues down the
    same calibrated ordering, which makes replenishment explicit and enables
    equal-recall ranking comparisons via the retained ``quality`` values.
    """
    if mode == "baseline":
        selected = resolve(xy, prob, base, cutoffs)
        if return_audit:
            candidates, tracks, clean_lines, context = (
                _prepared
                or _prepared_candidates(
                    xy, prob, base, cutoffs, image, lines, max_gap_samples
                )
            )
            selected_geometry = _event_geometry(
                selected, tracks, clean_lines, image, context, max_gap_samples
            )
            selected_xy = {(r["x"], r["y"], r["type"]) for r in selected}
            audit = [
                dict(
                    r,
                    selection=(
                        "selected_exact_baseline"
                        if (r["x"], r["y"], r["type"]) in selected_xy
                        else "not_selected_exact_baseline"
                    ),
                )
                for r in candidates
            ]
            return {"selected": selected_geometry, "candidates": audit}
        return selected
    if mode not in ("stem", "tracks", "both"):
        raise ValueError("mode must be 'baseline', 'stem', 'tracks', or 'both'")
    use_stem, use_tracks = mode in ("stem", "both"), mode in ("tracks", "both")
    candidates, tracks, clean_lines, context = _prepared or _prepared_candidates(
        xy, prob, base, cutoffs, image, lines, max_gap_samples
    )
    kept = []
    audit = [dict(candidate) for candidate in candidates]
    for candidate in audit:
        if candidate["type"] == "T":
            stem = candidate["stem_path"]
            if use_stem and stem["veto"]:
                candidate["selection"] = "rejected_contour_gap"
                continue
        close = [
            r
            for r in kept
            if np.hypot(r["x"] - candidate["x"], r["y"] - candidate["y"]) <= radius
        ]
        if (
            use_tracks
            and candidate["type"] == "X"
            and len(candidate["track_signature"]) == 2
            and all(
                r["type"] == "X"
                and len(r.get("track_signature", ())) == 2
                and r["track_signature"] != candidate["track_signature"]
                for r in close
            )
        ):
            candidate["selection"] = "selected_distinct_tracks"
            kept.append(candidate)
        elif not close:
            candidate["selection"] = "selected"
            kept.append(candidate)
        else:
            candidate["selection"] = "rejected_radius"
    return {"selected": kept, "candidates": audit} if return_audit else kept


def _event_geometry(events, tracks, clean_lines, image, context, max_gap_samples):
    out = []
    for event in events:
        record = dict(event)
        point = [record["x"], record["y"]]
        record["incident_tracks"] = incident_tracks(point, tracks, clean_lines)
        if record["type"] == "X":
            record["track_signature"] = _x_signature(record["incident_tracks"])
        elif image is not None:
            record["stem_path"] = t_stem_evidence(
                image,
                point,
                tracks,
                clean_lines,
                max_gap_samples=max_gap_samples,
                context=context,
            )
        out.append(record)
    return out


def resolve_all_modes(
    xy, prob, base, cutoffs, image, lines, *, radius=4, max_gap_samples=1
):
    """Return baseline/stem/tracks/both selections plus common redraw geometry."""
    prepared = _prepared_candidates(
        xy, prob, base, cutoffs, image, lines, max_gap_samples
    )
    audits = {
        name: resolve_contour_guided(
            xy,
            prob,
            base,
            cutoffs,
            image,
            lines,
            mode=name,
            radius=radius,
            max_gap_samples=max_gap_samples,
            return_audit=True,
            _prepared=prepared,
        )
        for name in ("baseline", "stem", "tracks", "both")
    }
    return {name: audit["selected"] for name, audit in audits.items()} | {
        "audit": audits
    }
