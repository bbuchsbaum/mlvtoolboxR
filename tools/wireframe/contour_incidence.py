"""Offline local-contour correspondence for independently traced drawings.

This deliberately does not compare global contour identifiers: tracing may split
one artist contour at a junction.  Instead, an artist junction's *named*
incident polylines and an R site's reported ``members`` are compared only in
the 3--9 pixel annulus around a same-type candidate.  Both directions of
coverage must hold, and exported arm directions are checked when both sides
provide them.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.spatial import cKDTree


def _xy(row):
    value = row["xy"] if "xy" in row else [row["x"], row["y"]]
    return np.asarray(value, dtype=float)


def _segments(contour):
    if "segments" in contour:
        return np.asarray(contour["segments"], dtype=float).reshape(-1, 4)
    points = np.asarray(contour.get("points", []), dtype=float).reshape(-1, 2)
    return (
        np.column_stack((points[:-1], points[1:]))
        if len(points) > 1
        else np.empty((0, 4))
    )


def _contour_map(contours):
    if isinstance(contours, dict):
        contours = contours.get("contours", contours.get("polylines", []))
    return {str(c["id"]): _segments(c) for c in contours}


def _member_ids(site):
    members = site.get("members", site.get("member_ids", []))
    ids = []
    for item in members:
        if isinstance(item, dict):
            value = item.get("contour_id", item.get("contourID", item.get("id")))
        else:
            value = item
        if value is not None:
            ids.append(str(value))
    return tuple(dict.fromkeys(ids))


def _annulus_points(segments, centre, inner=3.0, outer=9.0, step=0.5):
    """Deterministically sample observed geometry, retaining only an annulus."""
    points = []
    for x1, y1, x2, y2 in np.asarray(segments, float).reshape(-1, 4):
        a, b = np.array([x1, y1]), np.array([x2, y2])
        n = max(1, int(math.ceil(np.linalg.norm(b - a) / step)))
        sample = a + np.linspace(0, 1, n + 1)[:, None] * (b - a)
        d = np.linalg.norm(sample - centre, axis=1)
        points.extend(sample[(d >= inner) & (d <= outer)])
    return np.asarray(points, float).reshape(-1, 2)


def _directions(points, centre):
    if not len(points):
        return np.empty((0, 2))
    unit = points - centre
    unit /= np.linalg.norm(unit, axis=1, keepdims=True)
    # One representative per 15 degrees prevents dense sampling from changing
    # the direction test or making a curved/parallel line look like many arms.
    bins = (
        np.rint(np.arctan2(unit[:, 1], unit[:, 0]) * 180 / np.pi / 15).astype(int) % 24
    )
    return np.vstack([unit[np.flatnonzero(bins == b)[0]] for b in np.unique(bins)])


def _covered(source, target, tolerance):
    if not len(source) or not len(target):
        return False
    distance = np.linalg.norm(source[:, None] - target[None, :], axis=2).min(axis=1)
    return bool(np.all(distance <= tolerance))


def _direction_match(left, right, degrees_tolerance=25):
    if not len(left) or not len(right):
        return None
    # Each visible artist ray needs a unique observed ray; the symmetric
    # direction check follows the bidirectional coverage rule.
    cosine = math.cos(math.radians(degrees_tolerance))
    if len(left) != len(right):
        return False
    # Fixed-order augmenting paths provide a real one-to-one arm assignment.
    # Pairwise nearest-neighbour checks could let one observed bar arm explain
    # both artist bar arms.
    edges = [[j for j, b in enumerate(right) if float(a @ b) >= cosine] for a in left]
    assigned = {}

    def augment(i, seen):
        for j in edges[i]:
            if j in seen:
                continue
            seen.add(j)
            if j not in assigned or augment(assigned[j], seen):
                assigned[j] = i
                return True
        return False

    return all(augment(i, set()) for i in range(len(left)))


def _site_directions(site):
    values = site.get("directions", site.get("arms", []))
    if isinstance(values, str):
        values = [float(v) for v in values.split()]
    if isinstance(values, (list, tuple)):
        values = [v for v in values if v is not None]
    values = np.asarray(values, float)
    if values.size == 0:
        return np.empty((0, 2))
    if values.ndim == 1:
        # R exports degrees with mathematical y up; convert to image y down.
        radians = np.deg2rad(-values)
        return np.column_stack((np.cos(radians), np.sin(radians)))
    return values.reshape(-1, 2)


def pair_evidence(
    artist_site,
    artist_contours,
    observed_site,
    observed_contours,
    *,
    annulus=(3, 9),
    tolerance=2,
):
    """Explain whether one same-type spatial candidate has incident agreement."""
    if artist_site["type"] != observed_site["type"]:
        return {"eligible": False, "reason": "type_mismatch"}
    centre_a, centre_o = _xy(artist_site), _xy(observed_site)
    distance = float(np.linalg.norm(centre_a - centre_o))
    if distance > 4:
        return {"eligible": False, "reason": "outside_4px", "distance": distance}
    a_ids, o_ids = tuple(
        map(str, artist_site.get("incident_polyline_ids", []))
    ), _member_ids(observed_site)
    missing_a = [i for i in a_ids if i not in artist_contours]
    missing_o = [i for i in o_ids if i not in observed_contours]
    if not a_ids or missing_a:
        return {
            "eligible": False,
            "reason": "artist_members_missing",
            "missing": missing_a,
        }
    if not o_ids or missing_o:
        return {
            "eligible": False,
            "reason": "observed_members_missing",
            "missing": missing_o,
        }
    a_segments = np.vstack([artist_contours[i] for i in a_ids])
    o_segments = np.vstack([observed_contours[i] for i in o_ids])
    a_points = _annulus_points(a_segments, centre_a, *annulus)
    o_points = _annulus_points(o_segments, centre_o, *annulus)
    if not len(a_points) or not len(o_points):
        return {
            "eligible": False,
            "reason": "insufficient_annulus_geometry",
            "artist_points": len(a_points),
            "observed_points": len(o_points),
        }
    forward, reverse = _covered(a_points, o_points, tolerance), _covered(
        o_points, a_points, tolerance
    )
    geometry = forward and reverse
    if "incident_directions" in artist_site and "arms" in observed_site:
        a_dirs = _site_directions({"directions": artist_site["incident_directions"]})
        o_dirs = _site_directions(observed_site)
        directions = _direction_match(a_dirs, o_dirs)
        direction_status = (
            "checked" if directions is not None else "unchecked_empty_directions"
        )
    else:
        directions = None
        direction_status = (
            "unchecked_artist_directions_absent"
            if "incident_directions" not in artist_site
            else "unchecked_observed_arms_absent"
        )
    return {
        "eligible": True,
        "distance": distance,
        "artist_member_ids": a_ids,
        "observed_member_ids": o_ids,
        "artist_annulus_points": len(a_points),
        "observed_annulus_points": len(o_points),
        "artist_covered": forward,
        "observed_covered": reverse,
        "direction_agrees": directions,
        "direction_status": direction_status,
        "accepted": bool(geometry and directions is not False),
        "reason": (
            "accepted"
            if geometry and directions is not False
            else "direction_mismatch" if geometry else "annulus_coverage_mismatch"
        ),
    }


def incident_correspondence(
    artist, observed_sites, observed_contours, *, annulus=(3, 9), tolerance=2
):
    """Maximum-cardinality one-to-one local correspondence, split by T and X.

    ``artist`` is the staged Toronto geometry dictionary. ``observed_sites``
    accepts the R ``.contours.json`` object or its ``sites`` list, and
    ``observed_contours`` accepts that object's ``contours`` list.  The result
    contains stable counts plus every accepted/rejected local candidate.
    """
    artist_contours = _contour_map(artist["polylines"])
    if isinstance(observed_sites, dict):
        observed_sites = observed_sites.get("sites", [])
    observed_map = _contour_map(observed_contours)
    result = {
        "filters": {
            "spatial_radius": 4,
            "annulus": list(annulus),
            "coverage_tolerance": tolerance,
            "same_type": True,
            "direction_when_exported": True,
        },
        "matches": [],
    }
    for kind in ("T", "X"):
        left = [s for s in artist["junctions"] if s["type"] == kind]
        right = [s for s in observed_sites if s["type"] == kind]
        edges, evidence = [[] for _ in left], []
        tree = cKDTree(np.asarray([_xy(o) for o in right])) if right else None
        for i, a in enumerate(left):
            # This is the only potentially expensive geometry comparison:
            # query the <=4px spatial candidates, never all T/X pairs.
            nearby = tree.query_ball_point(_xy(a), 4) if tree is not None else []
            for j in sorted(nearby):
                o = right[j]
                item = pair_evidence(
                    a,
                    artist_contours,
                    o,
                    observed_map,
                    annulus=annulus,
                    tolerance=tolerance,
                )
                item.update(type=kind, artist_index=i, observed_index=j)
                evidence.append(item)
                if item.get("accepted"):
                    edges[i].append(j)
        # Fixed-order Kuhn search: maximum cardinality, deterministic mapping.
        assigned = {}

        def augment(i, seen):
            for j in edges[i]:
                if j in seen:
                    continue
                seen.add(j)
                if j not in assigned or augment(assigned[j], seen):
                    assigned[j] = i
                    return True
            return False

        for i in range(len(left)):
            augment(i, set())
        pairs = sorted((i, j) for j, i in assigned.items())
        result[kind] = {
            "artist": len(left),
            "observed": len(right),
            "eligible_pairs": sum(x.get("eligible", False) for x in evidence),
            "geometry_pairs": sum(x.get("accepted", False) for x in evidence),
            "matched": len(pairs),
            "pairs": pairs,
            "evidence": evidence,
        }
        result["matches"].extend(
            dict(type=kind, artist_index=i, observed_index=j) for i, j in pairs
        )
    return result
