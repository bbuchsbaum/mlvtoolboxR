"""Render locally connected contour segments without stamping T or X glyphs."""

import numpy as np

from score import cross


def nearest_member(track, point):
    best = None
    for seg in track["segments"]:
        a, b = np.asarray(seg["polyline"], float)
        u = b - a
        t = np.clip((point - a) @ u / (u @ u), 0, 1)
        d = np.linalg.norm(point - a - t * u)
        if best is None or d < best[0]:
            best = (d, seg)
    return best[1]


def extension_lines(point, pair):
    """Join two observed finite segments only at their nearby line intersection."""
    if len(pair) != 2:
        return []
    lines = [np.asarray(nearest_member(t, point)["polyline"], float) for t in pair]
    a, b = lines[0]
    c, d = lines[1]
    u, v = b - a, d - c
    den = cross(u, v)
    if abs(den) < 0.5 * np.linalg.norm(u) * np.linalg.norm(v):
        return []
    hit = a + u * cross(c - a, v) / den
    if np.linalg.norm(hit - point) > 4:
        return []
    ends = [
        a + np.clip((hit - a) @ u / (u @ u), 0, 1) * u
        for a, b in lines
        for u in [b - a]
    ]
    if any(np.linalg.norm(p - hit) > 4 for p in ends):
        return []
    return [[p, hit] for p in ends if np.linalg.norm(p - hit) > 0.01]


def selected_pair(event):
    tracks = event.get("incident_tracks", [])
    if event["type"] == "X":
        signature = event.get("track_signature") or []
        return [t for t in tracks if t["track_id"] in signature]
    stem_path = event.get("stem_path", {})
    tracks = stem_path.get("incident_tracks", tracks)
    stem = stem_path.get("track_id")
    member = next((t for t in tracks if t["track_id"] == stem), None)
    if member is None:
        return []
    bars = [
        t
        for t in tracks
        if t["track_id"] != stem
        and len(t["rays"]) == 2
        and abs(cross(member["direction"], t["direction"])) >= 0.5
    ]
    if not bars:
        return []
    return [member, min(bars, key=lambda t: t["perpendicular_distance"])]


def result_segments(segments, events, rejected):
    """Cut a measured gap only on the nominated stem, preserving other tracks."""
    cuts = {}
    for row in rejected:
        ev = row["stem_path"]
        if ev.get("reason") != "contour_gap":
            continue
        track = next(
            t for t in ev["incident_tracks"] if t["track_id"] == ev["track_id"]
        )
        for seg in track["segments"]:
            i = seg["segment_id"]
            a, b = np.asarray(seg["polyline"], float)
            direction = (b - a) / np.linalg.norm(b - a)
            samples = np.asarray(ev["sample_xy"])[np.asarray(ev["veto_mask"], bool)]
            for sample in samples:
                proj = (sample - a) @ direction
                if np.linalg.norm(sample - a - proj * direction) <= 2:
                    cuts.setdefault(i, []).append((proj - 0.7, proj + 0.7))
    paths = []
    for i, line in enumerate(segments):
        a, b = line.reshape(2, 2)
        length = np.linalg.norm(b - a)
        if length <= 0:
            continue
        if i not in cuts:
            paths.append([a, b])
            continue
        u = (b - a) / length
        start = 0.0
        for lo, hi in sorted(cuts[i]):
            stop = min(length, max(0, lo))
            if stop > start:
                paths.append([a + start * u, a + stop * u])
            start = max(start, min(length, hi))
        if start < length:
            paths.append([a + start * u, b])
    for row in events:
        paths.extend(
            extension_lines(np.array([row["x"], row["y"]]), selected_pair(row))
        )
    return np.asarray(paths).reshape(-1, 4)
