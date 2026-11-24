"""Synthetic controls for the experimental contour-guided selector."""

import unittest

import numpy as np

from contour_paths import (
    build_tracks,
    resolve_all_modes,
    resolve_contour_guided,
    t_stem_evidence,
)
from topology import resolve


def tee_image(broken=False):
    image = np.full((81, 81), 255, np.uint8)
    image[38:43, 12:69] = 0
    image[40:70, 38:43] = 0
    if broken:
        image[52:58, 38:43] = 255
    # One bar and one stem, all coordinates are one-based centres.
    lines = np.array([[13, 41, 69, 41], [41, 41, 41, 70]], float)
    return image, lines


class TestContourPaths(unittest.TestCase):
    def test_intact_stem_has_reconstructable_continuous_path(self):
        image, lines = tee_image()
        tracks, clean = build_tracks(lines)
        evidence = t_stem_evidence(image, [41, 41], tracks, clean)
        self.assertTrue(evidence["accepted"])
        self.assertEqual(evidence["reason"], "continuous")
        self.assertEqual(evidence["polyline"], [[41.0, 41.0], [41.0, 70.0]])
        self.assertEqual(evidence["segment_id"], 1)
        self.assertEqual(len(evidence["sample_xy"]), len(evidence["gap_mask"]))
        # Centre overlap can be low-contrast; it is evidence, not a veto.
        self.assertFalse(any(evidence["veto_mask"]))

    def test_broken_stem_is_vetoed_and_lower_candidate_replenishes(self):
        image, lines = tee_image(broken=True)
        # A second intact T is farther right and lower ranked.
        image[38:43, 50:78] = 0
        image[40:70, 62:67] = 0
        lines = np.vstack([lines, [[51, 41, 78, 41], [64, 41, 64, 70]]])
        xy = np.array([[41, 41], [64, 41]], float)
        # Both are calibrated T candidates; the broken higher-ranked T cannot
        # block the intact one because veto precedes spatial selection.
        prob = np.array([[0.90, 0.01], [0.80, 0.01]])
        selected = resolve_contour_guided(xy, prob, [], [0.5, 0.5], image, lines)
        self.assertEqual([(r["x"], r["y"]) for r in selected], [(64.0, 41.0)])
        self.assertEqual(selected[0]["stem_path"]["reason"], "continuous")
        all_modes = resolve_all_modes(xy, prob, [], [0.5, 0.5], image, lines)
        self.assertEqual(
            [(r["x"], r["y"]) for r in all_modes["baseline"]],
            [(41.0, 41.0), (64.0, 41.0)],
        )
        vetoed = all_modes["audit"]["both"]["candidates"][0]
        self.assertEqual(vetoed["selection"], "rejected_contour_gap")
        self.assertTrue(any(vetoed["stem_path"]["veto_mask"]))

    def test_visible_gap_between_a_finite_stem_and_the_bar(self):
        image = np.full((81, 81), 255, np.uint8)
        image[40, 8:68] = 0
        image[48:70, 40] = 0
        lines = np.array([[9, 41, 68, 41], [41, 49, 41, 70]], float)
        tracks, clean = build_tracks(lines)
        evidence = t_stem_evidence(image, [41, 41], tracks, clean)
        self.assertTrue(evidence["veto"])
        # A finite detector fragment can also stop short of a continuous photo edge.
        image[41:49, 40] = 0
        self.assertFalse(t_stem_evidence(image, [41, 41], tracks, clean)["veto"])

    def test_close_xs_same_tracks_are_duplicates_but_distinct_tracks_remain(self):
        image = np.full((90, 90), 255, np.uint8)
        # Pair 1 intersects at (31, 31); pair 2 intersects two pixels away but
        # uses two genuinely different tracks.  Thick ink only helps contour
        # rendering; X identity comes from observed global line tracks.
        lines = np.array(
            [[8, 31, 54, 31], [31, 8, 31, 54], [8, 33, 54, 33], [33, 8, 33, 54]],
            float,
        )
        for x1, y1, x2, y2 in lines.astype(int):
            if y1 == y2:
                image[y1 - 3 : y1 + 2, x1 - 1 : x2] = 0
            else:
                image[y1 - 1 : y2, x1 - 3 : x1 + 2] = 0
        xy = np.array([[31, 31], [31.6, 31], [33, 33]], float)
        prob = np.array([[0.01, 0.9], [0.01, 0.8], [0.01, 0.7]])
        selected = resolve_contour_guided(
            xy, prob, [], [0.5, 0.5], image, lines, radius=4
        )
        self.assertEqual(
            [(r["x"], r["y"]) for r in selected], [(31.0, 31.0), (33.0, 33.0)]
        )
        self.assertEqual(len(selected[0]["track_signature"]), 2)
        self.assertNotEqual(
            selected[0]["track_signature"], selected[1]["track_signature"]
        )

    def test_baseline_is_exact_topology_resolve(self):
        xy = np.array([[10, 10], [11, 10], [30, 30]], float)
        prob = np.array([[0.8, 0.01], [0.7, 0.01], [0.01, 0.8]])
        expected = resolve(xy, prob, [], [0.5, 0.5])
        got = resolve_contour_guided(
            xy,
            prob,
            [],
            [0.5, 0.5],
            np.zeros((40, 40)),
            np.empty((0, 4)),
            mode="baseline",
        )
        self.assertEqual(got, expected)
