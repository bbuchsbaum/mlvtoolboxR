"""Analytic controls for the learned-topology feature and selection paths."""
import unittest

import numpy as np

from fit_topology import negative_sites
from topology import (
    axis_features,
    geometry,
    intersections,
    model_features,
    probabilities,
    resolve,
)


class TestVectorTopology(unittest.TestCase):
    def test_context_uses_pixel_centers_and_separate_neighborhoods(self):
        import torch
        from context_features import SCENE_COLUMNS, descriptors, neighborhood
        from predict_hawp import sample_context

        # The map covers four photo pixels in each direction. The photo center
        # is equally far from all four map samples.
        maps = torch.tensor([[[[1.0, 3.0], [5.0, 7.0]]]])
        xy = np.array([[1, 1], [2.5, 2.5], [4, 4]])
        np.testing.assert_array_equal(sample_context(maps, xy, (4, 4))[:, 0], [1, 4, 7])
        empty = descriptors(np.empty((0, 140)), maps, xy[:0], (4, 4))
        self.assertEqual(empty.shape, (0, 230))
        xy = np.array([[1, 1], [5, 1], [24, 1], [70, 1]])
        geometry = np.zeros((4, 140), np.float32)
        geometry[:, SCENE_COLUMNS[0]] = [1, 3, 8, 12]
        n = len(SCENE_COLUMNS)
        scene = neighborhood(xy, geometry, (100, 10))
        np.testing.assert_allclose(
            scene[0, [0, n, 2 * n + 1, 3 * n + 1, 4 * n + 2]], [2, 3, 4, 8, 6]
        )

    def test_arm_geometry_and_gap(self):
        # Values are analytic for three or four perpendicular unit-strength rays.
        hist = np.zeros((2, 36))
        hist[:, [0, 9, 18]] = 1
        hist[1, 27] = 1
        for shift in (0, 2, 11):
            measured = axis_features(np.roll(hist, shift, axis=1))
            np.testing.assert_array_equal(measured[:, 2:5], [[1, 0, 1], [1, 1, 0]])
        lines = np.array([[10, 30, 50, 30, 1], [30, 7, 30, 50, 1]], float)
        np.testing.assert_allclose(intersections(geometry(lines)), [[30, 30]])
        lines[1, 3] = 24
        self.assertEqual(len(intersections(geometry(lines))), 0)
        np.testing.assert_allclose(
            negative_sites(lines[:, :4].reshape(-1, 2, 2), []), [[30, 30]]
        )
        self.assertEqual(
            len(negative_sites(lines[:, :4].reshape(-1, 2, 2), [("X", [31, 30])])), 0
        )

    def test_empty_and_score_export(self):
        wire = np.empty((0, 5))
        image = np.zeros((40, 60, 3), np.uint8)
        for variant, detail, n in [
            ("arms", None, 75),
            ("photoarms", None, 103),
            ("multiscale", wire, 140),
        ]:
            xy, f = model_features(image, wire, variant, detail)
            self.assertEqual(xy.shape, (0, 2))
            self.assertEqual(f.shape, (0, n))
        stump = dict(
            bias=0.0,
            trees=[
                dict(
                    feature=[0, 0, 0],
                    threshold=[2.0, 0.0, 0.0],
                    left=[1, 0, 0],
                    right=[2, 0, 0],
                    leaf=[0, 1, 1],
                    value=[0.0, 0.0, np.log(3)],
                )
            ],
        )
        np.testing.assert_allclose(
            probabilities(np.array([[2.0], [3.0]]), [stump])[:, 0], [0.5, 0.75]
        )

    def test_contour_and_new_site_have_distinct_thresholds(self):
        base = [dict(x=18, y=21, type="T")]
        xy = np.array([[19, 21], [40, 50], [44, 50], [70, 80]], float)
        prob = np.array([[0.2, 0.005], [0.02, 0.35], [0.01, 0.28], [0.18, 0.003]])
        hits = resolve(xy, prob, base, [0.1, 0.1, 0.3, 0.2])
        self.assertEqual(
            {(r["x"], r["y"], r["type"]) for r in hits},
            {(18.0, 21.0, "T"), (40.0, 50.0, "X")},
        )
