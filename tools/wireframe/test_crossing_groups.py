"""Geometric and budget controls for joint whole-path selection."""

import unittest
import numpy as np
from curve_rank import Selection, predict_paths
from crossing_groups import (
    Crossing,
    arm_support,
    force_then_fill,
    select_groups,
    take_group,
)
from test_curve_rank import fixture


class Density:
    def predict(self, features):
        return features[:, 0]


def group(pool, members, support):
    return Crossing(
        "mlsd",
        0,
        np.array(members),
        np.unique(np.concatenate([pool.pixels[i] for i in members])),
        np.unique(np.concatenate([pool.tubes[i] for i in members])),
        support,
    )


class CrossingTests(unittest.TestCase):
    def test_weakest_arm_distinguishes_crossing_from_t_and_collinear_rays(self):
        y, x = np.indices((31, 37))
        site = dict(x=18, y=15, arms=[0, 90, 180, 270])
        quadrant = np.logical_xor(x >= 17, y >= 14)
        rgb = np.repeat((quadrant * 240).astype(np.uint8)[..., None], 3, axis=2)
        self.assertGreater(arm_support(rgb, site), 0.9)
        # The upper half-plane has no vertical edge leaving the centre.
        rgb[y < 14] = 100
        self.assertEqual(arm_support(rgb, site), 0)
        for angles in ([0, 0, 180, 180], [0, 90, 180]):
            self.assertEqual(arm_support(np.ones_like(rgb), dict(site, arms=angles)), 0)
        self.assertEqual(arm_support(rgb, dict(site, x=1, y=1)), 0)

    def test_shared_pixel_is_charged_once(self):
        pool = fixture()
        pool.obs[:, 0] = [0.04, 0.04, 0.7]
        candidates = [group(pool, [0, 1], 1)]
        baseline = predict_paths(pool, Density())
        self.assertEqual(baseline.order[0], 2)
        state, n = select_groups(pool, Density(), candidates)
        self.assertEqual(n, 1)
        self.assertEqual(state.order, [0, 1])
        self.assertEqual(state.remaining, 0)
        self.assertEqual(state.used.sum(), 5)  # Six member pixels, one shared.
        unsupported, _ = select_groups(pool, Density(), [group(pool, [0, 1], 0)])
        self.assertEqual(unsupported.order, baseline.order)
        pool.budget = 4
        blocked, n = select_groups(pool, Density(), candidates)
        self.assertEqual(n, 0)
        self.assertEqual(blocked.order, [2])

    def test_forced_groups_are_atomic_and_deterministic(self):
        pool = fixture()
        crossing = group(pool, [0, 1], 0)
        state, n = force_then_fill(pool, Density(), [crossing])
        self.assertEqual(state.order, [0, 1])
        self.assertEqual(n, 1)
        pool.budget = 4
        forced, n = force_then_fill(pool, Density(), [crossing])
        self.assertEqual(n, 0)
        self.assertEqual(forced.order, predict_paths(pool, Density()).order)

    def test_covered_member_can_be_added_once_without_charging_ink(self):
        pool = fixture()
        state = Selection(pool)
        known = np.zeros(3, bool)
        # Model a path whose pixels have already been covered by other geometry.
        state.used[pool.pixels[0]] = True
        state.covered[pool.tubes[0]] = True
        state.ink[0] = 0
        state.area[0] = 0
        state.remaining = 0
        take_group(state, [0], known)
        take_group(state, [0], known)
        self.assertEqual(state.order, [0])
        self.assertEqual(state.remaining, 0)
        with self.assertRaises(ValueError):
            state.take(0, allow_covered=True)


if __name__ == "__main__":
    unittest.main()
