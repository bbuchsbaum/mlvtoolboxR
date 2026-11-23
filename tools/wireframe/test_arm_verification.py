"""Independent pixel and accounting controls for the fixed-location experiment."""
import unittest
import numpy as np
from arm_verification import arm_profiles, exchange
from score_arm_review import region_scores


class TestArmVerification(unittest.TestCase):
    def test_three_and_four_arms_and_a_gap_at_two_scales(self):
        for n in (65, 162):
            y, x = (np.indices((n, n)) + 0.5) / (n / 65) - 32.5
            for image, strong, absent in [
                (50 * (x > 0) + 70 * (y > 0), [0, 9, 18, 27], []),
                (np.where(y < 0, 100, 50 * (x > 0)), [0, 9, 18], [27]),
                (np.where(y < 0, 100, 50 * ((x > 0) & (y > 6))), [0, 18], [9, 27]),
            ]:
                arms = arm_profiles(
                    image.astype(np.float32),
                    np.array([[33.0, 33.0]]),
                    np.array([65, 65]),
                )
                length = arms[0, :, 4]
                self.assertTrue(all(length[strong] >= 17))
                self.assertTrue(all(length[absent] <= 1))

    def test_swaps_preserve_count_and_ignore_equal_advantage(self):
        labels = np.array([True] * 10 + [False] * 8)
        scores = np.arange(18)
        flips = exchange(labels, scores)
        np.testing.assert_array_equal(np.flatnonzero(labels != flips), [0, 1, 16, 17])
        np.testing.assert_array_equal(exchange(labels, np.ones(18)), labels)

    def test_region_coordinates_and_empty_references(self):
        key = [
            dict(
                token="R01",
                id="a",
                left=17,
                top=27,
                width=640,
                height=480,
                small_width=256,
                small_height=192,
            )
        ]
        r = {"R01": {"complete": True, "marks": [dict(x=20, y=30, label="X")]}}
        events = [
            dict(
                id="a", x=15.5, y=23.5, context="X", native="T", working="X", scores="X"
            )
        ]
        out = region_scores(key, r, events)
        self.assertEqual(
            [out["context"]["X"], out["native"]["X"], out["native"]["X_reference"]],
            [1, 0, 1],
        )
        self.assertEqual(out["native"]["gaps"], 0)
