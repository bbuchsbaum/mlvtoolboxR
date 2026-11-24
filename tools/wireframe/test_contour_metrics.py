import unittest

import numpy as np
from contour_metrics import covered, draw_paths, junction_counts
from contour_render import result_segments
from contour_york import drawings_bridge


class ContourMetrics(unittest.TestCase):
    def test_one_predicted_crossing_at_two_references(self):
        self.assertEqual(junction_counts([[3, 3], [4, 3]], [[4, 4]]), 1)
        self.assertEqual(junction_counts([[3, 3]], [[7, 3]]), 1)
        self.assertEqual(junction_counts([[3, 3]], [[7.1, 3]]), 0)

    def test_arrow_directions_and_wraparound(self):
        site = [[10, 10]]
        a = np.deg2rad([90, 180, 359])
        b = np.deg2rad([-1, 90, 180])
        self.assertEqual(junction_counts(site, site, [a], [b]), 1)
        # The bar aligns, but the stem points to the other side of the bar.
        self.assertEqual(
            junction_counts(site, site, [a], [np.deg2rad([180, 270, 0])]), 0
        )

    def test_stem_cut_preserves_bar_and_distant_stem(self):
        lines = np.array([[1.0, 5.0, 15.0, 5.0], [7.0, 5.0, 7.0, 15.0]])
        track = {
            "track_id": 1,
            "segments": [
                {"segment_id": 1, "polyline": lines[1].reshape(2, 2).tolist()}
            ],
        }
        rejected = {
            "stem_path": {
                "reason": "contour_gap",
                "track_id": 1,
                "incident_tracks": [track],
                "sample_xy": [[7, 7], [7, 8], [7, 12]],
                "veto_mask": [True, True, False],
            }
        }
        seg = result_segments(lines, [], [rejected])
        mask = draw_paths((16, 16), seg.reshape(-1, 2, 2))
        self.assertTrue(mask[4, :15].all())
        self.assertTrue(mask[11, 6])
        self.assertFalse(mask[7, 6])

    def test_mask_tolerance_and_connected_stem(self):
        stem = np.array([7.0, 15.0])
        end = np.array([7.0, 10.0])
        bar = np.array([[1.0, 5.0], [13.0, 5.0]])
        gap = draw_paths((16, 16), [[stem, end], bar])
        closed = draw_paths((16, 16), [[stem, [7, 5]], bar])
        self.assertFalse(drawings_bridge(gap, stem, end, bar))
        self.assertTrue(drawings_bridge(closed, stem, end, bar))
        self.assertEqual(covered(gap, closed), gap.sum())
        self.assertEqual(covered(closed, np.zeros_like(gap)), 0)


if __name__ == "__main__":
    unittest.main()
