"""Geometric controls for the experimental junction additions."""
import unittest

import cv2
import numpy as np

from evaluate_support import extend, rank_crops
from fragment_crossings import crossings


def site(x, y, kind, view):
    return dict(x=x, y=y, type=kind, method=view)


class TestPhotoJunctions(unittest.TestCase):
    def test_four_arm_support(self):
        horizontal = ((16, 64), (110, 64))
        cases = {
            "blank": [],
            "parallel": [horizontal, ((16, 80), (110, 80))],
            "T": [horizontal, ((64, 16), (64, 64))],
            "gap": [horizontal, ((64, 16), (64, 58))],
            "X": [horizontal, ((64, 16), (64, 110))],
        }
        for name, lines in cases.items():
            with self.subTest(name=name):
                image = np.full((128, 128, 3), 255, np.uint8)
                for a, b in lines:
                    cv2.line(image, a, b, (20, 20, 20), 1, cv2.LINE_AA)
                hits = crossings(image)
                if name == "X":
                    self.assertTrue(hits)
                    self.assertTrue(
                        all(np.hypot(x - 65, y - 65) <= 4 for _, x, y in hits)
                    )
                else:
                    self.assertEqual(hits, [])

    def test_support_needs_distinct_views_and_keeps_baseline_sites(self):
        primary = [site(x, 10, "T", "base") for x in (10, 20, 30, 40, 50)]
        crop = [
            site(10, 10, "X", "a"),
            site(30, 30, "T", "a"),
            site(31, 30, "T", "a"),
            site(45, 45, "T", "b"),
        ]
        enhanced = [site(10, 11, "X", "full"), site(45, 44, "T", "full")]
        ranked = rank_crops(enhanced, crop)
        self.assertEqual(ranked, [crop[0], crop[3]])
        result = extend(primary, ranked, "supported")
        self.assertEqual(
            [(r["x"], r["y"], r["type"]) for r in result],
            [(r["x"], r["y"], r["type"]) for r in primary + [crop[3]]],
        )
        # Five baseline sites allow only one addition even with two distinct proposals.
        self.assertEqual(len(extend(primary, crop[1::2], "supported")), 6)
        self.assertEqual(extend([], ranked, "supported"), [])
