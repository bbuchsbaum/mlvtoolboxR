"""Focused contract tests for fixed curve observations."""

import unittest

import numpy as np

from curve_observations import COLUMN_NAMES, observations


class CurveObservations(unittest.TestCase):
    def test_shape_names_and_fortran_pixel_contract(self):
        rgb = np.full((5, 6, 3), 255, np.uint8)
        rgb[:, 3:, :] = 0
        # x=3, y=1..3: zero-based Fortran indices are y + height * x.
        pixels = [np.array([16, 17, 18]), np.array([1, 6, 11])]
        paths = [np.array([[4, 2, 4, 4]]), np.array([[2, 2, 3, 2]])]
        got = observations(rgb, pixels, paths, ["fdog", "mlsd"])
        self.assertEqual(got.shape, (2, len(COLUMN_NAMES)))
        self.assertEqual(COLUMN_NAMES[0], "source_mlsd")
        self.assertEqual(got[0, 0], 0)
        self.assertEqual(got[1, 0], 1)
        self.assertEqual(got[0, 1], 3)
        self.assertGreater(got[0, 4], 100)  # normal contrast at the vertical edge

    def test_cross_source_support_and_endpoint_contact_are_geometric(self):
        rgb = np.full((8, 8, 3), 128, np.uint8)
        # Curves share two exact raster pixels; their first endpoints are near.
        pixels = [np.array([9, 17, 25]), np.array([17, 25, 33])]
        paths = [np.array([[2, 2, 2, 4]]), np.array([[3, 2, 3, 4]])]
        got = observations(rgb, pixels, paths, ["fdog", "mlsd"])
        support = COLUMN_NAMES.index("other_source_near_fraction")
        endpoint = COLUMN_NAMES.index("endpoint_near_other_curve")
        self.assertEqual(got[0, support], 1)
        self.assertEqual(got[0, endpoint], 1)

    def test_rejects_misaligned_inputs_and_bad_indices(self):
        rgb = np.zeros((3, 3, 3), np.uint8)
        with self.assertRaisesRegex(ValueError, "same length"):
            observations(rgb, [], [np.empty((0, 4))], [])
        with self.assertRaisesRegex(ValueError, "outside"):
            observations(rgb, [np.array([9])], [np.empty((0, 4))], ["fdog"])


if __name__ == "__main__":
    unittest.main()
