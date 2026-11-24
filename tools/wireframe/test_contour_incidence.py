"""Synthetic geometry controls for independent contour-incidence scoring."""

import unittest

from contour_incidence import incident_correspondence


def artist():
    return {
        "polylines": [
            {"id": 1, "points": [[1, 10], [19, 10]]},
            {"id": 2, "points": [[10, 10], [10, 19]]},
        ],
        "junctions": [
            {
                "type": "T",
                "xy": [10, 10],
                "incident_polyline_ids": [1, 2],
                "incident_directions": [[-1, 0], [1, 0], [0, 1]],
            }
        ],
    }


def observed(directions=None, shifted=False):
    x = 13 if shifted else 10
    return {
        "contours": [
            {"id": "bar", "segments": [[1, 10, 19, 10]]},
            {"id": "stem", "segments": [[x, 10, x, 19]]},
        ],
        "sites": [
            {
                "type": "T",
                "x": x,
                "y": 10,
                "members": [{"contour_id": "bar"}, {"contour_id": "stem"}],
                "arms": directions if directions is not None else [180, 0, 270],
            }
        ],
    }


class TestContourIncidence(unittest.TestCase):
    def test_named_incident_geometry_matches_once(self):
        data = observed()
        score = incident_correspondence(artist(), data, data)
        self.assertEqual(score["T"]["matched"], 1)
        self.assertEqual(score["X"]["matched"], 0)
        self.assertTrue(score["T"]["evidence"][0]["artist_covered"])
        self.assertTrue(score["T"]["evidence"][0]["observed_covered"])

    def test_opposite_stem_is_not_a_direction_match(self):
        data = observed(directions=[180, 0, 90])
        score = incident_correspondence(artist(), data, data)
        self.assertEqual(score["T"]["matched"], 0)
        self.assertEqual(score["T"]["evidence"][0]["reason"], "direction_mismatch")

    def test_nearby_parallel_member_geometry_does_not_count(self):
        data = observed(shifted=True)
        score = incident_correspondence(artist(), data, data)
        self.assertEqual(score["T"]["matched"], 0)
        self.assertEqual(
            score["T"]["evidence"][0]["reason"], "annulus_coverage_mismatch"
        )


if __name__ == "__main__":
    unittest.main()
