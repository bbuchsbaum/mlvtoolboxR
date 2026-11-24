import unittest

from score_crossing_groups import four_arm_pairs, gates


class CrossingGroupScoreTest(unittest.TestCase):
    @staticmethod
    def crossing(vectors, arms, observed_vectors=None):
        centre = [20, 20]
        artist_paths = [
            [centre, [20 + 12 * vector[0], 20 + 12 * vector[1]]] for vector in vectors
        ]
        observed_paths = [
            [centre, [20 + 12 * vector[0], 20 + 12 * vector[1]]]
            for vector in (observed_vectors or vectors)
        ]
        artist = {
            "junctions": [
                {
                    "type": "X",
                    "xy": centre,
                    "incident_directions": vectors,
                    "incident_polyline_ids": [str(i) for i in range(4)],
                }
            ],
            "polylines": [
                {"id": str(i), "points": path} for i, path in enumerate(artist_paths)
            ],
        }
        observed = {
            "sites": [
                {
                    "type": "X",
                    "x": 20,
                    "y": 20,
                    "arms": arms,
                    "members": [str(i) for i in range(4)],
                }
            ],
            "contours": [
                {"id": str(i), "points": path} for i, path in enumerate(observed_paths)
            ],
        }
        return artist, observed

    def test_four_arm_pairs_rejects_two_horizontal_rays_for_a_shallow_x(self):
        vectors = [[1, 0.2], [-1, -0.2], [1, -0.2], [-1, 0.2]]
        artist, observed = self.crossing(
            vectors,
            [0, 0, 180, 180],
            [[1, 0], [1, 0], [-1, 0], [-1, 0]],
        )
        self.assertEqual(four_arm_pairs(artist, observed), [])

    def test_four_arm_pairs_accepts_actual_orthogonal_x(self):
        artist, observed = self.crossing(
            [[1, 0], [0, -1], [-1, 0], [0, 1]], [0, 90, 180, 270]
        )
        self.assertEqual(four_arm_pairs(artist, observed), [(0, 0)])

    def test_gate_uses_four_arm_x_gain_on_two_photos(self):
        def row(method, ident, strict):
            return {
                "id": ident,
                "method": method,
                "output_status": "ok",
                "ink": 10,
                "artist_ink": 10,
                "pixel_precision": 0.9 if method == "fdog" else 1.0,
                "pixel_coverage": 0.9 if method == "fdog" else 1.0,
                "pixel_F1": 0.9 if method == "fdog" else 1.0,
                "T_ref": 1,
                "X_ref": 1,
                "T_position_pred": 1,
                "X_position_pred": 1,
                "T_position_match": 1,
                "X_position_match": 1,
                "T_incidence_match": 1,
                "X_incidence_match": 1,
                "X_four_arm_match": strict,
            }

        rows = []
        for ident in ("a", "b"):
            rows.extend(
                [
                    row("fdog", ident, 0),
                    row("baseline", ident, 0),
                    row("grouped", ident, 1),
                    row("diagnostic", ident, 9),
                ]
            )
        result = gates(rows)
        self.assertTrue(result["gate_pass"])
        self.assertFalse(result["diagnostic_qualifies"])

    def test_gate_rejects_a_two_match_gain_concentrated_on_one_photo(self):
        def row(method, ident, strict):
            return {
                "id": ident,
                "method": method,
                "output_status": "ok",
                "ink": 10,
                "artist_ink": 10,
                "pixel_precision": 0.9 if method == "fdog" else 1.0,
                "pixel_coverage": 0.9 if method == "fdog" else 1.0,
                "pixel_F1": 0.9 if method == "fdog" else 1.0,
                "T_ref": 1,
                "X_ref": 2,
                "T_position_pred": 1,
                "X_position_pred": 2,
                "T_position_match": 1,
                "X_position_match": 2,
                "T_incidence_match": 1,
                "X_incidence_match": 2,
                "X_four_arm_match": strict,
            }

        rows = []
        for ident, grouped_strict in (("a", 2), ("b", 0)):
            rows.extend(
                [
                    row("fdog", ident, 0),
                    row("baseline", ident, 0),
                    row("grouped", ident, grouped_strict),
                    row("diagnostic", ident, 9),
                ]
            )
        result = gates(rows)
        self.assertFalse(result["gate_pass"])
        self.assertFalse(
            result["checks"]["strict_X_four_arm_gain_on_at_least_2_photos"]
        )


if __name__ == "__main__":
    unittest.main()
