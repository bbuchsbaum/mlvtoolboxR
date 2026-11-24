"""Independent checks for the ink budget and topology training signal."""

import unittest

import numpy as np
from scipy.ndimage import distance_transform_edt

from curve_rank import Pool, Selection, best, membership, predict_paths, tube


def fixture():
    shape = (4, 7)
    pixels = [np.array(p) for p in ([0, 4, 8], [8, 9, 10], [18, 22, 26])]
    tubes = [tube(p, shape) for p in pixels]
    n = np.prod(shape)
    return Pool(
        shape,
        [None] * 3,
        pixels,
        tubes,
        np.array([[0, 8], [8, 10], [18, 26]]),
        np.ones((3, 1)),
        [("T", 3, 1)],
        membership([[0, 1]], 3).T.tocsr(),
        membership(pixels, n),
        membership(tubes, n),
        5,
    )


class TestCurveRank(unittest.TestCase):
    def test_training_rewards_against_dense_masks(self):
        p = fixture()
        artist = np.arange(np.prod(p.shape)) % 3 == 0
        state = Selection(p, artist, np.array([True]))
        dense = np.zeros(p.shape, bool, order="F")
        for j, credit in [(0, 2), (1, 6)]:
            before = (
                distance_transform_edt(~dense) <= 2 if dense.any() else dense.copy()
            )
            after = dense.copy(order="F")
            after.ravel(order="F")[p.pixels[j]] = True
            reached = distance_transform_edt(~after) <= 2
            self.assertEqual(state.ink[j], np.sum(after & ~dense))
            self.assertEqual(state.area[j], np.sum(reached & ~before))
            self.assertEqual(
                state.reward()[j],
                np.sum((reached & ~before).ravel(order="F") & artist) + credit,
            )
            state.take(j)
            dense = after
            np.testing.assert_array_equal(state.used, dense.ravel(order="F"))
        self.assertEqual(state.remaining, 0)
        self.assertFalse(state.eligible().any())
        with self.assertRaises(ValueError):
            state.take(2)

    def test_photo_only_prediction_and_ties(self):
        class Flat:
            def predict(self, features):
                return np.zeros(len(features))

        p = fixture()
        chosen = predict_paths(p, Flat())
        self.assertEqual(chosen.order, [0, 1])
        self.assertEqual(chosen.remaining, 0)
        state = Selection(p)
        state.take(0)
        # Equal gain per new ink favors larger incremental gain.
        state.remaining = 3
        self.assertEqual(best(state, np.array([0, 4, 6]), state.eligible()), 2)


if __name__ == "__main__":
    unittest.main()
