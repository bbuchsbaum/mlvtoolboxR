"""Focused checks for the Toronto Scenes private-staging contract."""

import hashlib
import importlib.util
from pathlib import Path

import numpy as np


MODULE = Path(__file__).with_name("toronto_scenes.py")
SPEC = importlib.util.spec_from_file_location("toronto_scenes", MODULE)
TORONTO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TORONTO)


def test_per_category_hash_split_is_complete_and_disjoint():
    records = []
    for category in ("cities", "forests"):
        for number in range(80):
            records.append(
                {
                    "category": category,
                    "row": {"originalImage": f"{category}_{number:05d}.jpg"},
                }
            )
    assigned = TORONTO.assign_splits(records)
    for category in ("cities", "forests"):
        subset = [row for row in assigned if row["category"] == category]
        assert [row["panel"] for row in subset] == ["development"] * 10 + [
            "confirmation"
        ] * 15 + ["reserve"] * 55
        assert [row["rank"] for row in subset] == list(range(1, 81))
        assert [row["split_key"] for row in subset] == sorted(
            row["split_key"] for row in subset
        )


def test_pixel_center_scaling_and_crossing_rays():
    assert TORONTO.scaled_point([1, 1], (800, 600), (256, 192)) == [0.66, 0.66]
    junction = {
        "position": np.array([5.0, 5.0]),
        "contourIDs": np.array([1, 2]),
        "type": "X",
    }
    contours = np.array(
        [
            np.array([[0.0, 5.0, 10.0, 5.0]]),
            np.array([[5.0, 0.0, 5.0, 10.0]]),
        ],
        dtype=object,
    )
    rays = TORONTO.junction_directions(junction, contours)
    assert len(rays) == 4
    assert sorted(round(abs(ray[0]) + abs(ray[1]), 8) for ray in rays) == [1.0] * 4


def test_split_key_is_explicit_and_stable():
    image_id = "cities_00123.jpg"
    assert (
        TORONTO.split_key("cities", image_id)
        == hashlib.sha256(
            (TORONTO.SPLIT_SALT + "cities/" + image_id).encode()
        ).hexdigest()
    )


if __name__ == "__main__":
    test_per_category_hash_split_is_complete_and_disjoint()
    test_pixel_center_scaling_and_crossing_rays()
    test_split_key_is_explicit_and_stable()
    print("Toronto Scenes focused checks passed")
