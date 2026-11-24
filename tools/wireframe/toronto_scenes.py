"""Stage a frozen, private Toronto Scenes photo-and-drawing benchmark panel.

The OSF archive and the local MATLAB vectors are inputs only.  This program
does not place either source asset under version control: the public evidence
file contains identities, hashes, dimensions and split assignments, while
pixels and vector geometry remain in the requested private staging directory.
"""

import argparse
import csv
import hashlib
import json
import math
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from PIL import Image, ImageOps
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[2]
MAT_DIR = ROOT / "reference-data/mlv/dataSets/TorontoScenes"
ARCHIVE_SHA256 = "c2f34ec5c7085baeb8c647a0886a74ee504cffb211ece6114b241d6127c396e2"
SPLIT_SALT = "toronto-scenes-photo-split-v1/"
SCHEMA_VERSION = 1
ARTIST_FILES = {
    "beaches_vecLD.mat": "9ca35525b267daaab189daa05d2a42e8ac335a97d2b3451d99e0799ed6a1949c",
    "cities_vecLD.mat": "dd840897655cd6385e4ac118b8ce0e01303291e5982d2efafadaae8c845d63db",
    "forests_vecLD.mat": "4a65aadbf4190736cf949f20baf5b2aaa2b57b3c2fe8cebbf5e7816bcc3c7445",
    "highways_vecLD.mat": "c8c040eb7f96264f0b9e115140396b1d95cc6d97fc9532f01bfd1a4fe4aaa886",
    "mountains_vecLD.mat": "d17a386c4e9263fbc85d6eb18e4a63d5349473d32ecd4631149b3f5c3ae538dd",
    "offices_vecLD.mat": "0655e33c8c8f064a3e0454f71cfe5561c5295429d3d4deea3d81822677a3424a",
}
PERCEPTUAL_HAMMING_THRESHOLD = 4


def resize_photo(rgb):
    """Use the Wireframe working-image convention: 256 pixels on long side."""
    size = tuple(max(1, round(d * 256 / max(rgb.size))) for d in rgb.size)
    return rgb.resize(size, Image.Resampling.LANCZOS)


def file_sha256(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def split_key(category, image_id):
    return hashlib.sha256((SPLIT_SALT + category + "/" + image_id).encode()).hexdigest()


def scaled_point(point, original_size, working_size):
    """Map one-based source pixel coordinates onto one-based working centers."""
    return [
        (float(point[0]) - 0.5) * working_size[0] / original_size[0] + 0.5,
        (float(point[1]) - 0.5) * working_size[1] / original_size[1] + 0.5,
    ]


def normalise(vector):
    length = math.hypot(*vector)
    if length == 0:
        return None
    return [vector[0] / length, vector[1] / length]


def point_segment_distance(point, left, right):
    segment = right - left
    denom = float(segment @ segment)
    if denom == 0:
        return float(np.linalg.norm(point - left))
    t = max(0.0, min(1.0, float((point - left) @ segment / denom)))
    return float(np.linalg.norm(point - (left + t * segment)))


def unique_directions(directions, radians=math.radians(2)):
    """Deduplicate coincident raw contour rays without inventing directions."""
    answer = []
    for direction in directions:
        if direction is None:
            continue
        if not any(
            float(np.dot(direction, old)) >= math.cos(radians) for old in answer
        ):
            answer.append(direction)
    return answer


def junction_directions(junction, contours):
    """Derive outward incident rays from the artist's contour segments.

    The MAT `angles` field contains angular gaps, not absolute directions.
    We therefore recover directions from each named contour at the stored
    junction position.  A contour crossing a junction supplies two rays.
    """
    center = np.asarray(junction["position"], dtype=float)
    rays = []
    for contour_id in np.asarray(junction["contourIDs"]).ravel():
        contour = np.atleast_2d(np.asarray(contours[int(contour_id) - 1], dtype=float))
        for x1, y1, x2, y2 in contour:
            left, right = np.array([x1, y1]), np.array([x2, y2])
            if point_segment_distance(center, left, right) > 1.5:
                continue
            # A segment containing the site contributes its two visible arms;
            # a segment terminating there contributes one.
            for endpoint in (left, right):
                ray = normalise(endpoint - center)
                if ray is not None:
                    rays.append(ray)
    return unique_directions(rays)


def polyline_points(contour):
    contour = np.atleast_2d(np.asarray(contour))
    if len(contour) == 0:
        return []
    return [contour[0, :2].astype(float).tolist()] + [
        segment[2:].astype(float).tolist() for segment in contour
    ]


def load_records(mat_dir=MAT_DIR):
    records = []
    for name, digest in ARTIST_FILES.items():
        path = Path(mat_dir) / name
        if file_sha256(path) != digest:
            raise ValueError(
                f"Artist vector file checksum differs from the freeze: {name}"
            )
        category = path.name.removesuffix("_vecLD.mat")
        data = loadmat(path, simplify_cells=True)["vecLD"]
        for row in data:
            records.append({"category": category, "row": row})
    return records


def assign_splits(records):
    by_category = {}
    for record in records:
        image_id = record["row"]["originalImage"]
        record["id"] = image_id
        by_category.setdefault(record["category"], []).append(record)
    for category, category_records in by_category.items():
        category_records.sort(key=lambda record: split_key(category, record["id"]))
        for rank, record in enumerate(category_records, start=1):
            record["rank"] = rank
            record["split_key"] = split_key(category, record["id"])
            record["panel"] = (
                "development"
                if rank <= 10
                else "confirmation" if rank <= 25 else "reserve"
            )
    return sorted(records, key=lambda record: (record["category"], record["rank"]))


def geometry(record, original_size, working_size, archive_member):
    row = record["row"]
    contours = row["contours"]
    junctions = row["junctions"]
    if isinstance(junctions, dict):
        junctions = [junctions]
    return {
        "schema_version": SCHEMA_VERSION,
        "id": record["id"],
        "category": record["category"],
        "panel": record["panel"],
        "source": {
            "archive_sha256": ARCHIVE_SHA256,
            "member": archive_member,
            "original_width": original_size[0],
            "original_height": original_size[1],
        },
        "working": {"width": working_size[0], "height": working_size[1]},
        "polylines": [
            {
                "id": index + 1,
                "points": [
                    scaled_point(point, original_size, working_size)
                    for point in polyline_points(contour)
                ],
            }
            for index, contour in enumerate(contours)
        ],
        "junctions": [
            {
                "id": index + 1,
                "type": str(junction["type"]),
                "xy": scaled_point(junction["position"], original_size, working_size),
                "incident_polyline_ids": [
                    int(value) for value in np.asarray(junction["contourIDs"]).ravel()
                ],
                "incident_directions": [
                    normalise(direction)
                    for direction in junction_directions(junction, contours)
                ],
            }
            for index, junction in enumerate(junctions)
        ],
    }


def contour_edge_score(image, contours, original_size, working_size, offset=(0, 0)):
    """Mean local RGB-gradient magnitude at sampled artist contour points."""
    pixels = np.asarray(image, dtype=float)
    gray = pixels @ np.array([0.2126, 0.7152, 0.0722])
    gy, gx = np.gradient(gray)
    gradient = np.hypot(gx, gy)
    values = []
    width, height = working_size
    for contour in contours:
        points = [
            scaled_point(point, original_size, working_size)
            for point in polyline_points(contour)
        ]
        for left, right in zip(points, points[1:]):
            # One sample per eight working pixels is enough for a coarse
            # registration sanity check and keeps this preparation bounded.
            length = max(1, math.ceil(math.dist(left, right) / 8))
            for t in np.linspace(0, 1, length + 1):
                x = round(left[0] * (1 - t) + right[0] * t - 1) + offset[0]
                y = round(left[1] * (1 - t) + right[1] * t - 1) + offset[1]
                if 1 <= x < width - 1 and 1 <= y < height - 1:
                    values.append(float(gradient[y - 1 : y + 2, x - 1 : x + 2].max()))
    return float(np.mean(values)) if values else float("nan")


def read_panel_hashes(evidence_dir):
    hashes = {}
    for panel in sorted(Path(evidence_dir).glob("*panel.csv")):
        with panel.open(newline="") as source:
            hashes[panel.name] = {row["pixel_sha256"] for row in csv.DictReader(source)}
    return hashes


def write_checks(path, registration, source_hashes, panel_hashes):
    overlaps = {
        panel: len(source_hashes & hashes) for panel, hashes in panel_hashes.items()
    }
    report = {
        "archive_sha256": ARCHIVE_SHA256,
        "registration": registration,
        "exact_resized_rgb_hash_overlap": overlaps,
        "perceptual_overlap": {
            "status": "not_computed",
            "reason": "Run --audit-inputs with the local comparison photos for the perceptual search.",
        },
    }
    Path(path).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


def difference_hash(rgb):
    """A fixed 64-bit dHash for input-overlap screening, not similarity scoring."""
    gray = rgb.convert("L").resize((9, 8), Image.Resampling.LANCZOS)
    pixels = np.asarray(gray, dtype=np.uint8)
    bits = pixels[:, 1:] > pixels[:, :-1]
    return sum(int(bit) << index for index, bit in enumerate(bits.ravel()))


def photo_fingerprint(image):
    rgb = ImageOps.exif_transpose(image).convert("RGB")
    return {
        "pixel_sha256": hashlib.sha256(rgb.tobytes()).hexdigest(),
        "dhash64": difference_hash(rgb),
    }


def hamming(left, right):
    return (left ^ right).bit_count()


def overlap_summary(left, right):
    """Summarise all cross-set pairs without retaining any source pixels."""
    exact = 0
    near = 0
    nearest = None
    for left_row in left:
        for right_row in right:
            distance = hamming(left_row["dhash64"], right_row["dhash64"])
            exact += left_row["pixel_sha256"] == right_row["pixel_sha256"]
            near += distance <= PERCEPTUAL_HAMMING_THRESHOLD
            candidate = (distance, left_row["id"], right_row["id"])
            if nearest is None or candidate < nearest:
                nearest = candidate
    return {
        "left_count": len(left),
        "right_count": len(right),
        "exact_rgb_pairs": exact,
        "dhash_pairs_at_or_below_threshold": near,
        "nearest_dhash_distance": nearest[0] if nearest else None,
        "nearest_pair": list(nearest[1:]) if nearest else None,
    }


def input_audit(
    archive, compare_root, output, public_evidence=None, registration_check=None
):
    """Audit only RGB inputs; does not write or modify the staged benchmark."""
    archive = Path(archive)
    compare_root = Path(compare_root)
    if file_sha256(archive) != ARCHIVE_SHA256:
        raise ValueError(
            "Toronto archive checksum differs from the official pinned archive"
        )
    records = assign_splits(load_records())
    toronto = []
    with ZipFile(archive) as source:
        for record in records:
            member = "CP_images/" + record["id"]
            with Image.open(BytesIO(source.read(member))) as image:
                small = resize_photo(ImageOps.exif_transpose(image).convert("RGB"))
                fingerprint = photo_fingerprint(small)
            toronto.append(
                {"id": record["id"], "panel": record["panel"], **fingerprint}
            )
    locations = {
        "wireframe-root": compare_root,
        "topology-train": compare_root / "topology-train",
        "prospective": compare_root / "prospective",
        "remainder": compare_root / "remainder",
        "confirmation": compare_root / "confirmation",
        "replication": compare_root / "replication",
        "york-fresh": compare_root / "arms-next/fresh",
        "york51": compare_root / "york51",
    }
    external = {}
    for name, folder in locations.items():
        if not folder.is_dir():
            raise ValueError(f"missing requested comparison folder: {folder}")
        rows = []
        for path in sorted(folder.glob("*.small.png")):
            with Image.open(path) as image:
                rows.append(
                    {
                        "id": path.name.removesuffix(".small.png"),
                        **photo_fingerprint(image),
                    }
                )
        external[name] = rows
    groups = {
        panel: [row for row in toronto if row["panel"] == panel]
        for panel in ("development", "confirmation", "reserve")
    }
    report = {
        "archive_sha256": ARCHIVE_SHA256,
        "dhash": {
            "bits": 64,
            "threshold": PERCEPTUAL_HAMMING_THRESHOLD,
            "algorithm": "Pillow LANCZOS 9x8 grayscale adjacent-horizontal comparison",
        },
        "toronto": {
            "total": len(toronto),
            "panels": {key: len(value) for key, value in groups.items()},
        },
        "source_paths": {key: str(value) for key, value in locations.items()},
        "intra_toronto": {
            "development_confirmation": overlap_summary(
                groups["development"], groups["confirmation"]
            ),
            "development_reserve": overlap_summary(
                groups["development"], groups["reserve"]
            ),
            "confirmation_reserve": overlap_summary(
                groups["confirmation"], groups["reserve"]
            ),
        },
        "external": {
            name: overlap_summary(toronto, rows) for name, rows in external.items()
        },
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    if public_evidence:
        if not registration_check or not Path(registration_check).is_file():
            raise ValueError(
                "public input evidence requires the private registration check"
            )
        registration = json.loads(Path(registration_check).read_text())["registration"]
        public = {
            "artist_file_sha256": ARTIST_FILES,
            "archive_sha256": ARCHIVE_SHA256,
            "registration": registration,
            "dhash": report["dhash"],
            "toronto": report["toronto"],
            "panels": {
                name: {
                    "relative_path": (
                        "."
                        if name == "wireframe-root"
                        else name.replace("york-fresh", "arms-next/fresh")
                    ),
                    "image_count": value["right_count"],
                    "exact_rgb_pairs": value["exact_rgb_pairs"],
                    "dhash_pairs_at_or_below_threshold": value[
                        "dhash_pairs_at_or_below_threshold"
                    ],
                    "nearest_dhash_distance": value["nearest_dhash_distance"],
                }
                for name, value in report["external"].items()
            },
            "intra_toronto": {
                name: {
                    key: value[key]
                    for key in (
                        "left_count",
                        "right_count",
                        "exact_rgb_pairs",
                        "dhash_pairs_at_or_below_threshold",
                        "nearest_dhash_distance",
                    )
                }
                for name, value in report["intra_toronto"].items()
            },
        }
        public_evidence = Path(public_evidence)
        public_evidence.parent.mkdir(parents=True, exist_ok=True)
        public_evidence.write_text(json.dumps(public, indent=2) + "\n")


def prepare(archive, stage, evidence):
    archive = Path(archive)
    stage = Path(stage)
    evidence = Path(evidence)
    if file_sha256(archive) != ARCHIVE_SHA256:
        raise ValueError(
            "Toronto archive checksum differs from the official pinned archive"
        )
    records = assign_splits(load_records())
    if len(records) != 475 or len({record["id"] for record in records}) != 475:
        raise ValueError("expected 475 unique Toronto Scenes image IDs")
    if sum(record["panel"] == "development" for record in records) != 60:
        raise ValueError("expected 60 development images")
    if sum(record["panel"] == "confirmation" for record in records) != 90:
        raise ValueError("expected 90 confirmation images")
    stage.mkdir(parents=True, exist_ok=True)
    evidence.parent.mkdir(parents=True, exist_ok=True)
    all_rows, selected_rows, registration = [], [], []
    with ZipFile(archive) as source:
        names = set(source.namelist())
        for record in records:
            member = "CP_images/" + record["id"]
            if member not in names:
                raise ValueError(f"archive does not contain {member}")
            with Image.open(BytesIO(source.read(member))) as image:
                rgb = ImageOps.exif_transpose(image).convert("RGB")
            original_size = rgb.size
            expected_size = tuple(
                int(value) for value in np.asarray(record["row"]["imsize"]).ravel()
            )
            if original_size != expected_size:
                raise ValueError(f"MAT/image dimensions disagree for {record['id']}")
            small = resize_photo(rgb)
            digest = hashlib.sha256(small.tobytes()).hexdigest()
            aligned = contour_edge_score(
                small, record["row"]["contours"], original_size, small.size
            )
            shifted = contour_edge_score(
                small,
                record["row"]["contours"],
                original_size,
                small.size,
                offset=(17, 11),
            )
            registration.append({"aligned": aligned, "shifted": shifted})
            common = dict(
                id=record["id"],
                panel=record["panel"],
                category=record["category"],
                rank=record["rank"],
                split_sha256=record["split_key"],
                original_width=original_size[0],
                original_height=original_size[1],
                width=small.width,
                height=small.height,
                pixel_sha256=digest,
            )
            all_rows.append(common)
            if record["panel"] != "reserve":
                small.save(stage / (record["id"] + ".small.png"))
                doc = geometry(record, original_size, small.size, member)
                doc["working"]["pixel_sha256"] = digest
                (stage / (record["id"] + ".json")).write_text(
                    json.dumps(doc, separators=(",", ":"), allow_nan=False) + "\n"
                )
                selected_rows.append(
                    {
                        key: common[key]
                        for key in (
                            "id",
                            "panel",
                            "category",
                            "width",
                            "height",
                            "pixel_sha256",
                        )
                    }
                )
    write_csv(
        evidence,
        all_rows,
        [
            "id",
            "category",
            "panel",
            "rank",
            "split_sha256",
            "original_width",
            "original_height",
            "width",
            "height",
            "pixel_sha256",
        ],
    )
    write_csv(
        stage / "panel.csv",
        selected_rows,
        ["id", "panel", "category", "width", "height", "pixel_sha256"],
    )
    aligned = np.asarray([row["aligned"] for row in registration])
    shifted = np.asarray([row["shifted"] for row in registration])
    write_checks(
        stage / "checks.json",
        {
            "dimension_matches": len(registration),
            "shift_pixels": [17, 11],
            "mean_aligned_gradient": float(aligned.mean()),
            "mean_shifted_gradient": float(shifted.mean()),
            "photos_aligned_gt_shifted": int((aligned > shifted).sum()),
        },
        {row["pixel_sha256"] for row in all_rows},
        read_panel_hashes(ROOT / "docs/evidence/wireframe"),
    )


def write_csv(path, rows, fields):
    with Path(path).open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument(
        "--stage", type=Path, default=ROOT / ".agent-work/toronto/staged"
    )
    parser.add_argument(
        "--evidence",
        type=Path,
        default=ROOT / "docs/evidence/wireframe/toronto/photos.csv",
    )
    parser.add_argument(
        "--audit-inputs",
        action="store_true",
        help="compare archive RGB against the supplied Wireframe/York input root without changing staged files; writes private and public-safe summaries",
    )
    parser.add_argument(
        "--compare-root",
        type=Path,
        help="root containing the named Wireframe/York input folders",
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=ROOT / ".agent-work/toronto/input-checks.json",
    )
    parser.add_argument(
        "--public-evidence",
        type=Path,
        default=ROOT / "docs/evidence/wireframe/toronto/input-checks.json",
    )
    parser.add_argument(
        "--registration-check",
        type=Path,
        default=ROOT / ".agent-work/toronto/staged/checks.json",
    )
    args = parser.parse_args()
    if args.audit_inputs:
        if not args.compare_root:
            parser.error("--audit-inputs requires --compare-root")
        input_audit(
            args.archive,
            args.compare_root,
            args.audit_output,
            args.public_evidence,
            args.registration_check,
        )
        return
    prepare(args.archive, args.stage, args.evidence)


if __name__ == "__main__":
    main()
