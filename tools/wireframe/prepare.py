"""Download pinned Wireframe test inputs and run the frozen NAVER M-LSD model."""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from tempfile import NamedTemporaryFile
from pathlib import Path
from urllib.request import urlopen
import cv2
import numpy as np
from ai_edge_litert.interpreter import Interpreter
from PIL import Image, ImageOps

DATA = "https://huggingface.co/datasets/lh9171338/Wireframe/resolve/81c53616b0a63cad5183f00345b640d346932653/test"
MODEL = "https://raw.githubusercontent.com/navervision/mlsd/453cafa09467d0272760578d35c1fda38e8895a5/tflite_models/M-LSD_512_large_fp32.tflite"
WITNESSES = ["00079423", "00078318"]


def enhance_contrast(rgb):
    """Equalize local lightness on the resized photo before model inference."""
    lab = cv2.cvtColor(rgb, cv2.COLOR_RGB2LAB)
    equalizer = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    lab[:, :, 0] = equalizer.apply(lab[:, :, 0])
    return cv2.cvtColor(lab, cv2.COLOR_LAB2RGB)


def download(path, url):
    if not path.exists():
        with urlopen(url, timeout=60) as response:
            data = response.read()
        with NamedTemporaryFile(dir=path.parent, delete=False) as temp:
            temp.write(data)
        Path(temp.name).replace(path)
    return path


def run_model(model, image):
    # Upstream utils.pred_lines input transform and default thresholds.
    resized = cv2.resize(image, (512, 512), interpolation=cv2.INTER_AREA)
    tensor = np.ones((1, 512, 512, 4), dtype=np.float32)
    tensor[0, :, :, :3] = resized
    model.set_tensor(model.get_input_details()[0]["index"], tensor)
    model.invoke()
    centers, scores, offsets = [
        model.get_tensor(d["index"])[0] for d in model.get_output_details()
    ]
    lines = []
    for (y, x), score in zip(centers, scores):
        ends = offsets[y, x].reshape(2, 2)
        if score > 0.10 and np.linalg.norm(ends[1] - ends[0]) > 20:
            lines.append((ends + [x, y]) * [image.shape[1] / 256, image.shape[0] / 256])
    return np.array(lines).reshape(-1, 4)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument(
        "--manifest", type=Path, help="Published issue #3 benchmark JSON"
    )
    inputs.add_argument("--roster", type=Path, help="Original issue panel images.csv")
    args = parser.parse_args()
    folder = args.directory
    folder.mkdir(parents=True, exist_ok=True)
    meta = download(folder / "metadata.jsonl", DATA + "/metadata.jsonl")
    if (
        hashlib.sha256(meta.read_bytes()).hexdigest()
        != "bad539fec59171893606bf883d80dcda6116a03a0a94d33ae87833387d203e95"
    ):
        raise ValueError("Wireframe metadata checksum mismatch")
    records = {
        Path(r["image_file"]).stem: r
        for r in map(json.loads, meta.read_text().splitlines())
    }
    ids = sorted(
        records.keys() - set(WITNESSES),
        key=lambda i: hashlib.sha256(("issue3-fixed-panel/" + i).encode()).digest(),
    )[:24]
    panel = {i: "sample" for i in ids}
    panel.update(dict.fromkeys(WITNESSES, "witness"))
    frozen = {}
    if args.manifest:
        photos = json.loads(args.manifest.read_text())["photos"]
        frozen = {Path(photo["image"]).stem: photo for photo in photos}
    elif args.roster:
        with args.roster.open() as roster:
            frozen = {Path(r["image"]).stem: r for r in csv.DictReader(roster)}
    for i, photo in frozen.items():
        if i not in records or photo["url"] != DATA + "/images/" + i + ".png":
            raise ValueError(f"Unknown image or different dataset revision for {i}")
        panel.setdefault(i, "issue")
    model_file = download(folder / "mlsd.tflite", MODEL)
    if (
        hashlib.sha256(model_file.read_bytes()).hexdigest()
        != "76b2ceae6b5e0000080e4846014dae1c0b653c0288cb70a852ccdc72e74c20b9"
    ):
        raise ValueError("M-LSD weight checksum mismatch")
    model = Interpreter(model_path=str(model_file), num_threads=1)
    model.allocate_tensors()

    def fetch(i):
        return download(folder / (i + ".png"), DATA + "/images/" + i + ".png")

    with ThreadPoolExecutor(max_workers=5) as pool:
        list(pool.map(fetch, panel))
    with (folder / "panel.csv").open("w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["id", "panel", "width", "height", "pixel_sha256"])
        for i, group in panel.items():
            img = ImageOps.exif_transpose(Image.open(folder / (i + ".png"))).convert(
                "RGB"
            )
            w, h = img.size
            assert [w, h] == records[i]["image_size"]
            size = tuple(round(d * 256 / max(w, h)) for d in (w, h))
            img = img.resize(size, Image.Resampling.LANCZOS)
            digest = hashlib.sha256(img.tobytes()).hexdigest()
            if i == WITNESSES[0]:
                assert (
                    digest
                    == "aebcda30a0f081e7f9a8c81cd24a50e0fb2fb8e4302eb494ea8c3fc03c8344e4"
                )
            if i in frozen and digest != frozen[i]["rgb_sha256"]:
                raise ValueError(f"Frozen RGB checksum mismatch for {i}")
            writer.writerow([i, group, *size, digest])
            img.save(folder / (i + ".small.png"))
            pixels = np.asarray(img)
            for suffix, model_input in [
                (".segments.csv", pixels),
                (".clahe.segments.csv", enhance_contrast(pixels)),
            ]:
                seg = run_model(model, model_input)
                np.savetxt(
                    folder / (i + suffix),
                    seg + 1,
                    delimiter=",",
                    header="x1,y1,x2,y2",
                    comments="",
                )
            # Map zero-based source pixel centers through the resize to one-based R coordinates.
            ref = (np.asarray(records[i]["lines"]) + 0.5) * (
                np.asarray(size) / [w, h]
            ) + 0.5
            np.savetxt(
                folder / (i + ".annot.csv"),
                ref.reshape(-1, 4),
                delimiter=",",
                header="x1,y1,x2,y2",
                comments="",
            )


if __name__ == "__main__":
    main()
