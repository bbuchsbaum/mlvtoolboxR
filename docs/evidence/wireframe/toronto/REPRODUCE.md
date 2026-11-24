# Reproducing the contour comparison

Use the Python environment, HAWPv2 checkout, CPU patch, and two model checkpoints from the [Wireframe setup](../TOPOLOGY-RUN.md). The commands below run from the package root. Set `work` to the directory containing `mlsd.tflite`, `hawp-src/`, `hawpv2-edb9b23f.pth`, and the Python environment. R needs the package's build dependencies, pkgload, png, and jsonlite.

## Toronto Scenes

Place the six authors’ `*_vecLD.mat` files in `reference-data/mlv/dataSets/TorontoScenes/`. Download `CP_images.zip` from the [authors’ archive](https://osf.io/9squn/). The preparation command verifies the archive and the six artist-vector checksums. Use a private staging directory for the photographs and artist vectors.

```sh
py="$work/venv/bin/python"
study=.agent-work/contour-study
stage=.agent-work/toronto/staged
out=docs/evidence/wireframe/toronto
"$py" tools/wireframe/toronto_scenes.py .agent-work/toronto/CP_images.zip --stage "$stage"
"$py" tools/wireframe/toronto_scenes.py .agent-work/toronto/CP_images.zip \
    --audit-inputs --compare-root "$work"
"$py" tools/wireframe/contour_inputs.py "$work" "$stage/panel.csv" "$stage" "$study/toronto"
"$py" tools/wireframe/contour_experiment.py "$stage" "$study/toronto" "$study/draw-toronto"
"$py" tools/wireframe/score_contour_drawings.py "$stage" "$study/draw-toronto" "$out"
```

Selection consumes resized RGB, the frozen context probabilities and full-image M-LSD segments. The R stage rasterizes the modified segments, runs `trace_contours(simplify = 1)`, and calls `detect_junctions()` on that drawing. It also traces two full-scene controls, full-image M-LSD segments with no joins added and the package's default FDoG mask. The artist ink reference uses the package rasterizer on the artist's vectors. Public tables contain photo IDs, hashes or aggregate measurements. The 325 reserve photos have no model scores.

## York gap replay

The scorer reads the frozen events and the reviewed gap ledger directly from the repository and verifies every context point, type, confidence and output order. Prepare the 51 photos and source segments using the official [York Urban archive](https://www.elderlab.yorku.ca/resources/york-urban-line-segment-database-information/).

```sh
york="$study/york-photos"
"$py" tools/wireframe/prepare_york.py "$work/YorkUrbanDB.zip" "$york" --split train --originals
"$py" tools/wireframe/contour_inputs.py "$work" "$york/panel.csv" "$york" "$study/york"
"$py" tools/wireframe/contour_experiment.py "$york" "$study/york" "$study/draw-york" --york
"$py" tools/wireframe/contour_york.py "$york" "$study/draw-york" "$out"
```

## Development drawing pilot

The following commands use the existing 60 development drawings. They preserve the selected contour IDs for the direct vector comparison and retrace identical segments for the raster comparison.

```sh
Rscript --vanilla tools/wireframe/contour_union.R "$stage" "$study/draw-toronto"
"$py" tools/wireframe/score_contour_drawings.py "$stage" "$study/draw-toronto" "$out" \
    --development-pilot --workers 4
```

For validation, the dense assignment algorithm in `score.match` independently reproduced all 720 point-count entries and 1,080 drawing-count entries on the 90 confirmation photos. A direct call to the single-photo context CLI reproduced the archived P1020177 CSV byte for byte. The comparison exercised the actual R package contour tracer, rasterizer and junction detector.
