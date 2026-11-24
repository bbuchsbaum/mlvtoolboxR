# Reproducing the fixed-ink contour comparison

Run these commands from the package root. Follow the [Toronto preparation instructions](../toronto/REPRODUCE.md) through the development drawing pilot to produce the private development assets and original vector blend. Use the same pinned RGB inputs and M-LSD weights. The recorded run used R 4.5.1, Python 3.12.11, scikit-learn 1.7.2, NumPy 2.5.5, SciPy 1.18.1, Pillow 12.1.1, and OpenCV 5.1.0. The public model file contains numeric splits and leaf values. The commands write photos and artist geometry only to the private trial directory.

## Development

```sh
py="$work/venv/bin/python"
trial=.agent-work/curve-rank
stage=.agent-work/toronto/staged
base=.agent-work/contour-study/draw-toronto
pub=docs/evidence/wireframe/curve-rank
gate="$pub/development/curve-rank-gate.json"
model="$pub/model.json"
Rscript --vanilla tools/wireframe/export_curve_pixels.R "$base" "$pub/development.csv" "$trial/pixels-dev"
"$py" tools/wireframe/curve_rank.py crossfit "$stage" "$base" "$trial/pixels-dev" "$trial/development"
Rscript --vanilla tools/wireframe/contour_drawings.R "$stage" "$trial/development" "$trial/development/r-plan.csv"
"$py" tools/wireframe/score_curve_rank.py "$stage" "$base" "$trial/development" "$pub/development.csv" "$trial/report-dev"
```

`selection.csv` records each drawing's exact ink, number of selected paths, and SHA-256 of the full grouped vector export. Compare those values with the [recorded 60 selections](development/selection.csv). A fresh refit checks its serialized predictor against the scikit-learn fit on every development training row. To reproduce the freeze, use the locally recomputed development gate and compare the model bytes to `model.json`.

```sh
"$py" tools/wireframe/curve_rank.py refit "$stage" "$base" "$trial/pixels-dev" "$trial/development" \
  --gate "$trial/report-dev/curve-rank-gate.json" --model "$trial/refit.json"
cmp "$trial/refit.json" "$model"
```

## Reserve test

The fixed test roster lists 20 photos in each of six categories. The prediction stage needs the photographs and photo-derived source drawings. It can run before the artist geometry is staged. The earlier blend uses the same selection order and drawing code as the development pilot.

```sh
archive=.agent-work/toronto/CP_images.zip
weights="$work/mlsd.tflite"
res="$trial/reserve"
"$py" tools/wireframe/prepare_curve_reserve.py photos "$archive" "$weights" "$res"
Rscript --vanilla tools/wireframe/export_curve_pixels.R "$res/base" "$pub/reserve-test.csv" "$res/pixels"
"$py" tools/wireframe/curve_rank.py predict "$res/staged" "$res/base" "$res/pixels" "$res/predicted" \
  --gate "$gate" --model "$model"
Rscript --vanilla tools/wireframe/contour_drawings.R "$res/staged" "$res/predicted" "$res/predicted/r-plan.csv"
Rscript --vanilla tools/wireframe/contour_union.R "$res/staged" "$res/base" "$pub/reserve-test.csv"
"$py" tools/wireframe/prepare_curve_reserve.py references "$archive" "$weights" "$res"
"$py" tools/wireframe/score_curve_rank.py "$res/staged" "$res/base" "$res/predicted" \
  "$pub/reserve-test.csv" "$res/report"
```

The score table counts each artist junction and each predicted junction at most once per photo and class. Artist pixel coverage counts traced pixels reached within two pixels of the prediction. Pixel precision uses the reverse comparison, and F1 is the harmonic mean of pooled precision and coverage. The stricter incidence count also checks that the connecting contours agree around the site. Paired intervals resample photographs within the six scene categories 2,000 times at seed 20261124.
