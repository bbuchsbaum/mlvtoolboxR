"""MLV agreement of the photo-only M3 source drawings on the 60 development photos.

Exploratory and outside the frozen screening rule. Reuses the measures of
tools/wireframe/mlv_agreement/analyze.py in the original frame, and checks the
M3 pixel coverage independently from the rendered PNGs.

Usage: python m3_agreement.py AGREEMENT_WORK M3_STATS_DIR DRAW_TORONTO OUT_CSV
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mlv_agreement"))
from analyze import JTYPES, features, js  # noqa: E402

work, stats_dir, draw, out_csv = map(Path, sys.argv[1:5])
PROPS = [
    "orientation",
    "length",
    "curvature",
    "horizontal_vertical",
    "junction_type",
    "junction_angle",
]
METHODS = ["fdog", "mlsd", "union_vector", "curve_rank", "fdog512_m3", "sed_m3"]

hist = defaultdict(list)
scal = {}
meta = {}
for hfile, sfile in (
    (work / "histograms.csv", work / "scalars.csv"),
    (stats_dir / "histograms-sources.csv", stats_dir / "scalars-sources.csv"),
):
    for r in csv.DictReader(open(hfile)):
        if r["config"] != "native":
            continue
        hist[(r["id"], r["method"], "native", r["property"])].append(
            (int(r["bin"]), float(r["value"]))
        )
        meta[r["id"]] = (r["category"], r["set"])
    for r in csv.DictReader(open(sfile)):
        if r["config"] == "native":
            scal[(r["id"], r["method"], "native", r["name"])] = (
                float(r["value"]) if r["value"] not in ("NA", "") else np.nan
            )
hist = {k: np.array([v for _, v in sorted(rows)]) for k, rows in hist.items()}
dev = sorted(i for i in meta if meta[i][1] == "development")
assert len(meta) == 180 and len(dev) == 60
# The screen's plan excludes the test photos, so the classifier's training rows
# and the identification candidates are the 60 development artist drawings.
ids = dev


# Independent pixel check of M3 drawings (coverage and precision within 2 px).
def mask(path):
    return np.asarray(Image.open(path).convert("L")) > 0


pix = defaultdict(lambda: np.zeros(4))
for i in dev:
    art = mask(draw / f"{i}.artist.drawing.png")
    near_art = distance_transform_edt(~art) <= 2
    for m in ("fdog512", "sed"):
        pred = mask(stats_dir / f"{i}.{m}.m3.png")
        near_pred = distance_transform_edt(~pred) <= 2
        pix[m] += [
            (art & near_pred).sum(),
            art.sum(),
            (pred & near_art).sum(),
            pred.sum(),
        ]

# Artist-trained classifier for every development photo (train on the other 179 artist drawings).
y = np.array([meta[i][0] for i in ids])
XA = np.vstack([features(hist, scal, i, "artist", "native") for i in ids])
transfer = {m: [] for m in METHODS}
for i in dev:
    k = ids.index(i)
    keep = np.arange(len(ids)) != k
    sc = StandardScaler().fit(XA[keep])
    clf = LogisticRegression(C=1.0, solver="lbfgs", max_iter=10000).fit(
        sc.transform(XA[keep]), y[keep]
    )
    for m in METHODS:
        x = features(hist, scal, i, m, "native")
        transfer[m].append(clf.predict(sc.transform(x[None, :]))[0] == y[k])

rows = []
for m in METHODS:
    row = dict(method=m, photos=len(dev), transfer_accuracy=float(np.mean(transfer[m])))
    for p in PROPS:
        A = {j: hist[(j, "artist", "native", p)] for j in ids}
        d_own, pct = [], []
        for i in dev:
            h = hist[(i, m, "native", p)]
            d = {j: js(h, A[j]) for j in ids}
            own = d[i]
            others = np.array([v for j, v in d.items() if j != i and np.isfinite(v)])
            d_own.append(own)
            if np.isfinite(own) and len(others):
                pct.append(
                    (np.sum(others > own) + 0.5 * np.sum(others == own)) / len(others)
                )
        row[f"js_{p}"] = float(np.nanmean(d_own))
        row[f"pct_{p}"] = float(np.mean(pct))
    row["js_mean"] = float(np.mean([row[f"js_{p}"] for p in PROPS]))
    for name in [
        "total_length",
        "n_contours",
        "density_T",
        "density_X",
        "mean_abs_curvature",
    ]:
        art = np.array([scal[(i, "artist", "native", name)] for i in dev])
        met = np.array([scal[(i, m, "native", name)] for i in dev])
        ok = np.isfinite(art) & np.isfinite(met)
        row[f"rho_{name}"] = float(spearmanr(met[ok], art[ok]).statistic)
        pos = ok & (art > 0)
        row[f"ratio_{name}"] = float(np.median(met[pos] / art[pos]))
    if m in ("fdog512_m3", "sed_m3"):
        c, a, pc, p = pix[m.replace("_m3", "")]
        cov, prec = c / a, pc / p
        row.update(
            pixel_coverage=cov,
            pixel_precision=prec,
            pixel_f1=2 * cov * prec / (cov + prec),
            ink=int(p),
        )
    rows.append(row)
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=sorted(
            {k for r in rows for k in r}, key=lambda k: (k != "method", k)
        ),
    )
    w.writeheader()
    for r in rows:
        w.writerow(
            {k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()}
        )
for r in rows:
    print(
        r["method"],
        {
            k: (round(v, 3) if isinstance(v, float) else v)
            for k, v in r.items()
            if k.startswith(
                (
                    "transfer",
                    "js_mean",
                    "pixel",
                    "ratio_density_T",
                    "rho_density_T",
                    "ratio_n_contours",
                    "ratio_mean_abs",
                    "pct_junction_type",
                    "pct_curvature",
                    "js_length",
                    "pct_length",
                )
            )
        },
    )
