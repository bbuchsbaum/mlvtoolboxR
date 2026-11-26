"""Score MLV-statistic agreement between automatic and artist drawings.

Frozen plan: docs/evidence/wireframe/mlv-agreement/PLAN.md.
Usage: python analyze.py WORK_DIR PUBLIC_DIR
WORK_DIR holds histograms.csv and scalars.csv from compute.R.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

SEED = 20261126
BOOT = 2000
METHODS = ["fdog", "mlsd", "union_vector", "curve_rank"]
CONFIGS = ["native", "working", "published"]
PROPS = {
    "native": [
        "orientation",
        "length",
        "curvature",
        "horizontal_vertical",
        "junction_type",
        "junction_angle",
    ],
    "working": [
        "orientation",
        "length",
        "curvature",
        "horizontal_vertical",
        "junction_type",
        "junction_angle",
    ],
    "published": [
        "orientation",
        "length",
        "curvature",
        "junction_type",
        "junction_angle",
    ],
}
JTYPES = ["T", "Y", "X", "Arrow", "Star"]
SCALARS = (
    ["total_length", "n_contours"]
    + [f"density_{t}" for t in JTYPES]
    + ["mean_abs_curvature", "hv_fraction"]
)
CATEGORIES = ["beaches", "cities", "forests", "highways", "mountains", "offices"]


def load(work, control=False):
    hist = defaultdict(list)
    meta = {}
    suffixes = [""]
    if control:
        suffixes += sorted(
            p.name[len("histograms") : -4] for p in work.glob("histograms-control*.csv")
        )
    for suffix in suffixes:
        with open(work / f"histograms{suffix}.csv") as f:
            for r in csv.DictReader(f):
                hist[(r["id"], r["method"], r["config"], r["property"])].append(
                    (int(r["bin"]), float(r["value"]))
                )
                meta[r["id"]] = (r["category"], r["set"])
    hist = {k: np.array([v for _, v in sorted(rows)]) for k, rows in hist.items()}
    scal = {}
    for suffix in suffixes:
        with open(work / f"scalars{suffix}.csv") as f:
            for r in csv.DictReader(f):
                scal[(r["id"], r["method"], r["config"], r["name"])] = (
                    float(r["value"]) if r["value"] not in ("NA", "") else np.nan
                )
    return hist, scal, meta


def js(p, q):
    sp, sq = p.sum(), q.sum()
    if sp == 0 and sq == 0:
        return np.nan
    if sp == 0 or sq == 0:
        return 1.0
    # Round so that distances equal up to floating-point noise (for example
    # 0.9999999999999999 and 1.0 for disjoint supports) tie when ranked.
    return round(float(jensenshannon(p / sp, q / sq, base=2)), 12)


def boot_indices(ids, meta, rng):
    by = defaultdict(list)
    for k, i in enumerate(ids):
        by[meta[i][0]].append(k)
    out = []
    for _ in range(BOOT):
        out.append(
            np.concatenate(
                [rng.choice(v, size=len(v), replace=True) for v in by.values()]
            )
        )
    return out


def interval(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    return (
        (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
        if len(v)
        else (np.nan, np.nan)
    )


def features(hist, scal, i, m, config):
    vec = []
    for p in ["orientation", "length", "curvature", "junction_angle"]:
        h = hist[(i, m, config, p)]
        vec.extend(h / h.sum() if h.sum() > 0 else np.zeros_like(h))
    sconf = "native" if config == "published" else config
    for t in JTYPES:
        if config == "published":
            jt = hist[(i, m, config, "junction_type")]
            total = scal[(i, m, sconf, "total_length")]
            d = jt[JTYPES.index(t)] / total * 1e4 if total > 0 else 0.0
        else:
            d = scal[(i, m, sconf, f"density_{t}")]
        vec.append(np.log1p(d if np.isfinite(d) else 0.0))
    vec.append(np.log1p(scal[(i, m, sconf, "total_length")]))
    return np.array(vec, float)


def loo_predictions(X_train_all, y, X_test_all):
    pred = np.empty(len(y), dtype=object)
    for k in range(len(y)):
        keep = np.arange(len(y)) != k
        sc = StandardScaler().fit(X_train_all[keep])
        clf = LogisticRegression(C=1.0, solver="lbfgs", max_iter=10000)
        clf.fit(sc.transform(X_train_all[keep]), y[keep])
        pred[k] = clf.predict(sc.transform(X_test_all[k : k + 1]))[0]
    return pred


def main():
    work, pub = Path(sys.argv[1]), Path(sys.argv[2])
    control = "--control" in sys.argv[3:]
    if control:
        # Post hoc representation control: score the retraced artist drawing
        # like a method. The primary frozen outputs are written without it.
        for p in sorted(work.glob("histograms-control*.csv")):
            res = p.name[len("histograms-control") : -4]
            METHODS.append("artist_retraced" + (f"_{res}" if res else ""))
    pub.mkdir(parents=True, exist_ok=True)
    hist, scal, meta = load(work, control)
    ids = sorted(meta)
    assert len(ids) == 180
    subsets = {"all": ids, "test": [i for i in ids if meta[i][1] == "test"]}
    rng = np.random.default_rng(SEED)
    boots = {name: boot_indices(sub, meta, rng) for name, sub in subsets.items()}

    # O1 distances and O2 identification, per photo.
    dist = {}  # (config, method, prop) -> array over ids
    pct = {}
    top1 = {}
    excluded = defaultdict(int)
    for config in CONFIGS:
        for p in PROPS[config]:
            A = [hist[(i, "artist", config, p)] for i in ids]
            for m in METHODS:
                M = [hist[(i, m, config, p)] for i in ids]
                D = np.array(
                    [[js(M[a], A[b]) for b in range(len(ids))] for a in range(len(ids))]
                )
                own = np.diag(D).copy()
                excluded[(config, m, p)] = int(np.isnan(own).sum())
                dist[(config, m, p)] = own
                pr, t1 = [], []
                for a in range(len(ids)):
                    others = np.delete(D[a], a)
                    others = others[np.isfinite(others)]
                    if not np.isfinite(own[a]) or not len(others):
                        pr.append(np.nan)
                        t1.append(np.nan)
                        continue
                    pr.append(
                        (np.sum(others > own[a]) + 0.5 * np.sum(others == own[a]))
                        / len(others)
                    )
                    t1.append(float(np.all(others > own[a])))
                pct[(config, m, p)] = np.array(pr)
                top1[(config, m, p)] = np.array(t1)
    # Reference: artist-to-artist distance between different photos of the same category.
    ref = {}
    for config in CONFIGS:
        for p in PROPS[config]:
            vals = []
            for a in range(len(ids)):
                for b in range(a + 1, len(ids)):
                    if meta[ids[a]][0] == meta[ids[b]][0]:
                        vals.append(
                            js(
                                hist[(ids[a], "artist", config, p)],
                                hist[(ids[b], "artist", config, p)],
                            )
                        )
            ref[(config, p)] = float(np.nanmedian(vals))

    rows = []
    for sub, sid in subsets.items():
        idx = np.array([ids.index(i) for i in sid])
        for config in CONFIGS:
            for p in PROPS[config] + ["mean_over_properties"]:
                for m in METHODS:
                    if p == "mean_over_properties":
                        d = np.nanmean(
                            np.vstack(
                                [dist[(config, m, q)][idx] for q in PROPS[config]]
                            ),
                            axis=0,
                        )
                        pc = np.nanmean(
                            np.vstack(
                                [pct[(config, m, q)][idx] for q in PROPS[config]]
                            ),
                            axis=0,
                        )
                        t1 = np.nanmean(
                            np.vstack(
                                [top1[(config, m, q)][idx] for q in PROPS[config]]
                            ),
                            axis=0,
                        )
                        excl = sum(excluded[(config, m, q)] for q in PROPS[config])
                        r = np.nan
                    else:
                        d, pc, t1 = (
                            dist[(config, m, p)][idx],
                            pct[(config, m, p)][idx],
                            top1[(config, m, p)][idx],
                        )
                        excl = int(np.isnan(d).sum())
                        r = ref[(config, p)]
                    rows.append(
                        dict(
                            subset=sub,
                            config=config,
                            property=p,
                            method=m,
                            photos=len(idx),
                            excluded_pairs=excl,
                            median_distance=float(np.nanmedian(d)),
                            mean_distance=float(np.nanmean(d)),
                            same_category_artist_median=r,
                            mean_identification_percentile=float(np.nanmean(pc)),
                            top1_rate=float(np.nanmean(t1)),
                        )
                    )
    write(pub / "distribution-agreement.csv", rows)

    # O3 scalar summaries.
    rows = []
    for sub, sid in subsets.items():
        for config in ["native", "working"]:
            for name in SCALARS:
                art = np.array([scal[(i, "artist", config, name)] for i in sid])
                for m in METHODS:
                    met = np.array([scal[(i, m, config, name)] for i in sid])
                    ok = np.isfinite(art) & np.isfinite(met)
                    rho = (
                        spearmanr(met[ok], art[ok]).statistic
                        if ok.sum() > 2
                        else np.nan
                    )
                    pos = ok & (art > 0)
                    ratio = (
                        float(np.median(met[pos] / art[pos])) if pos.any() else np.nan
                    )
                    rows.append(
                        dict(
                            subset=sub,
                            config=config,
                            scalar=name,
                            method=m,
                            photos=int(ok.sum()),
                            spearman=float(rho),
                            median_ratio=ratio,
                            artist_median=float(np.nanmedian(art)),
                            method_median=float(np.nanmedian(met)),
                        )
                    )
    write(pub / "scalar-agreement.csv", rows)

    # O4 category decoding, leave one photo out over all 180.
    y = np.array([meta[i][0] for i in ids])
    correct = {}
    for config in CONFIGS:
        XA = np.vstack([features(hist, scal, i, "artist", config) for i in ids])
        correct[(config, "artist", "artist")] = loo_predictions(XA, y, XA) == y
        for m in METHODS:
            XM = np.vstack([features(hist, scal, i, m, config) for i in ids])
            correct[(config, m, "method")] = loo_predictions(XM, y, XM) == y
            correct[(config, m, "transfer")] = loo_predictions(XA, y, XM) == y
    rows = []
    for sub, sid in subsets.items():
        idx = np.array([ids.index(i) for i in sid])
        for (config, m, kind), c in correct.items():
            cc = c[idx].astype(float)
            ci = interval([cc[b].mean() for b in boots[sub]])
            label = {
                "artist": "artist -> artist",
                "method": "method -> method",
                "transfer": "artist -> method",
            }[kind]
            rows.append(
                dict(
                    subset=sub,
                    config=config,
                    method=m,
                    comparison=label,
                    photos=len(idx),
                    accuracy=float(cc.mean()),
                    ci_low=ci[0],
                    ci_high=ci[1],
                )
            )
    write(pub / "category-decoding.csv", rows)

    # Paired bootstrap: selector minus FDoG and minus the earlier blend.
    rows = []
    for sub, sid in subsets.items():
        idx = np.array([ids.index(i) for i in sid])
        for config in CONFIGS:
            for comp in ["fdog", "union_vector"]:

                def add(outcome, a, b, stat=np.nanmean):
                    diff = stat(a) - stat(b)
                    ci = interval([stat(a[bb]) - stat(b[bb]) for bb in boots[sub]])
                    rows.append(
                        dict(
                            subset=sub,
                            config=config,
                            comparison=f"curve_rank - {comp}",
                            outcome=outcome,
                            difference=float(diff),
                            ci_low=ci[0],
                            ci_high=ci[1],
                        )
                    )

                for p in PROPS[config]:
                    add(
                        f"mean JS distance: {p}",
                        dist[(config, "curve_rank", p)][idx],
                        dist[(config, comp, p)][idx],
                    )
                md = lambda m: np.nanmean(
                    np.vstack([dist[(config, m, q)][idx] for q in PROPS[config]]),
                    axis=0,
                )
                add(
                    "mean JS distance: mean over properties", md("curve_rank"), md(comp)
                )
                mp = lambda m: np.nanmean(
                    np.vstack([pct[(config, m, q)][idx] for q in PROPS[config]]), axis=0
                )
                add(
                    "identification percentile: mean over properties",
                    mp("curve_rank"),
                    mp(comp),
                )
                for kind in ["method", "transfer"]:
                    add(
                        f"decoding accuracy: {kind}",
                        correct[(config, "curve_rank", kind)][idx].astype(float),
                        correct[(config, comp, kind)][idx].astype(float),
                    )
                if config != "published":
                    for name in SCALARS:
                        art = np.array([scal[(i, "artist", config, name)] for i in sid])
                        a = np.array(
                            [scal[(i, "curve_rank", config, name)] for i in sid]
                        )
                        b = np.array([scal[(i, comp, config, name)] for i in sid])

                        def rho(x, sel=None):
                            xs, ys = (x, art) if sel is None else (x[sel], art[sel])
                            ok = np.isfinite(xs) & np.isfinite(ys)
                            return (
                                spearmanr(xs[ok], ys[ok]).statistic
                                if ok.sum() > 2
                                else np.nan
                            )

                        diff = rho(a) - rho(b)
                        ci = interval([rho(a, bb) - rho(b, bb) for bb in boots[sub]])
                        rows.append(
                            dict(
                                subset=sub,
                                config=config,
                                comparison=f"curve_rank - {comp}",
                                outcome=f"spearman: {name}",
                                difference=float(diff),
                                ci_low=ci[0],
                                ci_high=ci[1],
                            )
                        )
    write(pub / "paired-bootstrap.csv", rows)

    # Per category O1 (mean distance) and O3 (median ratio), native frame, all photos.
    rows = []
    for cat in CATEGORIES:
        idx = np.array([k for k, i in enumerate(ids) if meta[i][0] == cat])
        for m in METHODS:
            row = dict(category=cat, method=m, photos=len(idx))
            for p in PROPS["native"]:
                row[f"js_{p}"] = float(np.nanmean(dist[("native", m, p)][idx]))
            for name in SCALARS:
                art = np.array([scal[(ids[k], "artist", "native", name)] for k in idx])
                met = np.array([scal[(ids[k], m, "native", name)] for k in idx])
                pos = np.isfinite(art) & np.isfinite(met) & (art > 0)
                row[f"ratio_{name}"] = (
                    float(np.median(met[pos] / art[pos])) if pos.any() else np.nan
                )
            rows.append(row)
    write(pub / "by-category.csv", rows)

    # Per-photo native distances (aggregates only, no geometry).
    rows = []
    for k, i in enumerate(ids):
        for m in METHODS:
            row = dict(id=i, category=meta[i][0], set=meta[i][1], method=m)
            for p in PROPS["native"]:
                row[f"js_{p}"] = dist[("native", m, p)][k]
                row[f"pct_{p}"] = pct[("native", m, p)][k]
            rows.append(row)
    write(pub / "photo-distances.csv", rows)
    print("written", pub)


def write(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow(
                {k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()}
            )


if __name__ == "__main__":
    main()
