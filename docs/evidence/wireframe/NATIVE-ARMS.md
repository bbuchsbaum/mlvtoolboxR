# Verifying junction arms at original resolution

The original-resolution arm verifier failed the fresh York comparison. The existing context scorer is the better-supported choice. [Results](#results)

## Comparison fixed before evaluation

Start with the frozen context detector at commit 5088d3a. Preserve every predicted location and the per-photo T and X counts. Train three otherwise identical binary tree models on the 23 development, 24 replication and 48 confirmation Wireframe photos. Use unambiguous one-to-one matches to reference Ts or Xs as training labels. The three inputs are the frozen type and confidence alone, those inputs plus arm measurements from the 256-pixel working image, and those inputs plus the same arm measurements from the original image. The third is the proposed change. Select no model or cutoff on the evaluation photos.

Trace 36 oriented rays per site, in the same physical field of view for both resolutions. On each ray, choose one transverse offset and width, require consistent contrast polarity, and estimate how far the arm runs before two consecutive weak samples. Fit two axes at least 30 degrees apart and keep the measurements of all four arms, including the weak side of a possible T. Within each image, pair the highest-ranked existing Ts with the lowest-ranked existing Xs when the X score favors the swap. Allow at most one fifth of the original X predictions to change. All methods share exactly the same positions, so the flag count for a reviewed gap near *either* type is invariant. X-specific gap flags can change.

Use the 342 previously evaluated Wireframe photos as a retrospective check. Reserve the 51 photos from York Urban’s official training half for the new external comparison. Neither the context scorer nor the new verifiers trained on York. The earlier York experiment used the dataset’s separate official test half. Match by type within four working pixels and resample whole photos in pairs for uncertainty. For promotion, require more annotated X matches than both the frozen and working-resolution comparators on fresh York, no more than a one-percentage-point T recall loss against the frozen scorer, and a visible-X gain in blinded review. Use independent human confirmation before recommending a new default.

Sample shuffled native-pixel review cards from new York photos, stratified by unchanged Xs, promoted Ts, demoted Xs and unchanged Ts, with separate reference and geometric near-miss controls. Keep sampling weights and type identities hidden during review. Independently sample unmarked regions for two people to mark all visible Ts and Xs, enabling recall measurement outside the detector's proposals. Distinguish human ratings from a two-Astra pilot. Record all comparisons, including a negative outcome.

## Results

**Keep the existing context scorer.** On the 51 newly evaluated York photos, original-resolution arms lost four annotated X matches and eleven T matches. The working-resolution version lost one X and eight Ts. On the reused Wireframe photos, original resolution gained six X matches and lost 65 Ts. The original pixels occasionally changed the right label, but those gains did not carry over to the fresh comparison. [Per-photo counts](native-arms/scores.csv)

| Scorer | Wireframe T / 6,292 | Wireframe X / 1,243 | York T / 1,975 | York X / 1,297 | T/X pairs exchanged, Wireframe / York |
| --- | ---: | ---: | ---: | ---: | ---: |
| Frozen context | 3,781 | 658 | 683 | 373 | 0 / 0 |
| Type and confidence alone | 3,781 | 658 | 683 | 373 | 0 / 0 |
| Working-resolution arms | 3,747 | 660 | 675 | 372 | 222 / 72 |
| Original-resolution arms | 3,716 | 664 | 672 | 369 | 329 / 109 |

Each scorer emitted 12,444 Ts and 6,064 Xs on Wireframe and 3,441 Ts and 2,523 Xs on York. Against the frozen context scorer, the native-arm change on York was −0.56 percentage points of T recall (95% paired photo-bootstrap interval −1.09 to −0.06) and −0.31 points of X recall (−0.96 to 0.00). Against working-resolution arms it was −0.15 points for T (−0.75 to +0.46) and −0.23 for X (−0.83 to 0.00). On Wireframe the corresponding changes against context were −1.03 for T (−1.38 to −0.70) and +0.48 for X (−0.08 to +1.09). Each interval uses 10,000 paired resamples by photograph, seed 71323. [Frozen photo and pixel roster](native-arms/photos.csv) · [York detections](native-arms/york-events.csv.gz)

**Visible crossings.** Two independent, blinded Astra sessions rated 108 cards and agreed on the exact class on 83. For 20 sampled promotions to X out of 109, reviewer A saw six Xs and reviewer B saw two. For 20 of the 109 demotions from X, each reviewer saw three. Weighting the four sampling groups by their population, estimated X precision changed from 38.9% to 39.6% for A and from 34.1% to 33.9% for B. They both saw visible Xs in six of the eight reference X controls. The reviewers used the same model family. These are pilot estimates. [Sampling key](native-arms/cards.csv) · [Individual ratings](native-arms/astra-ratings.csv)

**Visible recall on unmarked regions.** As a supplementary pilot, two additional blinded Astra sessions marked every visible T and X in the first eight of the 24 randomly sampled regions. Reviewer A marked 19 Xs and reviewer B marked 30. The frozen and native methods each recovered 8/19 against A and 11/30 against B. Native T matches were 86/188 versus 87/188 for the frozen scorer against A, and both recovered 78/158 against B. Neither reviewer marked a clear gap in those eight regions. The differing reference counts illustrate the difficulty of judging dense architectural patterns. [Independent region marks](native-arms/astra-region-marks.csv)

**False bridges.** Neither reviewer called any of the 12 randomly interspersed geometrically proposed gaps a gap on the ordinary target cards. We therefore ran a supplementary blinded screen using the earlier protocol with a separate overlay identifying the intended stem and target. The same reviewers accepted 29 and 25 of 45 candidates, with 24 accepted by both. Each method flagged ten of those 24 with either type, including one with an X. On the earlier 49 reviewed Wireframe gaps each method flagged five with either type and two with an X. Unchanged all-type gap flags follow from preserving every position. [Fresh gap decisions and individual flags](native-arms/gap-review.csv)

**Interpretation.** The training set supplied 1,165 unambiguous T locations and 246 Xs matched to line-derived references. This training target rewards geometric intersections even when an edge disappears behind another surface. A different training target using independently marked visible arms is still plausible. Native resolution on its own did not produce the proposed improvement with this fitter, angular sampling and training set. Near the image boundary the sampler uses edge-pixel replication. The protected human review consists of 24 randomly placed unmarked regions and 108 masked cards. Human-rated precision and recall are pending.

## Reproduce the comparison and human packet

Use the existing Wireframe environment and pinned photographs. The evidence directory contains the exact frozen detector outputs and reference coordinates used for fitting and evaluation. The three fitted experimental models are stored with the evidence. The single-photo inference command continues to use the published scorers. [Setup](TOPOLOGY-RUN.md) · [Archived inputs](native-arms/inputs.csv.gz) · [Experimental models](native-arms/models.json.gz) · [York Urban data, Denis, Elder and Estrada, 2008](https://www.elderlab.yorku.ca/resources/york-urban-line-segment-database-information/)

From the repository root, the following stages the 51 new York photographs, applies the frozen models, and rebuilds the private human packet. Use a new destination for each staging command. The packet contains only photos, target identifiers and rating instructions. Keep the sampling keys separate during review.

```sh
py="$work/venv/bin/python"
e=docs/evidence/wireframe/native-arms
"$py" tools/wireframe/prepare_york.py "$work/YorkUrbanDB.zip" "$work/york-originals" --split train --originals
"$py" tools/wireframe/arm_experiment.py stage fresh "$work/york-originals" "$work/arm-york"
"$py" tools/wireframe/arm_experiment.py score "$work/arm-york" "$e/models.json.gz" "$work/york-comparison"
"$py" tools/wireframe/arm_review.py "$work/arm-york" "$work/york-comparison.sites.csv" "$work/packet" "$work/keys"
```

Each human reviewer opens their own copy of `packet/index.html`, completes the unmarked regions first, rates the cards, then exports a JSON file. Score each review independently. The scoring command checks for missing cards and unfinished regions. It estimates card precision using the sampling weights, and measures visible recall by one-to-one matching against the independent region marks. [Card key](native-arms/cards.csv) · [Region coordinates](native-arms/regions.csv)

```sh
"$py" tools/wireframe/score_arm_review.py "$e" "$e/york-events.csv.gz" reviewer-a.json
"$py" tools/wireframe/score_arm_review.py "$e" "$e/york-events.csv.gz" reviewer-b.json
```

To reproduce the fit, use the prior run's original and resized Wireframe photo directories. The stage command verifies the pixel hashes and recreates the exact 95-photo fit inputs. The same command accepts `prospective` and `remainder` for the retrospective scores.

```sh
for split in development replication confirmation; do
  source_photos="$work/$split"
  if [ "$split" = development ]; then source_photos="$work"; fi
  "$py" tools/wireframe/arm_experiment.py stage "$split" "$source_photos" "$work/arm-$split"
done
"$py" tools/wireframe/arm_experiment.py fit "$work/refitted.gz" \
  "$work/arm-development" "$work/arm-replication" "$work/arm-confirmation"
```

A fresh fit reproduced the frozen model file byte-for-byte. Replaying the new York set from the archived inputs reproduced the per-photo results and individual decisions byte-for-byte. Pixel controls tested two image scales, complete three- and four-arm patterns, and a six-pixel gap. The offline human form was exercised in a temporary browser profile through opening an image, marking a location, saving, reloading, undoing and navigating to the final card.
