# Where the photographic junction detector can improve

The most promising new scorer improved T recall by **2.72 percentage points** across the 342 previously evaluated Wireframe photographs. It lost **0.80 points** of X recall and flagged a gap on the original issue panel. Using its T scores while keeping the existing X locations and photo-level detection budgets gained just **15 T matches** and flagged one additional gap across the 342 photos. The shipped scorer remains the reference for both types. [Prior held-out evaluation](TOPOLOGY.md)

## Start with the locations

The table isolates each restriction in the existing scorer. Each entry counts the maximum possible one-to-one matches of that type within four pixels on the same 342 photographs. The raw and anchored rows give upper bounds with no limit on the number of eligible candidates.

| Candidates eligible to match | T matches / 6,292 | X matches / 1,243 |
| --- | ---: | ---: |
| Raw candidates | 6,129 | 1,130 |
| Candidate locations after anchoring | 6,027 | 1,105 |
| Locations clearing the class cutoff | 4,103 | 750 |
| Locations assigned that class before suppression | 3,750 | 686 |
| Final output | 3,540 | 649 |

The biggest reduction comes from scoring. On the original floor photo, 23 of the 33 annotated Xs could match a candidate before scoring, but only six could match candidates that passed the cutoff. The four-pixel suppression step alone accounts for 37 lost X matches across the 342 photos. [Scoring code](../../../tools/wireframe/topology.py)

## Controlled comparisons

All fitted variants used the same 320 training photographs and the same candidate positions. Training weights either favored better-centered positives, gave equal weight to each represented reference, or did both. Two feature variants measured matching edge polarity on opposing rays or measured through-support within individual finite segments. The neural geometry scorer used the existing features with two hidden layers. The image variant also saw a 25-pixel grayscale patch around each candidate, augmented by rotations and reflection. Both used the same hard gap negatives and eight training passes. Each comparison below retuned the score thresholds on the same 95 calibration photos. The search favored improving both types with at most 5,079 total detections and four reviewed gap flags. Each type could use up to 105% of its baseline count.

| Variant | T matches / 1,688 | X matches / 363 | Detections | Gap flags / 42 |
| --- | ---: | ---: | ---: | ---: |
| Shipped scorer | 1,002 | 215 | 5,079 | 4 |
| Retune thresholds only | 1,007 | 215 | 5,077 | 4 |
| Retune with 3 px suppression | 995 | 217 | 4,903 | 4 |
| Retune with 5 px suppression | 964 | 209 | 4,936 | 4 |
| Favor centered candidates | 1,036 | 216 | 5,072 | 4 |
| Equal weight per reference | 997 | 217 | 4,972 | 4 |
| Equal weight per reference, favor centered candidates | 1,022 | 218 | 4,954 | 4 |
| Signed gradient continuity | 985 | 213 | 5,075 | 3 |
| Finite-segment support | 992 | 216 | 4,770 | 3 |
| More and larger trees | 978 | 215 | 4,880 | 4 |
| Neural geometry scorer | **1,029** | **218** | 5,003 | 4 |
| Neural geometry and image patches | 997 | 214 | 5,070 | 3 |

The neural geometry scorer ranked first on the smaller recall gain, with mean gain breaking the tie. On the 342-photo retrospective check it recovered 3,711 of 6,292 Ts and 639 of 1,243 Xs, with 18,235 detections and five gap flags. The 95% paired photo-bootstrap interval for the T gain was **1.92 to 3.53 points**. For X it was **−3.18 to +1.33 points**. On the original panel it recovered 36/62 Ts and 58/101 Xs and flagged one of 15 reviewed gaps. The floor X count rose from six to nine. [Original benchmark](https://github.com/bbuchsbaum/mlvtoolboxR/issues/3#issuecomment-7101349061)

We then retained every shipped X location and reranked candidate Ts using the neural geometry scores. We capped Ts at the original count for each photo and chose one anchor-versus-new-site score offset on the 95 calibration photos. On the retrospective check it matched 3,555 Ts and the same 649 Xs at the same total detection count. The T gain was **0.24 points** (95% paired interval **−0.26 to +0.75**). Gap flags rose from five to six of 49. Both matching implementations agreed on the per-photo counts for the last two comparisons.

The promising change improved T ranking across the tested photos. The X experiments gave little gain under comparable detection and gap budgets. The prior masked review also found that 23 of 36 missed reference Xs looked like visible Ts. A test of better visual crossing recovery would benefit from explicit labels for visible crossings, occluded intersections and genuine gaps on fresh photos. This round reused all evaluation photos and adds no independent estimate of human-rated precision. [Earlier reference audit](prospective-reference-audit.csv)
