# Toronto development crossing opportunity check

This exploratory, reference-guided analysis uses the 60 development photos at repository commit `22609d5`. Each photo is evaluated at its 256-pixel working resolution. Neither probe predicts junctions from a new photograph. Both ask whether existing source branches can trace four directions when the artist crossing location and directions are already known.

Of 1,207 artist X labels, 729 have four finite, distinct exported directions. In the traced-branch probe, 84 of those allow four different branches of the named artist polylines, each at least nine pixels long, to be assigned to the four rays. The artist branches are sampled at arc lengths 3, 5, 7 and 9 from the point on each path nearest the crossing. A source branch must approach the centre within four pixels, pass within two pixels of each checkpoint in sequence, and agree with the artist's direction from first to last checkpoint within 25 degrees. Each outward branch can explain at most one arm. Distinct fragments from the two sources can overlap, and this screen does not establish physical connectivity or enforce an ink budget.

All 84 eligible references match themselves by construction. Five have four supported arms in the full combined FDoG and M-LSD path pools (all on one photo), zero in the previous selector’s final drawings. One of the five lacks any source X declaration within four pixels. Rotating the requested branches by 45 degrees around the same sites gives zero matches. Of the three strict source X matches from the previous experiment, two have measurable templates and both pass. The third is excluded by the branch definition. The 79 misses among 84 are specific to this screen and its length and continuity requirements. The other 645 four-ray references are unclassified by this test. A separate straight-ray probe gives five combined opportunities among 118 artist-supported cases. Its five are across two photos and need not identify the same sites.

The next comparison should record recoverable reference counts as a function of arm length, then test the same extractor at two input resolutions in common image coordinates. Count newly available arms at fixed references, near-miss errors, and retained contour/T performance before using fresh photos. Compare method results on common eligible references and report denominators.

## Reproduction

Use the same private input and Python environment as the [crossing-group comparison](../crossing-groups/REPRODUCE.md). Set `baseline` to the previously rendered held-out drawings from the single-path selector. Run the following commands from the package root. They print the aggregate counts and write per-reference tables to the specified private output directory.

```sh
py="$work/venv/bin/python"
staged=.agent-work/toronto/staged
base=.agent-work/contour-study/draw-toronto
baseline=.agent-work/crossing-groups/reproduce/baseline
out=.agent-work/arm-opportunity-replay
"$py" tools/wireframe/arm_opportunity.py "$staged" "$base" "$baseline" "$out"
"$py" tools/wireframe/traced_arm_opportunity.py "$staged" "$base" "$baseline" "$out"
```
