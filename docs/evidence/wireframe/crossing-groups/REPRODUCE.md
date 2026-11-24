# Reproduce the development experiment

Use the [whole-contour study environment and private inputs](../curve-rank/REPRODUCE.md). The commands below require its development staging directory, FDoG and M-LSD source drawings, and their per-curve pixel exports. The 60 photo IDs, RGB hashes and fold assignments come from the original [development roster](../curve-rank/development.csv).

```sh
py="$work/venv/bin/python"
staged=.agent-work/toronto/staged
base=.agent-work/contour-study/draw-toronto
pixels=.agent-work/curve-rank/pixels-dev
trial=.agent-work/crossing-groups/reproduce
"$py" tools/wireframe/crossing_groups.py "$staged" "$base" "$pixels" "$trial"
for method in baseline grouped diagnostic; do
  Rscript --vanilla tools/wireframe/contour_drawings.R "$staged" "$trial/$method" "$trial/$method/r-plan.csv"
done
"$py" tools/wireframe/score_crossing_groups.py "$staged" "$base" "$trial/baseline" \
  "$trial/grouped" "$trial/diagnostic" "$trial/report"
```

The selector checks the SHA-256 of every baseline vector export against the earlier five-fold run. `selections.csv` contains the vector hashes and the number of complete source crossings for all three selection policies. For the reference-assisted diagnostic, `joint_actions` counts the reference groups admitted before filling the rest of the ink budget. The scorer checks each baseline's artist and detector counts and pixel metrics against the earlier report. Its paired intervals use 2,000 photo resamples stratified by scene category, with seed 20261124.
