# MLV port assessment: executed evidence

Date: 2026-10-25. Inputs were synthetic. No R port or performance
benchmark exists yet.

## Source revisions and scope

| Reference | Observed local checkout | Git HEAD |
|----|----|----|
| MATLAB MLV toolbox | `MLV_toolbox/` | `f12418a5e03c2f38e222fed1dffccf42b227cd29` |
| Python MLV port | `vendor/pyMLV/` | `14b3183b91146b1021eb7ef97db75ad703b06832` |

These are the paths observed during the assessment. The user is
arranging the vendoring layout; refresh the paths before implementation.
Both tracked-file status checks were empty after the probes. The probes
did not modify either checkout. Git emitted sandbox temporary-cache
warnings but returned status 0. No restricted image datasets were read.

## Octave availability

Octave 11.1.0 and image 2.20.0 were found locally. After
`pkg load image`:

- Available: `bwdist`, `bwmorph`, `bwareaopen`, `rgb2gray`,
  `imbinarize`.
- Missing from this installed environment: `smoothdata`, `knnsearch`,
  `insertShape`, `xmlread`, `imsegsam`.
- A 5-by-7 single-feature `bwdist` probe returned single-precision
  distances and uint32 nearest-feature indices. The feature distance was
  0 and a diagonal neighbor distance was approximately sqrt(2).
- Installed `bwmorph.m`, lines 191–197, explicitly documents that its
  `skel` result differs from MATLAB on the documented circles example.
  Therefore Octave execution is valuable evidence, with backend-specific
  qualification for morphology and tie-sensitive downstream
  calculations.

The initial sandboxed help-rendering command ended in a temporary-file
error after the capability and distance checks completed. A later
adapted MAT attempt using the `octave` GUI wrapper failed in Qt startup
(missing cocoa platform plugin, child exit -6); `octave-cli` avoided
that failure.

## Unmodified MATLAB source run under Octave

The following synthetic setup was passed through upstream functions:

``` matlab
pkg load image;
addpath("MLV_toolbox/MLVcode");
v = struct();
v.imsize = [30,20];
v.numContours = 1;
v.contours = {[1,1,11,1]};
v.originalImage = "synthetic";
v.lineMethod = "synthetic";
v = computeContourProperties(v, {"orientation","length","curvature"});
[v,h,b] = getHorizontalVerticalStats(v);
v.junctions = [];
[v,jh,jb] = getJunctionStats(v);

w = v;
w.imsize = [30,30];
w.numContours = 3;
w.contours = {[1,1,5,1], [1,10,5,10], [1,20,3,20]};
w = computeContourProperties(w);
[top,bottom] = splitLDbyProperties(w, "length", 0.5);

m = ones(25,35);
m(4:22,8) = 0;
m(4:22,27) = 0;
% The original computeMAT call was wrapped in try/catch to retain its error.
mat = computeMAT(m,28);
```

Observed transcript (process exit 0 because the MAT exception was
caught):

``` text
OMP: Warning #179: Function Can't set size of /tmp file failed:
core: length=10 orientation=0 curvature=0
horizontal histogram: total=0 expected_length=10 bins=8
empty junction angles: histogram_length=5 centers_length=8
split: source_length=10 top_length=8 bottom_length=2
MAT_BLOCKED: syntax error near line 67, column 67 in file /Users/bbuchsbaum/code/mlvtoolbox/MLV_toolbox/MLVcode/computeGradientVectorField.m
```

Interpretation:

- The analytical contour case succeeds: length 10, orientation 0,
  curvature 0.
- The horizontal/vertical histogram loses the full segment length at its
  upper boundary.
- Empty junction-angle values and bin centers have inconsistent lengths.
- Sorting the three contours by length and splitting at 0.5 retains 80%
  in the top group, exceeding the nominal 50% budget. The earlier
  arithmetic-only illustration did not apply this sort; this
  function-level result supersedes it.
- The original distance-field source contains two `jj++m_Neighbors8`
  expressions. Octave rejects this syntax before evaluating the
  function.

## Minimal Octave adaptation

A temporary copy of `computeGradientVectorField.m` replaced only the two
`jj++m_Neighbors8` expressions with `jj + m_Neighbors8`. In MATLAB, the
former expression denotes addition with unary plus; the explicit
addition avoids Octave’s parsing ambiguity. The original vendor file was
unchanged. The temporary directory was removed by Python’s
`TemporaryDirectory` context.

Reproduction driver:

``` python
from pathlib import Path
import tempfile, subprocess
src=Path('MLV_toolbox/MLVcode/computeGradientVectorField.m')
text=src.read_text()
assert text.count('jj++m_Neighbors8')==2
with tempfile.TemporaryDirectory(prefix='mlv-octave-assessment-') as td:
    (Path(td)/src.name).write_text(text.replace('jj++m_Neighbors8','jj + m_Neighbors8'))
    expr='pkg load image; addpath("'+str(src.parent.resolve())+'"); addpath("'+td+'","-begin"); fprintf("runtime=%s image=2.20.0 adaptation=two jj++ to jj + occurrences\\n",version); m=ones(25,35); m(4:22,8)=0; m(4:22,27)=0; try, mat=computeMAT(m,28); fprintf("adapted MAT: skeleton_pixels=%d finite_dist=%d finite_AOF=%d shape=%dx%d\\n",nnz(mat.skeleton),all(isfinite(mat.distance_map(:))),all(isfinite(mat.AOF(:))),size(mat.skeleton,1),size(mat.skeleton,2)); catch err, fprintf("ADAPTED_MAT_BLOCKED: %s\\n",err.message); exit(2); end;'
    proc=subprocess.run(['octave-cli','--quiet','--no-history','--eval',expr],capture_output=True,text=True)
    print(proc.stdout,end=''); print(proc.stderr,end='')
    print('process_exit='+str(proc.returncode))
    raise SystemExit(proc.returncode)
```

Observed transcript:

``` text
runtime=11.1.0 image=2.20.0 adaptation=two jj++ to jj + occurrences
adapted MAT: skeleton_pixels=44 finite_dist=1 finite_AOF=1 shape=25x35
process_exit=0
```

This demonstrates that the AOF/MAT construction path executes on one
small synthetic image with a documented syntax adaptation. It does not
validate the skeleton against MATLAB, exercise branch scores/projection,
or establish performance.

## Python assessment

Python 3.14.7 and NumPy 2.4.3 were used for the probe. SciPy and tqdm
were discoverable; cv2 and skimage were not installed in this Python
environment. No dependencies were installed and the full Python pipeline
was not run.

Static tree/module inspection found 60 top-level Python modules, 52 with
same-named MATLAB core counterparts, plus eight Python helpers. No
`pyproject.toml`, `setup.py`, `setup.cfg`, automated test directory or
GitHub workflow was identified in the pinned tree. Name overlap does not
establish functional coverage.

Missing local imports include `smoothData`, `diff`, `fitLineSegments`,
and `GetConSeg`. The Python gradient-field function loads
`outerBoundaryOriginal.mat` and opens diagnostic plots. These are
source-level findings, not a failed full-pipeline run.

An actual call to `computeOrientation` with this subpixel segment:

``` python
contours = np.empty((1, 1), dtype=object)
contours[0, 0] = np.array([[1.1, 1.1, 1.9, 1.9]], dtype=float)
v = {"numContours": np.array([[1]]), "contours": contours}
observed = computeOrientation(v)["orientations"][0][0]
expected = np.degrees(np.arctan2(-0.8, 0.8)) % 360
```

returned 0 degrees; the analytical result is 315 degrees. The
implementation casts all endpoints to int32 before subtraction. It
therefore discards the subpixel geometry in this example.

Probe receipt:

``` json
{"python": "3.14.7", "numpy": "2.4.3", "core_py_files": 60, "missing_local_imports": [{"file": "getBranchDerivative.py", "module": "smoothData"}, {"file": "getBranchDerivative.py", "module": "diff"}, {"file": "computeMATpropertyPerBranch.py", "module": "fitLineSegments"}, {"file": "computeMATpropertyPerBranch.py", "module": "smoothData"}, {"file": "computeMATpropertyPerBranch.py", "module": "diff"}, {"file": "traceSkeleton.py", "module": "GetConSeg"}], "subpixel_orientation": {"input": [[1.1, 1.1, 1.9, 1.9]], "observed": 0.0, "expected": 315.0}}
```

Further source-level risks, not executed here:

- `extract2DSkeletonFromBinaryImage` aliases `skeletonImage = fluxImage`
  and then thresholds in place, overwriting the flux values returned as
  AOF.
- It overrides the calculated component threshold with the constant 100.
- `computeMATpropertyPerBranch` changes MATLAB local derivative
  magnitudes to cumulative sums, and its taper branch references
  `newddR` before assignment.
- Several functions depend on MATLAB-file nesting conventions and mix
  coordinate/index representations.

These examples make pyMLV useful for ideas, adapters, and selected
differential tests, while each numerical path needs its own correctness
evidence.

## Boundaries

No package tests, cross-language end-to-end parity checks, speedups,
memory benchmarks, scientific method improvements, hosted CI or
distribution checks are claimed. The demonstrated differences should
become regression tests in the R implementation. Passing the adapted MAT
smoke is evidence of executability on this input only.
