# Frozen estimator acceptance contract, version 1

This contract precedes new-estimator implementation. The executable candidate
gate is `expect_curvature_contract(profile, turns)` in
`tests/testthat/helper-invariance.R`. Existing `mlv2023` fixtures remain separate.
New numerical methods must pass the gate; an absent method is not a skipped pass.

Contour geometry uses directed image coordinates, converted to Cartesian
orientation for signed turns. Only forward collinear runs are merged, including
the cyclic seam. Exactly closed contours wrap once. Exact 180-degree retracing
has no intrinsic signed turn and is rejected by new estimators. Open endpoints
carry no artificial turn.

`turning` distributes each vertex's signed turn uniformly over its dual
arc-length cell, bounded by neighboring segment midpoints. Open endpoint
half-segments have zero density. It supplies polygonal turning density and exact
signed/absolute vertex totals; it is not promised to resist geometric noise.

`arclength` convolves signed turning atoms with a triangular (tent) kernel of
physical half-width h. Each atom's kernel is normalized over its actual support
after clipping at open ends. Closed contours use circular distance and effective
h = min(h, L/2). The resulting field is piecewise linear. Partition it at all
kernel knots, source vertices and zero crossings. Integrate means and histogram
length weights analytically rather than choosing an arbitrary sample count.
This smooths turning density, not curve coordinates. The signed integral equals
total signed turning; the unsigned integral can decrease when opposite turns
cancel. Exact absolute vertex turning remains available separately.

Analytical gates include an unequal-arm 90-degree elbow, a subdivided copy,
reversal, rotation, scale transformations, an asymmetric closed triangle with
every cyclic starting vertex, a regular circle polygon, and invalid retracing.
An elbow tent with h=2 has height 45 degrees/unit, integral 90 degrees, and
squared integral 2700 degrees squared/unit. Its dual cell spans 15 units at
density 6. A simple closed contour has signed total plus or minus 360 degrees.
Signed curvature reverses sign; its absolute distribution does not.

Floating-point invariance tolerances are 1e-10 on these bounded open fixtures,
1e-9 for cyclic polygons, and 1e-8 for a 512-segment circle. They allow accumulated
double-precision coordinate arithmetic, not changed estimator definitions. At
radius 10 and h=2, mean absolute circle curvature error must be below 4% relative
to 180/(pi*10). Scaling lengths and h by c scales curvature by 1/c.

Tracing and simplification alter geometry and have sensitivity gates rather
than exact invariance. Freeze a radius-40 circle on a 101-square grid, centers
at pixel offsets 0 and 0.25, resolutions 1x and 2x, and simplification tolerances
0, 0.25, 0.5, 1 pixel. Use a physical h of 8 original pixels. Report contour
topology, relative length error, signed turn, mean absolute curvature and its
distribution for every variant. For preserved single-loop topology, require
absolute signed turn 360 within 1e-7 and mean absolute curvature within 25% of
180/(pi*40). Preserve failures and changed topology explicitly; do not tune the
grid or h after looking at results. Also retain a small-amplitude zigzag to
show why simplification tolerance alone cannot bound total absolute turning.
Smoothing does not correct staircase arc-length inflation or ambiguous traced
endpoint tangents. Convergence to the underlying Euclidean smooth curve is not
promised without a separate extraction/dequantization contract.

Branch gates use exact piecewise-linear segment integrals over fixed physical
windows: a straight axis, constant radius, affine radius, an absolute-value
radius profile, a right-angle centerline, nonuniform collinear subdivision with
interpolated radii, and zero-length/duplicate-point cases. Require pointwise
agreement at common physical positions, not equality of differently sampled
unweighted averages. Constant radius has parallelism and taper consistency 1;
an affine radius slope v has parallelism (1+v^2)^(-power/2) and taper consistency
1. Straightness is a geometric proxy, not evidence of reflection symmetry.

All gates apply to explicit new estimators. Raster MAT topology and perceptual
utility are separate claims; these tests do not establish them.
