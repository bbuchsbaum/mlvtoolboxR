# MATLAB capability map

All 62 top-level MATLAB functions are accounted for below. Functions are
consolidated into R workflows rather than exported as one-for-one wrappers.
Mapping a capability does not assert identical morphology, rasterization,
tracing or model behavior; see DIFFERENCES.md and IMPLEMENTATION_STATUS.md.

| MATLAB function | R implementation / scope |
|---|---|
| `InitializeNeighborhoods` | Native Pratt LUT and pixel-graph kernels |
| `MATpropertiesToContours` | map_mat_properties |
| `allLDHistogramsToTable` | batch_statistics / feature_table |
| `applyCircularAperture` | clip_circle |
| `averageProperty` | average_property |
| `cleanupJunctions` | detect_junctions (transitive union-find) |
| `computeAOF` | average_outward_flux |
| `computeAllMATfromVecLD` | analyse_mat |
| `computeAllMATproperties` | compute_mat_properties |
| `computeColorIndex` | plot / render_drawing palettes |
| `computeContourProperties` | compute_contour_properties |
| `computeCurvature` | as.data.frame / compute_contour_properties |
| `computeGradientVectorField` | compute_mat / distance_transform |
| `computeJunctionAnglesTypes` | detect_junctions |
| `computeJunctions` | detect_junctions |
| `computeLength` | as.data.frame / compute_contour_properties |
| `computeMAT` | compute_mat |
| `computeMATproperty` | compute_mat_properties |
| `computeMATpropertyPerBranch` | score_branch |
| `computeOrientation` | as.data.frame / compute_contour_properties |
| `detectJunctions` | detect_junctions (sweep or exhaustive) |
| `drawAllProperties` | plot_properties |
| `drawContoursByStatsModel` | plot selected contours / render_drawing(values=predict_contours(...)) |
| `drawJunctions` | plot(property='junctions') |
| `drawLinedrawing` | plot |
| `drawLinedrawingProperty` | plot(property=...) |
| `drawMATproperty` | plot(property=...) with local pixel scores |
| `extract2DSkeletonFromBinaryImage` | compute_mat / skeletonize |
| `fitLineSegments` | fit_line_segments |
| `generateFeatureDensityMap` | feature_density |
| `getContourPropertiesStats` | contour_statistics |
| `getCurvatureStats` | contour_statistics(properties='curvature') |
| `getDistanceFromLineSegment` | fit_line_segments internal distance predicate |
| `getHorizontalVerticalStats` | contour_statistics(properties='horizontal_vertical') |
| `getJunctionStats` | contour_statistics(properties=c('junction_type','junction_angle')) |
| `getLengthStats` | contour_statistics(properties='length') |
| `getMATpropertyStats` | contour_statistics(properties=...) |
| `getOrientationStats` | contour_statistics(properties='orientation') |
| `histogramToTable` | contour_statistics / feature_table |
| `importSVG` | read_mlv_svg (explicit support/sampling contract) |
| `labelMatrix2edges` | labels_to_edges (explicit two-sided boundary convention) |
| `lineIntersection` | Native exact intersection predicate in detect_junctions |
| `mapMATtoContour` | map_mat_properties |
| `mergeLineSegments` | simplify_contours (no silent removal of merely nearby contours) |
| `predictContoursByStatsModel` | predict_contours |
| `randomlyShiftContours` | shift_contours |
| `removeZeroLengthContours` | mlv_drawing(zero_length='drop') |
| `renderContoursByStatsModel` | render_drawing(values=predict_contours(...)) |
| `renderJunctions` | render_drawing(property='junctions') |
| `renderLinedrawing` | rasterize / render_drawing |
| `renderLinedrawingProperty` | render_drawing(property=...) |
| `renderMATproperty` | render_drawing(property=...) with local score image |
| `rotateLinedrawing` | rotate_drawing |
| `saveSceneLDsToTable` | batch_statistics / write_feature_csv |
| `segmentContoursAtJunctions` | segment_at_junctions |
| `splitLDbyHistogramWeights` | split_by_histogram |
| `splitLDbyProperties` | split_contours |
| `splitLDbyStatsModel` | split_contours(scores=predict_contours(...)) |
| `splitLDmiddleSegmentsVsJunctions` | split_junction_regions |
| `traceLineDrawingFromRGB` | trace_rgb (explicit provider; bundled SAM/legacy forest not included) |
| `traceLinedrawingFromEdgeMap` | trace_contours |
| `traceSkeleton` | trace_skeleton |
