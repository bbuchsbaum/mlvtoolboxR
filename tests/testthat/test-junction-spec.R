test_that("feature specifications carry the junction rule and gap metric", {
  s <- feature_spec("junction_type")
  expect_equal(s$settings$junction$rule, "exit")
  expect_equal(s$settings$junction$gap_metric, "euclidean")
  expect_equal(unname(s$definitions["junction"]), "exit-neighbourhood-1-euclidean")
  m <- feature_spec("junction_type", junction = list(rule = "mlv2023"))
  expect_equal(m$settings$junction$gap_metric, "chebyshev")
  expect_equal(unname(m$definitions["junction"]), "mlv2023-port-f12418a")
  expect_error(feature_spec(junction = list(rule = "mlv2023", gap_metric = "euclidean")), "chebyshev")
  expect_error(feature_spec(junction = list(rule = "other")))
  expect_silent(validate_feature_spec(m))
  old <- s; old$definitions <- old$definitions[names(old$definitions) != "junction"]
  expect_error(validate_feature_spec(old), "noncanonical")
  # The specification's rule reaches detect_junctions.
  d <- mlv_drawing(list(matrix(c(0, 0, 10, 10), 1), matrix(c(5, 5, 10, 0), 1)), 20, 20)
  a <- contour_statistics(d, spec = feature_spec("junction_angle", bins = 4, normalize = FALSE))
  b <- contour_statistics(d, spec = feature_spec("junction_angle", bins = 4, normalize = FALSE, junction = list(rule = "mlv2023")))
  # Diagonal T: minimum arm gap 90 (exit) vs 0 under MLV-M007; the
  # junction_angle histogram bins the minimum angle over [0, 120].
  expect_equal(a$value, c(0, 0, 0, 1))
  expect_equal(b$value, c(1, 0, 0, 0))
})
