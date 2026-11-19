# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

## Stack

Local single-page app shipped inside the mlvtoolbox R package (`inst/studio/`), served by `httpuv` from `mlv_studio()`. TypeScript, React, Vite, Tailwind (matches psycloud). Not Shiny. See docs/plans/mlv-studio-spec.md.

## Users

Vision scientists and their students, at a desk, building and measuring line-drawing stimuli for psychophysics and neuroimaging experiments, moving back and forth between R scripts and the UI, sometimes for hours. Equally used to explain and teach how mid-level vision measurements work: in talks, papers, methods sections, and courses.

## Purpose

Make contour geometry (orientation, length, curvature, junctions) and medial-axis symmetry (parallelism, separation, mirror, taper) of line drawings visible, inspectable, and controllable, and let researchers construct matched stimulus sets whose differences they can quantify.

## Mechanism

R is the only compute engine; the UI renders and drives it. Every number shown equals what the R function returns. Every UI action is an R call, and the session exports as a runnable R script.

## Modules

Library/import, tracing workbench, drawing canvas, medial-axis view, statistics, stimulus lab, provenance/R code, estimator lab (later).

## Constraints

- Stimuli are judged on a white ground; the drawing surface must render stimuli faithfully (no tinting of the stimulus itself unless a colour encoding is chosen).
- Measurement colour encodings need legends with units; orientation is cyclic (mod 180°).
- Built bundle ≤ 1.5 MB; no external network at runtime.
- Implements methods of Walther et al. (2023); published estimators must remain visible as such.

## Voice

Precise, technical, peer-to-peer. Units always shown. No hype.
