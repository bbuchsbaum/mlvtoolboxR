% P1 version drift: recompute geometry and junctions with the pinned MLV
% MATLAB source (vendor/MLV_toolbox, f12418a) on stored contours.
% Usage: octave-cli --quiet --no-history tools/mlvref/drift_pinned.m [per_dataset] [max_segments]
% For every dataset .mat file, selects per_dataset drawings (default 4) at
% evenly spaced indices round(linspace(1,n,k)); drawings with more than
% max_segments segments (default 1500) are replaced by the nearest smaller
% index below the cap, so the sample is deterministic and bounded in cost.
% Writes .agent-work/track-a/pinned/<dataset>/{segments,junctions,sample}.csv
% in the export.m schema. Natural data stay local; only aggregates are committed.
args = argv();
k = 4; cap = 1500;
if numel(args) >= 1, k = str2double(args{1}); end
if numel(args) >= 2, cap = str2double(args{2}); end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'vendor', 'MLV_toolbox', 'MLVcode'));
src = fullfile(root, 'reference-data', 'mlv', 'dataSets');
out = fullfile(root, '.agent-work', 'track-a', 'pinned');
files = dir(fullfile(src, '*', '*.mat'));
strip = {'orientations', 'lengths', 'contourLengths', 'curvatures', 'junctions'};
for i = 1:numel(files)
  p = fullfile(files(i).folder, files(i).name);
  [~, stem] = fileparts(p);
  d = fullfile(out, stem);
  if exist(fullfile(d, 'COMPLETE'), 'file'), printf('%s: done\n', stem); continue; end
  if ~exist(d, 'dir'), mkdir(d); end
  x = load(p); n = fieldnames(x); V = x.(n{1});
  nseg = arrayfun(@(v) sum(cellfun(@(c) size(c, 1), v.contours)), V);
  want = unique(round(linspace(1, numel(V), k)));
  pick = [];
  for w = want
    ok = find(nseg <= cap);
    ok = setdiff(ok, pick);
    if isempty(ok), continue; end
    [~, m] = min(abs(ok - w) + 0.5 * (ok > w)); % prefer the lower index on ties
    pick(end+1) = ok(m);
  end
  fs = fopen(fullfile(d, 'segments.csv'), 'w');
  fprintf(fs, 'drawing,contour,segment,orientation,length,curvature\n');
  fj = fopen(fullfile(d, 'junctions.csv'), 'w');
  fprintf(fj, 'drawing,junction,x,y,type,contours,segments,angles,min_angle,max_angle\n');
  fm = fopen(fullfile(d, 'sample.csv'), 'w');
  fprintf(fm, 'drawing,segments,seconds\n');
  for kk = pick
    v = V(kk);
    for f = strip
      if isfield(v, f{1}), v = rmfield(v, f{1}); end
    end
    t0 = tic;
    v = computeContourProperties(v, {'orientation', 'length', 'curvature', 'junctions'});
    el = toc(t0);
    for c = 1:v.numContours
      ns = size(v.contours{c}, 1);
      M = [repmat([kk c], ns, 1), (1:ns)', v.orientations{c}(:), v.lengths{c}(:), v.curvatures{c}(:)];
      fprintf(fs, '%d,%d,%d,%.17g,%.17g,%.17g\n', M.');
    end
    J = v.junctions;
    for j = 1:numel(J)
      fprintf(fj, '%d,%d,%.17g,%.17g,%s,%s,%s,%s,%.17g,%.17g\n', kk, j, J(j).position(1), J(j).position(2), ...
        char(J(j).type), strjoin(arrayfun(@num2str, J(j).contourIDs, 'UniformOutput', false), ';'), ...
        strjoin(arrayfun(@num2str, J(j).segmentIDs, 'UniformOutput', false), ';'), ...
        strjoin(arrayfun(@(a) sprintf('%.17g', a), J(j).angles, 'UniformOutput', false), ';'), ...
        J(j).minAngle, J(j).maxAngle);
    end
    fprintf(fm, '%d,%d,%.3f\n', kk, nseg(kk), el);
    printf('%s drawing %d (%d segments): %.1fs\n', stem, kk, nseg(kk), el); fflush(stdout);
  end
  fclose(fs); fclose(fj); fclose(fm);
  fc = fopen(fullfile(d, 'COMPLETE'), 'w'); fprintf(fc, '%s\n', p); fclose(fc);
end
