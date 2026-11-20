% Export stored MATLAB MLV outputs from reference-data/mlv to flat CSV files.
% Usage: octave-cli --quiet --no-history tools/mlvref/export.m
% Writes gzip CSVs .agent-work/mlvref/<dataset>/{drawings,segments,junctions,mat_samples}.csv.gz;
% datasets with a COMPLETE marker are skipped. MAT scores are per pixel in MLV
% (MATpropertiesToContours rasterizes each contour), so they are exported only
% as the *_allX/_allY/_allScores samples in mat_samples.csv.gz. Natural data never leave the local disk.
gzw = @(path) popen(['gzip -1 -c > ''' path ''''], 'w');
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
src = fullfile(root, 'reference-data', 'mlv', 'dataSets');
out = fullfile(root, '.agent-work', 'mlvref');
files = dir(fullfile(src, '*', '*.mat'));
for i = 1:numel(files)
  p = fullfile(files(i).folder, files(i).name);
  [~, stem] = fileparts(p);
  d = fullfile(out, stem); if ~exist(d, 'dir'), mkdir(d); end
  if exist(fullfile(d, 'COMPLETE'), 'file'), printf('%s: done\n', stem); continue; end
  x = load(p); n = fieldnames(x); V = x.(n{1});
  hasmat = isfield(V, 'parallelism');
  fs = gzw(fullfile(d, 'segments.csv.gz'));
  fprintf(fs, 'drawing,contour,segment,x1,y1,x2,y2,orientation,length,curvature');
  props = {'parallelism', 'separation', 'mirror'};
  fprintf(fs, '\n');
  fj = gzw(fullfile(d, 'junctions.csv.gz'));
  fprintf(fj, 'drawing,junction,x,y,type,contours,segments,angles,min_angle,max_angle\n');
  fd = gzw(fullfile(d, 'drawings.csv.gz'));
  fprintf(fd, 'drawing,width,height,num_contours,line_method,original_image\n');
  if hasmat
    fa = gzw(fullfile(d, 'mat_samples.csv.gz'));
    fprintf(fa, 'drawing,property,x,y,score\n');
  end
  for k = 1:numel(V)
    v = V(k);
    fprintf(fd, '%d,%d,%d,%d,%s,%s\n', k, v.imsize(1), v.imsize(2), v.numContours, ...
            strrep(char(v.lineMethod), ',', ';'), strrep(char(v.originalImage), ',', ';'));
    for c = 1:v.numContours
      z = v.contours{c}; ns = size(z, 1);
      o = v.orientations{c}(:); l = v.lengths{c}(:); cu = v.curvatures{c}(:);
      M = [repmat([k c], ns, 1), (1:ns)', z(:, 1:4), o, l, cu];
      fmt = [strjoin(repmat({'%.17g'}, 1, size(M, 2)), ','), '\n'];
      fprintf(fs, fmt, M.');
    end
    J = v.junctions;
    for j = 1:numel(J)
      fprintf(fj, '%d,%d,%.17g,%.17g,%s,%s,%s,%s,%.17g,%.17g\n', k, j, J(j).position(1), J(j).position(2), ...
        char(J(j).type), strjoin(arrayfun(@num2str, J(j).contourIDs, 'UniformOutput', false), ';'), ...
        strjoin(arrayfun(@num2str, J(j).segmentIDs, 'UniformOutput', false), ';'), ...
        strjoin(arrayfun(@(a) sprintf('%.17g', a), J(j).angles, 'UniformOutput', false), ';'), ...
        J(j).minAngle, J(j).maxAngle);
    end
    if hasmat
      for q = 1:3
        X = v.([props{q} '_allX'])(:); Y = v.([props{q} '_allY'])(:); S = v.([props{q} '_allScores'])(:);
        fprintf(fa, ['%d,' props{q} ',%.17g,%.17g,%.17g\n'], [repmat(k, numel(X), 1), X, Y, S].');
      end
    end
  end
  pclose(fs); pclose(fj); if hasmat, pclose(fa); end; pclose(fd);
  fc = fopen(fullfile(d, 'COMPLETE'), 'w'); fprintf(fc, '%s\n', p); fclose(fc);
  printf('%s: %d drawings\n', stem, numel(V)); fflush(stdout);
end
