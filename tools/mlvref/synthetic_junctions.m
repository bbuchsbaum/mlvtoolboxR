% Run the pinned MLV junction pipeline on synthetic cases.
% Usage: octave-cli --quiet --no-history tools/mlvref/synthetic_junctions.m <cases.csv> <out.csv>
% cases.csv columns: case,contour,x1,y1,x2,y2 (numeric, header row).
args = argv();
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'vendor', 'MLV_toolbox', 'MLVcode'));
M = dlmread(args{1}, ',', 1, 0);
fo = fopen(args{2}, 'w');
fprintf(fo, 'case,junction,x,y,type,contours,segments,angles,min_angle,max_angle,error\n');
for cs = unique(M(:, 1))'
  R = M(M(:, 1) == cs, :);
  v = struct();
  ids = unique(R(:, 2), 'stable');
  v.numContours = numel(ids);
  v.contours = cell(1, numel(ids));
  for c = 1:numel(ids)
    v.contours{c} = R(R(:, 2) == ids(c), 3:6);
  end
  v.imsize = [100000 100000];
  try
    v = computeContourProperties(v, {'orientation', 'length', 'curvature', 'junctions'});
    J = v.junctions;
    if isempty(J)
      fprintf(fo, '%d,0,NaN,NaN,none,,,,NaN,NaN,\n', cs);
    end
    for j = 1:numel(J)
      fprintf(fo, '%d,%d,%.17g,%.17g,%s,%s,%s,%s,%.17g,%.17g,\n', cs, j, J(j).position(1), J(j).position(2), ...
        char(J(j).type), strjoin(arrayfun(@num2str, J(j).contourIDs, 'UniformOutput', false), ';'), ...
        strjoin(arrayfun(@num2str, J(j).segmentIDs, 'UniformOutput', false), ';'), ...
        strjoin(arrayfun(@(a) sprintf('%.17g', a), J(j).angles, 'UniformOutput', false), ';'), ...
        J(j).minAngle, J(j).maxAngle);
    end
  catch err
    fprintf(fo, '%d,0,NaN,NaN,error,,,,NaN,NaN,%s\n', cs, strrep(err.message, ',', ';'));
  end
end
fclose(fo);
