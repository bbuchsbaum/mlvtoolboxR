% Track B: pinned MLV GetConSeg (Kovesi edgelink + cleanedgelist, minlength 10) on
% synthetic skeleton masks. Usage: octave-cli mat_edgelink_fixture.m <adapterdir>
% The adapter directory replaces bwmorph(...,'skel',Inf) inside edgelink/findendsjunctions
% by identity: the masks are fixed points of MATLAB-equivalent thinning (checked in R).
args = argv(); here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
pkg load image
addpath(genpath(fullfile(root,'vendor/MLV_toolbox/3rdParty/contour_tracing')));
addpath(args{1}, '-begin');
d = fullfile(root,'inst','fixtures','edgelink');
for nm = {'grid','fan','kite'}
  m = logical(dlmread(fullfile(d,[nm{1} '.mask.csv']),','));
  L = GetConSeg(m);
  f = fopen(fullfile(d,[nm{1} '.octave.branches.csv']),'w');
  for i = 1:numel(L), fprintf(f,'%d,%d,%d\n',[repmat(i,size(L{i},1),1), L{i}]'); end
  fclose(f); printf('%s: %d branches\n', nm{1}, numel(L));
end
