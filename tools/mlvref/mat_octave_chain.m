% Track B: run pinned MLV tracing/scoring/mapping on supplied MAT stage fields.
% Usage: octave-cli mat_octave_chain.m <indir> ; indir holds h,w,D.bin,AOF.bin,S.bin,pix.bin (int32 x,y pairs)
% Writes branches.csv (branch,row,col,radius) and mapped_<prop>.bin (double at pix order).
args = argv(); d = args{1};
here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
pkg load image
addpath(fullfile(root,'vendor/MLV_toolbox/MLVcode'));
addpath(genpath(fullfile(root,'vendor/MLV_toolbox/3rdParty/contour_tracing')));
addpath(args{2}, '-begin');   % adapter directory (knnsearch, smoothdata, gradient field)
hw = dlmread(fullfile(d,'hw.csv')); h = hw(1); w = hw(2);
rd = @(f,t) reshape(fread(fopen(fullfile(d,f)),Inf,t), h, w);
MAT.distance_map = single(rd('D.bin','double'));; MAT.AOF = rd('AOF.bin','double'); MAT.skeleton = logical(rd('S.bin','uint8'));
fid = fopen(fullfile(d,'pix.bin')); P = reshape(fread(fid,Inf,'int32'),2,[])'; fclose(fid);
img = ones(h,w); img(sub2ind([h w], P(:,2), P(:,1))) = 0;
B = traceSkeleton(MAT);
fb = fopen(fullfile(d,'branches.csv'),'w');
for i=1:numel(B), fprintf(fb,'%d,%d,%d,%.17g\n',[repmat(i,numel(B(i).X),1), B(i).X, B(i).Y, double(B(i).Radius)]'); end; fclose(fb);
props = {'parallelism','separation','mirror'};
for q = 1:3
  [sk, B] = computeMATproperty(MAT, props{q}, B);
  C = mapMATtoContour(B, img, sk);
  v = C(sub2ind([h w], P(:,2), P(:,1)));
  fo = fopen(fullfile(d,['mapped_' props{q} '.bin']),'w'); fwrite(fo, v, 'double'); fclose(fo);
  fs = fopen(fullfile(d,['skel_' props{q} '.bin']),'w'); fwrite(fs, sk, 'double'); fclose(fs);
end
printf('branches %d\n', numel(B));
