% Track B: pinned computeGradientVectorField + computeAOF under Octave on a recovered thin-ink raster.
% Usage: octave-cli mat_octave_aof.m <indir> <adapterdir>; indir has hw.csv, D.bin (stored distance), AOF.bin.
args = argv(); d = args{1}; here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
pkg load image; addpath(fullfile(root,'vendor/MLV_toolbox/MLVcode')); addpath(args{2}, '-begin');
hw = dlmread(fullfile(d,'hw.csv')); h = hw(1); w = hw(2);
rd = @(f) reshape(fread(fopen(fullfile(d,f)),Inf,'double'), h, w);
Ds = rd('D.bin'); As = rd('AOF.bin');
thin_boundary = double(~(Ds <= 0));
[D, IDX] = computeGradientVectorField(thin_boundary);
alpha = (2*pi)/60; s = (0:59)'; sp = [cos(alpha*s), sin(alpha*s)];
A = computeAOF(D, IDX, sp, 1);
printf('D class %s; max|D - stored| %g; AOF max|diff| %g; pixels |diff|>1e-9: %d\n', class(D), max(abs(double(D(:)) - Ds(:))), max(abs(A(:) - As(:))), nnz(abs(A(:) - As(:)) > 1e-9));
fid = fopen(fullfile(d,'octave_idx.bin'),'w'); fwrite(fid, IDX, 'double'); fclose(fid);
