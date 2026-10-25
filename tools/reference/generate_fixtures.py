#!/usr/bin/env python3
"""Regenerate small synthetic fixtures with the pinned Octave reference."""
from pathlib import Path
import hashlib,json,subprocess,tempfile,datetime
root=Path(__file__).resolve().parents[2]
source=root/'vendor/MLV_toolbox/MLVcode'
out=root/'inst/fixtures/octave';out.mkdir(parents=True,exist_ok=True)
text=(source/'computeGradientVectorField.m').read_text()
assert text.count('jj++m_Neighbors8')==2
with tempfile.TemporaryDirectory(prefix='mlv-reference-') as td:
    (Path(td)/'computeGradientVectorField.m').write_text(text.replace('jj++m_Neighbors8','jj + m_Neighbors8'))
    expr=f"""
pkg load image;
addpath('{source}');addpath('{td}','-begin');
out='{out}';
m=false(13,17);m(3:11,4:14)=true;
dlmwrite(fullfile(out,'pratt_input_box.csv'),m,',','precision',17);
dlmwrite(fullfile(out,'pratt_output_box.csv'),bwmorph(m,'skel',Inf),',','precision',17);
m=false(13,17);m(2:12,2:16)=true;m(5:9,5:13)=false;
dlmwrite(fullfile(out,'pratt_input_ring.csv'),m,',','precision',17);
dlmwrite(fullfile(out,'pratt_output_ring.csv'),bwmorph(m,'skel',Inf),',','precision',17);
rand('seed',41);m=rand(15,19)>.55;
dlmwrite(fullfile(out,'pratt_input_random.csv'),m,',','precision',17);
dlmwrite(fullfile(out,'pratt_output_random.csv'),bwmorph(m,'skel',Inf),',','precision',17);
[rr,cc]=ndgrid(1:15,1:19);D=sqrt((rr-8).^2+(cc-10).^2);
IDX=repmat(sub2ind(size(D),8,10),size(D));
circle=[cos(2*pi*(0:59)'/60),sin(2*pi*(0:59)'/60)];
F=computeAOF(D,IDX,circle,1);
dlmwrite(fullfile(out,'aof_distance.csv'),D,',','precision',17);
dlmwrite(fullfile(out,'aof_index.csv'),IDX,',','precision',17);
dlmwrite(fullfile(out,'aof_output.csv'),F,',','precision',17);
img=ones(25,35);img(4:22,8)=0;img(4:22,27)=0;
MAT=computeMAT(img,28);
[DD,II]=computeGradientVectorField(double(~bwmorph(~img,'skel',Inf)));
dlmwrite(fullfile(out,'mat_index.csv'),II,',','precision',17);
dlmwrite(fullfile(out,'mat_input.csv'),img,',','precision',17);
dlmwrite(fullfile(out,'mat_distance.csv'),MAT.distance_map,',','precision',17);
dlmwrite(fullfile(out,'mat_aof.csv'),MAT.AOF,',','precision',17);
dlmwrite(fullfile(out,'mat_skeleton.csv'),MAT.skeleton,',','precision',17);
vecLD(1).imsize=[20,15];vecLD(1).numContours=2;
vecLD(1).contours={{[1,1,11,1],[2,2,2,7;2,7,8,7]}};vecLD(1).originalImage='synthetic';
vecLD(2).imsize=[10,12];vecLD(2).numContours=1;vecLD(2).contours={{[1,1,5,5]}};vecLD(2).originalImage='second';
save('-mat-binary',fullfile(out,'synthetic_drawings.mat'),'vecLD');
fprintf('Octave %s; synthetic fixtures written.\\n',version);
"""
    subprocess.run(['octave-cli','--quiet','--no-history','--eval',expr],check=True)
files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!='MANIFEST.json'}
(root/'inst/fixtures/octave/MANIFEST.json').write_text(json.dumps({
 'generated':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'source':'f12418a5e03c2f38e222fed1dffccf42b227cd29','runtime':'Octave 11.1.0; image 2.20.0',
 'adaptation':'Two jj++m_Neighbors8 expressions changed to jj + m_Neighbors8 in a temporary copy.',
 'scope':'Synthetic Pratt, scalar AOF, signed-distance/MAT and MAT import fixtures; no natural-image data.',
 'sha256':files},indent=2)+'\n')
