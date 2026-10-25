% Minimal source-level reproducers. Run with: octave-cli --quiet --no-history tools/upstream/reproduce_matlab.m
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'vendor','MLV_toolbox','MLVcode'));
v.imsize=[30,30];v.numContours=1;v.originalImage='synthetic';v.lineMethod='synthetic';
v.contours={[1,1,11,1]};
v=computeContourProperties(v,{'orientation','length','curvature'});
[v,h,b]=getHorizontalVerticalStats(v);
fprintf('MLV-M001 horizontal histogram total: observed=%g expected=10\n',sum(h));
v.junctions=[];[v,h,b]=getJunctionStats(v);
fprintf('MLV-M002 empty junction angle bins: observed=%d expected=8\n',numel(v.junctionAngleHistogram));
w=v;w.numContours=3;w.contours={[1,1,5,1],[1,10,5,10],[1,20,3,20]};
w=computeContourProperties(w);[top,bottom]=splitLDbyProperties(w,'length',.5);
fprintf('MLV-M003 top split length: observed=%g maximum=5\n',sum(top.contourLengths));
v.imsize=[10,10];v.numContours=1;v.contours={[0,5,10,5]};
clipped=applyCircularAperture(v,2);clipped=computeLength(clipped);
fprintf('MLV-M005 outside-outside aperture length: observed=%g expected=4\n',sum(clipped.contourLengths));
for i=1:3
 j(i).position=[0,0];j(i).contourIDs=[2*i-1,2*i];j(i).segmentIDs=[1,1];
end
j(2).position=[4,0];j(3).position=[2,0];
clean=cleanupJunctions(j,2);
fprintf('MLV-M006 transitive junction groups: observed=%d expected=1\n',numel(clean));
v.imsize=[20,20];v.numContours=2;v.contours={[0,0,10,10],[5,5,10,0]};
v=computeContourProperties(v);
fprintf('MLV-M007 diagonal T minimum angle: observed=%g expected=90\n',v.junctions(1).minAngle);
position=lineIntersection([-5,-1,-1,-5],[-5,-5,-1,-1],0,0);
fprintf('MLV-M008 negative-coordinate intersection count: observed=%d expected=2 (coordinates -3,-3)\n',numel(position));
try
 [lines,dists,scores]=fitLineSegments([0,0;1,0]);
 fprintf('MLV-M009 two-point line fitting: unexpectedly succeeded\n');
catch err
 fprintf('MLV-M009 two-point line fitting error: %s\n',err.message);
end
fprintf('MLV-M004 threshold conversion: omitted=0.3 explicit28=%.17g\n',2/pi*sind(28));
fprintf('MLV-M010 parse incompatibility: computeGradientVectorField contains two jj++m_Neighbors8 expressions; see reference adapter.\n');
