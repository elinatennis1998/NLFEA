%

NodesOnElement = NodesOnElementg;
ix = [NodesOnElement RegionOnElement];
MatTypeTable=MatTypeTableg;
nen = neng;
iedof2 = zeros(2,nen,30);
iedof2(:,1:nen/2,:) = iedof(:,1:nen/2,:);
iedof2(:,nen/2+1:nen,:) = iedof(:,1:nen/2,:);
iedof = iedof2;

% Assemble CG K
FormCGK
GrainInteg

InterInteg
% ModelDx'*Kdd11*ModelDx/GrainVol(1,5)-GrainSig(:,5)'*GrainEps(:,5)/GrainVol(1,5)^2
% interior=GrainVol(2,5)/GrainVol(1,5)-GrainSig(1:3,5)'*GrainEps(:,5)/GrainVol(1,5)^2
inters = [8 9 12 16];
signs = [1 1 2 2];[2 2 1 1];
boundterm = 0;
for l = 1:4;
    i = inters(l);
    j = signs(l);
    boundterm = boundterm + BoundEps(1+j,i);
end
boundary=boundterm/GrainVol(5)
boundary2=BoundEps(2,10)/GrainVol(5)
boundary3=BoundEps(3,10)/GrainVol(5)