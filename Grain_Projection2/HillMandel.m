%
batchinter = 'batch';
batchname = 'SquareRVE.m';
fprintf('starting file %s\n',batchname)
FEA_Program

NodesOnElement = NodesOnElementg;
ix = [NodesOnElement RegionOnElement];
MatTypeTable=MatTypeTableg;
nen = neng;
iedof2 = zeros(2,8,21);
iedof2(:,1:4,:) = iedof(:,1:4,:);
iedof2(:,5:8,:) = iedof(:,1:4,:);
iedof = iedof2;

% Assemble CG K
RegList = 5;
FormCGK
GrainInteg

InterInteg
% ModelDx'*Kdd11*ModelDx/GrainVol(1,5)-GrainSig(:,5)'*GrainEps(:,5)/GrainVol(1,5)^2
GrainVol(2,5)/GrainVol(1,5)-GrainSig(:,5)'*GrainEps(:,5)/GrainVol(1,5)^2
inters = [4 5 7 9];
signs = [2 2 1 1];
boundterm = 0;
for l = 1:4;
    i = inters(l);
    j = signs(l);
    boundterm = boundterm + BoundEps(1+j,i);
end
boundterm/GrainVol(5)