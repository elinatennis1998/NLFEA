GrainSig(1:length(ElemS),grain) = GrainSig(1:length(ElemS),grain) + ElemS;
GrainEps(1:length(ElemE),grain) = GrainEps(1:length(ElemE),grain) + ElemE;
GrainVol(:,grain) = GrainVol(:,grain) + ElemV;

%Modified 6/27/2019
GrainDisp(1:length(ElemDisp),grain) = GrainDisp(1:length(ElemDisp),grain) + ElemDisp;