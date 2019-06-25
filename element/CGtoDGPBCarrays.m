% 04/19/2013
% Extract Left/Right arrays from coupled array

maL = mateprop(1);
maR = mateprop(2);
nelL = mateprop(3);
nelR = mateprop(4);
Coeffs = sign(mateprop(5:4+ndm));
nstL = nelL*ndf;
nstR = nelR*ndf;
nstM = ndm*ndf;
ElemFlagL = ElemFlag(1:nelL);
ElemFlagR = ElemFlag(nelL+1:nelL+nelR);
xlL = xl(1:ndm,1:nelL);
xlR = xl(1:ndm,nelL+1:nelL+nelR);
xlM = xl(1:ndm,nelL+nelR+1:nelL+nelR+ndm); 
%  xlM = xl(1:ndm,15:16); 
if isw > 1 && isw ~= 53 %Modificaion of 4/17/2019
ulL = ul(1:ndf,1:nelL);
ulR = ul(1:ndf,nelL+1:nelL+nelR);
ulM = ul(1:ndf,nelL+nelR+1:nelL+nelR+ndm);
end
matepropL = MateT(maL,:);
matepropR = MateT(maR,:);