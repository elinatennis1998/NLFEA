function tau = RFBEfuncE(ECoordinates,Emu,Elam,minc,submeshtype)
%
% Tim Truster
% 10/13/2013
%
% Function to compute residual-free bubble for elasticity on linear
% triangular element with coordinates xl=ECoordinates and material
% coefficients Emu,Elam, with a subgrid size of minc.
% Adapted from IVMS folder.

if nargin == 4
    submeshtype = 1;
end

bmat = zeros(4,2);
tauRFB = zeros(4,4);

NCR = 1; % Input file
batchinter = 'batch';
if submeshtype == 1
batchname = 'RFBM_Elasticity'; % uniform submesh
else
batchname = 'RFBM_Elasticity2'; % graded submesh
end
NL_FEA_Program

ModelDx1 = ModelDx;
len = ECoordinates(:,1)-ECoordinates(:,2);
blen = sqrt(len'*len);

bmat(1,1) = sum(Node_U_V(4:2+minc,1))/minc*blen;
bmat(2,2) = sum(Node_U_V(4:2+minc,2))/minc*blen;
ModelDxx1 = zeros(neq,1);
ModelDxx1(1:2:neq-1) = ModelDx1(1:2:neq-1);
ModelDxy1 = zeros(neq,1);
ModelDxy1(2:2:neq) = ModelDx1(2:2:neq);
% ModelDxy1 = diag([0; xind(1:neq-1)])*ModelDx1;
tauRFB(1,1) = ModelDxx1'*Kdd11*ModelDxx1;
tauRFB(1,2) = ModelDxx1'*Kdd11*ModelDxy1;
tauRFB(2,1) = ModelDxy1'*Kdd11*ModelDxx1;
tauRFB(2,2) = ModelDxy1'*Kdd11*ModelDxy1;

% SurfacesL = [4 1 1 1 0 1 0
%              5 4 3 1 0 1 0
%              6 5 5 1 0 1 0
%              7 6 7 1 0 1 0
%              8 7 9 1 0 1 0
%              9 8 11 1 0 1 0
%              10 9 13 1 0 1 0
%              2 10 15 1 0 1 0];


NCR = 2; % Continuation file
batchinter = 'batch';
if submeshtype == 1
batchname = 'RFBC_Elasticity'; % uniform submesh
else
batchname = 'RFBC_Elasticity2'; % graded submesh
end
NL_FEA_Program

ModelDx2 = ModelDx;
len = ECoordinates(:,1)-ECoordinates(:,2);
blen = sqrt(len'*len);
% xind = zeros(neq,1);
% for i = 1:2:neq-1
% xind(i) = 1;
% end
bmat(3,1) = sum(Node_U_V(4:2+minc,1))/minc*blen;
bmat(4,2) = sum(Node_U_V(4:2+minc,2))/minc*blen;
% ModelDxx2 = diag(xind)*ModelDx2;
% ModelDxy2 = diag([0; xind(1:neq-1)])*ModelDx2;
ModelDxx2 = zeros(neq,1);
ModelDxx2(1:2:neq-1) = ModelDx2(1:2:neq-1);
ModelDxy2 = zeros(neq,1);
ModelDxy2(2:2:neq) = ModelDx2(2:2:neq);
tauRFB(3,3) = ModelDxx2'*Kdd11*ModelDxx2;
tauRFB(3,4) = ModelDxx2'*Kdd11*ModelDxy2;
tauRFB(4,3) = ModelDxy2'*Kdd11*ModelDxx2;
tauRFB(4,4) = ModelDxy2'*Kdd11*ModelDxy2;
tauRFB(1,3) = ModelDxx1'*Kdd11*ModelDxx2;
tauRFB(1,4) = ModelDxx1'*Kdd11*ModelDxy2;
tauRFB(2,3) = ModelDxy1'*Kdd11*ModelDxx2;
tauRFB(2,4) = ModelDxy1'*Kdd11*ModelDxy2;
tauRFB(3,1) = ModelDxx2'*Kdd11*ModelDxx1;
tauRFB(3,2) = ModelDxx2'*Kdd11*ModelDxy1;
tauRFB(4,1) = ModelDxy2'*Kdd11*ModelDxx1;
tauRFB(4,2) = ModelDxy2'*Kdd11*ModelDxy1;

tau = bmat'*tauRFB^-1*bmat;

% Verification:
% xlintL =
%                    0                   0  -4.000000000000000
%    1.875000000000000   2.000000000000000   1.875000000000000
% tauL = RFBEfuncE(xlintL,40,40,4,1)
% tauL = %% PLANE STRESS
%   1.0e-005 *
%    0.945582094585593   0.009222147354223
%    0.009222147354223   0.357271222860414
% tauL = %% PLANE STRAIN
%   1.0e-005 *
%    0.944612785597053   0.009828225567741
%    0.009828225567741   0.317634104118477
   