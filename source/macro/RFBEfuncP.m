function tau = RFBEfuncP(ECoordinates,EA,minc,submeshtype)
%
% Tim Truster
% 10/13/2013
%
% Function to compute residual-free bubble for Poisson equation on linear
% triangular element with coordinates xl=ECoordinates and material
% coefficients EA, with a subgrid size of minc.
% Adapted from IVMS folder.

if nargin == 3
    submeshtype = 1;
end

bmat = zeros(4,2);
tauRFB = zeros(4,4);

NCR = 1; % Input file
batchinter = 'batch';
if submeshtype == 1
batchname = 'RFBM_Poisson'; % uniform submesh
else
batchname = 'RFBM_Poisson2'; % graded submesh
end
NL_FEA_Program

% plotModelContP(Coordinates, Node_U_V, ix, numel, ndm, 1, 1, 1, '',0,[3 4 6 9 0])
ModelDx1 = ModelDx;
len = ECoordinates(:,1)-ECoordinates(:,2);
blen = sqrt(len'*len);
% xind = zeros(neq,1);
% for i = 1:2:neq-1
% xind(i) = 1;
% end
if submeshtype == 1
bmat = sum(Node_U_V(4:2+minc,1))/minc*blen;
else
bmat = sum(Node_U_V(2:minc,1))/minc*blen;
end
tauRFB = ModelDx'*Kdd11*ModelDx;

tau = bmat'*tauRFB^-1*bmat;

% Verification:
% AvL =
%      1     1     0
% xlintL =
%    0.100000000000000   0.125000000000000   0.200000000000000
%    0.250000000000000   0.250000000000000   0.312500000000000
% tauL = RFBEfuncP(xlintL,AvL,8,1)
% tauL =
%     7.231408809117955e-005
% tauL = RFBEfuncP(xlintL,AvL,8,2)
% tauL =
%     9.413775286993840e-005
