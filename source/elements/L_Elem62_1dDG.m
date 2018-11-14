%9/14/2016
%

%MateT=[1 1 1
%     1 1 1
%     1 1 1];
% xl=[1 2]
%
switch isw %Task Switch
%%

    case {3,6,21}
%Modified Subroutine to fit the input file
nelL=2;
nelR=2;
ndf=1;
p=1;
nelL; %number of elemental nodes in the left element
nelR; %number of elemental nodes in the right element
% KLL=zeros(nelL);
% KLR=zeros(nelL,nelR);
% KRL=zeros(nelR,nelL);
% KRR=zeros(nelR,nelR);

PatchEL = mateprop(1);
PatchER = mateprop(2);
PatchAL = mateprop(3);
PatchAR = mateprop(4);
  
DmatL = PatchEL * PatchAL;
DmatR = PatchER * PatchAR;
NmatL = zeros(1,nelL);
BmatL = zeros(1,nelL);
NmatR = zeros(1,nelR);
BmatR = zeros(1,nelR);
%Gauss guadrature for linear element

basisL = zeros(1,p+1);
basisR = zeros(1,p+1);
derbasis = zeros(1,p+1);
nP=1; %number of points
xL = 1; %SamplePoint
xR = -1; %SamplePoint
nL = 1; % outward normal on left element
nR = -1; % outward normal on right element
w = 1; %weight
basisL(1,1) = (1-xL)/2;
basisL(1,2) = (1+xL)/2;
basisR(1,1) = (1-xR)/2;
basisR(1,2) = (1+xR)/2;
derbasis(1,1) = -1/2;
derbasis(1,2) = 1/2;
derbasisL=derbasis;
derbasisR=derbasis;
JacbnL = derbasisL * xl(1,1:nelL)';
JacbnR = derbasisR * xl(1,nelL+1:nelL+nelR)';
NmatL = basisL;
BmatL =1/JacbnL*(derbasisL);
NmatR = basisR;
BmatR = derbasisR/JacbnR;
% ElemK = zeros(nel)
% ElemF = zeros(nel+1,1)
% KvLL = w*BmatL'*DmatL*BmatL*Jacbn;
% KvLR = zeros(nelL);
% KvRL = zeros(nelR);
% KvRR = KvLL;
% ElemK=[KvLL KvLR
%        KvRL KvRR];
%% for the interface stiffness
%% for the first integral term
KLL1 = -1/2*w*NmatL'*DmatL*BmatL;
KLR1 = -1/2*w*NmatL'*DmatR*BmatR;
KRL1 = 1/2*w*NmatR'*DmatL*BmatL;
KRR1 = 1/2*w*NmatR'*DmatR*BmatR;
Ks1 = [KLL1 KLR1
       KRL1 KRR1];
%% for the second integral
theta= -1; % symmetric
% theta= 1; % nonsymmetric 
Ks2=-theta*(Ks1');
%% for the third integral
delGoS=10000;1/2;0; % YOU NEED A BETTER FORMULA FOR DELGOS
KLL3 = 1/2*delGoS*w*NmatL'*NmatL;
KLR3 = -1/2*delGoS*w*NmatL'*NmatR;
KRL3 = -1/2*delGoS*w*NmatR'*NmatL;
KRR3 = 1/2*delGoS*w*NmatR'*NmatR;
Ks3 = [KLL3 KLR3
       KRL3 KRR3];
ElemK = Ks1+Ks2+Ks3;


end










