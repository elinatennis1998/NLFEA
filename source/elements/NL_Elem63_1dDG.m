%9/14/2016

switch isw %Task Switch
%%
case {3,6,21}
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
%zeta =  mateprop(5);
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
% xl=Coordinates
JacbnL = derbasisL * xl(1,1:nelL)';
JacbnR = derbasisR * xl(1,nelL+1:nelL+nelR)';
NmatL = basisL;
BmatL =1/JacbnL*(derbasisL);
NmatR = basisR;
BmatR = derbasisR/JacbnR;
% ElemK = zeros(nel);
% ElemF = zeros(nel);
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
theta=-1; % symmetric
% theta= 1; % nonsymmetric 
Ks2=-theta*(Ks1');
%% for the third integral
sigmaC=6.25; %maximum stress at the interface
delC=0.1; % maximum separation in the TSL, the rato should be related to r
Hc=sigmaC/delC; % slope of TSL
r=1000*Hc; %penalty parameter
%r= 10000000;1/2;0; % YOU NEED A BETTER FORMULA FOR r
KLL3 = 1/2*r*w*NmatL'*NmatL;
KLR3 = -1/2*r*w*NmatL'*NmatR;
KRL3 = -1/2*r*w*NmatR'*NmatL;
KRR3 = 1/2*r*w*NmatR'*NmatR;
Ks3 = [KLL3 KLR3
       KRL3 KRR3];
ElemK = Ks1+Ks2+Ks3;
%ElemF=-ElemK*ul
%Bmat=[BmatL BmatR]
% epsil = Bmat*ul;
% sigma = Dmat*epsil
%% Displace on left and right element
ulL=ul(1,1:nelL); % displacement at the left element
ulR=ul(1,nelL+1:nelL+nelR); % displacement of the right element
%%        
%% obtain stress
sigmaL = DmatL*BmatL*ulL'; %stress at the leftelement
sigmaR = DmatR*BmatR*ulR'; %stress at the right element
avg_sigma=1/2*(sigmaL+sigmaR)% average stress at the interface
U_Jump=ul(nelL+1)-ul(nelL)%displacement jump
tau= subplus(avg_sigma*nL+r*U_Jump);% traction input
% calculate the fracture energy
sigmaC=6.25; %maximum stress at the interface
delC=0.1; % maximum separation in the TSL, the rato should be related to r
Gc=1/2*sigmaC*delC % fracture Energy
K=sigmaC*tau %max(last yield criterion)
if zeta == 0 % perfect adhesion
    zeta1=0 % zeta1 is the new zeta
    elseif K >0 && tau > 2*Gc % crumpling
    zeta1=tau/r
    elseif zeta >0 && K <= 2*Gc % damage
    A=(2*Gc)/(2*Gc*r-sigmaC^2);
    zeta1=(tau-sigmaC)*A;
end
SmallPhi=sigmaC*zeta1 % new yield criterion
if zeta1> 0 & SmallPhi<K % unloading
    zeta1=(tau/r)-(zeta1*sigmaC/(zeta*r))-(sigmaC/(2*Gc*r))
    % linearization using Newton_Raphson
    zeta1=( (tau-r*zeta1+(sigmaC-sigmaC^2*zeta1/(2*Gc)))+(r+sigmaC^2/(2*Gc))*(zeta1))/(r+sigmaC^2/(2*Gc))
end
zeta1
%%
%%%%%%%% obtain the force terms
av_tr=1/2*(DmatL*BmatL*ulL'+DmatR*BmatR*ulR');% Average traction
%% obtain the first force term
F1L=-w*NmatL'*av_tr; %force on the left element
F1R=w*NmatR'*av_tr; %force on the second element
F_1term=[F1L;F1R]; % arrange the force term
%% Obtain the second force term
U_Jump=ul(nelL+1)-ul(nelL); % obtain displacement jump [[u]]
Jump_Gap=U_Jump-zeta1;% displacement jump - gap ([[u]]-z)
F2L=-1/2*w*BmatL'*DmatL*Jump_Gap; % force on the left element
F2R=1/2*w*BmatR'*DmatR*Jump_Gap;% force on the left element
F_2term=[F2L;F2R]; % % arrange the force term
%% obtain the third force term ( the penalty term)
F3L=-1/2*r*w*NmatL'*Jump_Gap;
F3R=1/2*r*w*NmatR'*Jump_Gap;
F_3term=[F3L;F3R];% third force term
Fint=F_1term+F_2term+F_3term;
ElemF=-Fint
disp=ul
end
