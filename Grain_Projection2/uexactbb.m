function [ue,duex,duey] = uexactbb(XX,YY,EX,VV,D)
%
% Beam bending for nonconforming meshes, exact solution

% % conversion from plane stress to plane strain
% EX = EX/(1-VV^2);
% VV = VV/(1-VV);

SL=72;
Z1 = SL/D;
UEX= 2*Z1/EX*XX*YY;
UEY= -Z1/EX*(XX^2+VV*YY^2);

EE(1)= 2*Z1/EX*YY;
EE(2)= -2*Z1/EX*VV*YY;
EE(3)= 2*Z1/EX*XX;
EE(4)= -2*Z1/EX*XX;

% PE= 2*VV*Z1*YY/(1+VV);
% PEX=0;
% PEY= 2*VV*Z1/(1+VV);

ue(1) = UEX;
ue(2) = UEY;
% ue(3) = PE;


duex(1) = EE(1);
duex(2) = EE(4);
% duex(3) = PEX;

duey(1) = EE(3);
duey(2) = EE(2);
% duey(3) = PEY;