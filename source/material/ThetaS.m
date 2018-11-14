function [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop)
%
% Volumetric constitutive functions

matv = mateprop(2);

if matv == 1

%           U(J) = 1/2*(J-1)^2;
    theta  = JxX-1;
    theta1 = 1;
    theta2 = 0;
    theta3 = 0;

else % matv == 2

%           U(J) = 1/2*(ln(J))^2;
    theta  = log(JxX);
    theta1 = 1/JxX;
    theta2 = -1/JxX^2;
    theta3 = 2/JxX^3;

end