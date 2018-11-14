function [theta1,theta2,theta3] = ThetaNSIE(JxX,mateprop)
%
% Volumetric constitutive functions

matv = mateprop(2);

if matv <= 4
    if matv <= 2
        if matv == 1

%           U(J) = 1/2*(ln(J))^2;
            theta1 = log(JxX)/JxX;
            theta2 = (1-log(JxX))/JxX^2;
            theta3 = (2*log(JxX)-3)/JxX^3;
            
        else % matv == 2

%           U(J) = 1/2*(ln(J))^2;
            theta1 = log(JxX)/JxX;
            theta2 = (1-log(JxX))/JxX^2;
            theta3 = (2*log(JxX)-3)/JxX^3;

        end
    elseif matv == 3
        
%       U(J) = 1/2*(J-1)^2;
        theta1 = JxX-1;
        theta2 = 1;
        theta3 = 0;
        
    else % matv == 4

%           U(J) = 2*((J-1)-ln(J));
            theta1 = 1/2*2*(JxX-1)/JxX;
            theta2 = 1/2*2/JxX^2;
            theta3 = -1/2*4/JxX^3;
        
    end
else
    if matv <= 6
        if matv == 5

%           U(J) = 1/4*(ln(J))^2+1/4*(J-1)^2;
            theta1 = 1/2*log(JxX)/JxX+1/2*(JxX-1);
            theta2 = 1/2*(1-log(JxX))/JxX^2+1/2;
            theta3 = 1/2*(2*log(JxX)-3)/JxX^3;
            
        else % matv == 6
            
%           U(J) = 1/beta^2*(J^(-beta)-1+beta*ln(J));
            beta = mateprop(5);
            b2 = beta^2;
            theta1 = 1/b2*(-beta*JxX^(-beta-1)+beta/JxX);
            theta2 = 1/b2*(-beta*(-beta-1)*JxX^(-beta-2)-beta/JxX^2);
            theta3 = 1/b2*(-beta*(-beta-1)*(-beta-2)*JxX^(-beta-3)+beta*2/JxX^3);
            
        end
    elseif matv == 7
        
    else % matv == 8
        
    end
end