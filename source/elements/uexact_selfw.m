function [ue,duex,duey] = uexact_selfw(XX,YY,EX,VV,rho) %,PSPS

% %Plane Stress
% if PSPS == 's'
EX = EX/(1-VV^2);
VV = VV/(1-VV);
% end

EE = zeros(4,1);
ue = zeros(3,1);

D = 1.d0;
	L = 5.d0;
	grav = 9.81;
	XX1 = XX - L;
	XX2 = XX1*XX1;
	XX4 = XX2*XX2;
	YY2 = YY*YY;
	YY4 = YY2*YY2;

% 	VV2 = VV*VV;
	Z1 = rho*grav/(10.d0*EX*D*D);

	UEX = Z1*XX1*YY*(5.d0*XX2 + (5.d0*VV + 6.d0)*D*D - 15.d0*L*L - ...
          (10.d0 + 5.d0*VV)*YY2);
      
	UEY = Z1/4.d0*((5.d0 + 10.d0*VV)*YY4 - (10.d0 + 12.d0*VV)*YY2*D*D - ...
            30.d0*VV*YY2*XX2 + 30.d0*VV*YY2*L*L - 5.d0*XX4 + ...
            (48.d0 + 50.d0*VV)*XX2*D*D + 30.d0*XX2*L*L - 25.d0*L*L*L*L - ...
            (48.d0 + 50.d0*VV)*L*L*D*D);

	EE(1) = Z1*YY*(15.d0*XX2 - 15.d0*L*L + (5.d0*VV + 6.d0)*D*D - ...
            (5.d0*VV + 10.d0)*YY2);
	EE(2) = Z1*YY*(-15.d0*VV*XX2 + 15.d0*VV*L*L - (5.d0 + 6.d0*VV)*D*D + ...
            (10.d0*VV + 5.d0)*YY2);
	EE(3) = Z1*XX1*(5.d0*XX2 - 15.d0*L*L + (5.d0*VV + 6.d0)*D*D - ...
            15.d0*(VV + 2.d0)*YY2);                            
	EE(4) = Z1*XX1*(-5.d0*XX2 + 15.d0*L*L + (25.d0*VV + 24.d0)*D*D - ...
            15.d0*VV*YY2);
	
	Z1 = rho*grav*VV/(10.d0*(1.d0+VV)*D*D);
	
	PE = Z1*YY*(15.d0*(XX2 - L*L) + D*D - 5.d0*YY2);
	PEX = 30.d0*Z1*XX1*YY;
	PEY = Z1*(15.d0*(XX2 - L*L) + D*D - 15.d0*YY2);
	
     	ue(1) = UEX;
	ue(2) = UEY;
	ue(3) = PE;

      duex(1) = EE(1);
	duex(2) = EE(4);
      duex(3) = PEX;
	
      duey(1) = EE(3);
	duey(2) = EE(2);
	duey(3) = PEY;