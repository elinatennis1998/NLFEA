function [ue,duex,duey] = uexact_beam(XX,YY,EX,VV,PP1) %,PSPS

% %Plane Stress
% if PSPS == 's'
EX = EX/(1-VV^2);
VV = VV/(1-VV);
% end

EE = zeros(4,1);
ue = zeros(3,1);

C=1.d0;
	SL=10.D0;
	WI=2.D0*C^3/3.D0;
% 	PP1=2560.D0;
	Z1=PP1/(6.D0*EX*WI);
	UEX=-Z1*YY*((6.D0*SL-3.D0*XX)*XX+(2.D0+VV)*(YY*YY-C*C));
	UEY=Z1*(3.D0*VV*YY*YY*(SL-XX)+(4.D0+5.D0*VV)*C*C*XX+...
        (3.D0*SL-XX)*XX*XX);

	EE(1)=-Z1*YY*(6.D0*SL-6.D0*XX); %diff(UEX,x)
	EE(2)=Z1*6.D0*VV*YY*(SL-XX); %diff(UEy,y)
	EE(3)=-Z1*(6.D0*SL-3.D0*XX)*XX-Z1*(2.D0+VV)*(3.D0*YY^2-C^2); %diff(UEX,y)
	EE(4)=Z1*(-3.D0*VV*YY*YY+6.D0*SL*XX-3.D0*XX*XX+(4.D0+5.D0*VV)*C*C); %diff(UEY,x)

	PE=-PP1*(SL-XX)*YY*VV/(WI*(1.d0+VV));
	PEX=PP1*YY*VV/(WI*(1.d0+VV));
	PEY=-PP1*(SL-XX)*VV/(WI*(1.d0+VV));
	
     	ue(1) = UEX;
	ue(2) = UEY;
	ue(3) = PE;
	

      duex(1) = EE(1);
	duex(2) = EE(4);
      duex(3) = PEX;
	
      duey(1) = EE(3);
	duey(2) = EE(2);
	duey(3) = PEY;