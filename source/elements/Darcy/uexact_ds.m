function [ue,duex,duey] = uexact_ds(XX,YY,kappamu,mu)

EE = zeros(4,1);
ue = zeros(3,1);

L=1.d0;
	SL=2.d0*pi/L;
    
	UEX=-SL*kappamu*cos(SL*XX)*sin(SL*YY);
	UEY=-SL*kappamu*sin(SL*XX)*cos(SL*YY);

	EE(1)=SL*SL*kappamu*sin(SL*XX)*sin(SL*YY); %diff(UEX,x)
	EE(2)=SL*SL*kappamu*sin(SL*XX)*sin(SL*YY); %diff(UEy,y)
	EE(3)=-SL*SL*kappamu*cos(SL*XX)*cos(SL*YY); %diff(UEX,y)
	EE(4)=-SL*SL*kappamu*cos(SL*XX)*cos(SL*YY); %diff(UEY,x)

	PE=(1+16*kappamu*mu*pi^2)*sin(SL*XX)*sin(SL*YY);
	PEX=(1+16*kappamu*mu*pi^2)*SL*cos(SL*XX)*sin(SL*YY);
	PEY=(1+16*kappamu*mu*pi^2)*SL*sin(SL*XX)*cos(SL*YY);
	
    ue(1) = UEX;
	ue(2) = UEY;
	ue(3) = PE;
	

    duex(1) = EE(1);
	duex(2) = EE(4);
    duex(3) = PEX;
	
    duey(1) = EE(3);
	duey(2) = EE(2);
	duey(3) = PEY;