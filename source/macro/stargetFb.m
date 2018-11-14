% Tim Truster
% 09/10/2013
%
% Library for getFbcell routines

	if nonlin == 0 %linear analysis
        switch iel
            case 3 %Stabilized Mixed Pressure-Displacement Element
                [eF] = getFbcellL3(xc,d(:,mat),xl,ul,ndfs,ndm,nele,nen,nummat,strong,MR,BR,MS,BS,sslote,tractyn);
        end
    else %nonlinear analysis
        switch iel
            case 5 % Mixed form
                [eF] = getFbcellNL5(xc,d(:,mat),xl,ul,ndfs,ndm,nele,nen,nummat,strong,MR,BR,MS,BS,sslote,tractyn,lamda);
        end
	end