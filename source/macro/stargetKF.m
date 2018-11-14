% Tim Truster
% 09/10/2013
%
% Library for getKFcell routines

	if nonlin == 0 %linear analysis
        switch iel
            case 3 %Stabilized Mixed Pressure-Displacement Element
                [eK,eF] = getKFcellL3(xc,d(:,mat),xe,ue,ndfs,ndm,lint,nele,nen,nummat,strong,MR,BR,MS,BS,sslot,t11,t12,t21,t22);
        end
    else %nonlinear analysis
        switch iel
            case 5 % Mixed form
                [eK,eF] = getKFcellNL5(xc,d(:,mat),xe,ue,ndfs,ndm,lint,nele,nen,nummat,strong,MR,BR,MS,BS,sslot,t11,t12,t21,t22,symmns);
        end
	end