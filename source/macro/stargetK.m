% Tim Truster
% 09/10/2013
%
% Library for getKcell routines

	if nonlin == 0 %linear analysis
        switch iel
            case 3 %Stabilized Mixed Pressure-Displacement Element
                eK = getKcellL3(xc,d(:,mat),ndfs,ndm,lint,nele,nen,nummat);
        end
    else %nonlinear analysis
        switch iel
            case 5 % Mixed form
                eK = getKcellNL5(xc,d(:,mat),xe,ue,MR,BR,MS,BS,ndfs,ndm,lint,nele,nen,nummat,symmns);
        end
	end