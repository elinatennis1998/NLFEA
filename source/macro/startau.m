% Tim Truster
% 09/10/2013
%
% extract tau matrix for current parent element of cells in submesh

	if nonlin == 0 %linear analysis
        switch iel
            case 3 %Stabilized Mixed Pressure-Displacement Element
                [t11,t12,t21,t22] = Tau3_2d(xe,mu,nele,nen,lint);
        end
    else %nonlinear analysis
        
        ht3 = 1 + ixFEAP(3,elem) + ieFEAP(nie-4,mat);

%             If Element variables exist move into nh3

        if(ieFEAP(nie-5,mat) > 0) %then
            for i = 0:ieFEAP(nie-5,mat)-1
              hr(nhc+i) = hrvec(ht3+i);
            end % i
        end
        
        switch iel
            case 5 % Mixed form
                t11 = hr(nhc+0);
                t12 = hr(nhc+1);
                t21 = hr(nhc+2);
                t22 = hr(nhc+3);
        end
	end