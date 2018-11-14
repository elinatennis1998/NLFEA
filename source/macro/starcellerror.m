% Tim Truster
% 09/10/2013
%
% Library for getFbcell routines

	if nonlin == 0 %linear analysis
        switch iel
            case 3 %Stabilized Mixed Pressure-Displacement Element
                cellerrorL3
        end
    else %nonlinear analysis
        switch iel
            case 5 % Mixed form
                cellerrorNL5
        end
	end