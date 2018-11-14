% Tim Truster
% 10/06/2014
%
% Function to assemble the total Neumann boundary conditions from all
% types of sources:

% User loads
puserforce
Fext1 = lamda*Fc1+Fc1np+FcU; % standard contribution

if numSLmt > 0
    
    for load_entry = 1:numSLmt % loop over each of the sets of tabulated proportional multipliers for BCs
        multfact = SurfacesLmt{1,load_entry}; % multiplier table ID
        % get the current value of the multiplier for this Load table
        if multfact < 0 % has an extra row for the initial displacements
            if step == 0
                Loadmult = BCLoadFacts{multfact}(step+1);
            elseif numincrem1 > 1 % Interpolate between steps for current adapted increment
                Loadmult = BCLoadFacts{multfact}(step)+(BCLoadFacts{multfact}(step+1)-BCLoadFacts{multfact}(step))*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
            else
                Loadmult = BCLoadFacts{multfact}(step+1);
            end
        else
            if step == 0
                Loadmult = 0;
            else
              if numincrem1 > 1 % Interpolate between steps for current adapted increment
                if step > 1
                  Loadmult = BCLoadFacts{multfact}(step-1)+(BCLoadFacts{multfact}(step)-BCLoadFacts{multfact}(step-1))*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
                else
                  Loadmult = BCLoadFacts{multfact}(step)*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
                end
              else
                Loadmult = BCLoadFacts{multfact}(step);
              end
            end
        end
        Loadvals = Fc1mt{load_entry};
        Fext1 = Fext1 + Loadmult*Loadvals;
    end
    
end