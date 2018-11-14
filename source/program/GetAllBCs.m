% Tim Truster
% 10/06/2014
%
% Function to assemble the total Dirichlet boundary conditions from all
% types of sources:
% ModelDc = global proportional BCs
% ModelDcnp = global non-proportional BCs
% NodeBCmt = tabulated multipliers for regional BCs

gBC = lamda*ModelDc + ModelDcnp; % standard contribution

if numBCmt > 0
    
    for bc_entry = 1:numBCmt % loop over each of the sets of tabulated proportional multipliers for BCs
        multfact = NodeBCmt{1,bc_entry}; % multiplier table ID
        % get the current value of the multiplier for this BC table
        if multfact < 0 % has an extra row for the initial displacements
            if step == 0
                BCmult = BCLoadFacts{multfact}(step+1);
            elseif numincrem1 > 1 % Interpolate between steps for current adapted increment
                BCmult = BCLoadFacts{multfact}(step)+(BCLoadFacts{multfact}(step+1)-BCLoadFacts{multfact}(step))*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
            else
                BCmult = BCLoadFacts{multfact}(step+1);
            end
        else
            if step == 0
                BCmult = 0;
            else
              if numincrem1 > 1 % Interpolate between steps for current adapted increment
                if step > 1
                  BCmult = BCLoadFacts{multfact}(step-1)+(BCLoadFacts{multfact}(step)-BCLoadFacts{multfact}(step-1))*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
                else
                  BCmult = BCLoadFacts{multfact}(step)*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
                end
              else
                BCmult = BCLoadFacts{multfact}(step);
              end
            end
        end
        nodes = NodeBCmt{2,bc_entry}(:,1);
        dirs = NodeBCmt{2,bc_entry}(:,2);
        linearInd = sub2ind([numnp,ndf], nodes, dirs);
        dofs = NDOFT(linearInd) - neq;
        BCvals = NodeBCmt{2,bc_entry}(:,3);
        gBC(dofs) = gBC(dofs) + BCmult*BCvals;
    end
    
end