%% A script to compute the mass matrix for R only once per element
if ~(exist('CountElem_omar','var'))
     CountElem_omar = 0; 
end
if ~(exist('Mass_omar','var'))
    Mass_omar = zeros(numnp,numnp,nummat_omar);
end
if ~(exist('DiagOne_omar','var'))
     DiagOne_omar = 1; 
end
if ~(exist('SparseIt_omar','var'))
     SparseIt_omar = 1; 
end

if ((CountElem_omar<numel) && (ma2 <= nummat_omar))
    Mass_omar(ElemFlag,ElemFlag,ma2) = Mass_omar(ElemFlag,ElemFlag,ma2)+ ElemM;
    CountElem_omar = CountElem_omar + 1;
    if ((CountElem_omar == numel) && (ma2 ~= nummat_omar))
        CountElem_omar = 0;
    end     
end
if ((CountElem_omar == numel) && (ma2 == nummat_omar) && (DiagOne_omar))
    for CountNodeC = 1:numnp
        for CountMatC = 1:nummat_omar
            if (Mass_omar(CountNodeC,CountNodeC,CountMatC) ==0)
                Mass_omar(CountNodeC,CountNodeC,CountMatC) =1;
            end
        end
    end
    DiagOne_omar = 0;
end

if ((CountElem_omar == numel) && (ma2 == nummat_omar) && (SparseIt_omar))
    MassCell_omar = cell(1,nummat_omar);
    for CountMatC = 1:nummat_omar
        MassCell_omar{CountMatC} = sparse(Mass_omar(:,:,CountMatC));
    end
    SparseIt_omar = 0;
end
