%FormKs
% Create sparsity pattern for stiffness matrix to conserve memory

% Determine number of linear and nonlinear materials
numberlinear = 0;
numbernonlinear = 0;
for ma = 1:nummat
    if MatTypeTable(3,ma) == 0
        numberlinear = numberlinear + 1;
    elseif MatTypeTable(3,ma) == 1
        numbernonlinear = numbernonlinear + 1;
    else
        disp('invalid definition of L/NL material flag')
        MatTypeTable(3,ma)
        return
    end
end

Ksdd1 = zeros(numel*(ndf*nen)^2,2);
Ksdf1 = zeros(numel*(ndf*nen)^2,2);
Ksff1 = zeros(numel*(ndf*nen)^2,2);
sinddd = 0;
sinddf = 0;
sindff = 0;

for elemn = 1:length(reglist1)
    
    elem = reglist1(elemn);
    
    %Determine element size parameters
    nel = nnz(NodesOnElement(elem,1:nen));
    nst = nel*ndf;
    
    %Extract element nodal coordinates
    ElemFlag = NodesOnElement(elem,1:nel);
    actnode = find(ElemFlag>0);

    %Compute and Assemble Patch Stiffness
%     EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT1, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq1);

    %Assemble Element contribution to Model Quantity

    AssemSparse1
    
end

Ksdd1 = sortrows(Ksdd1(1:sinddd,:),[1 2]);
Ksdd1 = unique(Ksdd1,'rows');
Ksdf1 = sortrows(Ksdf1(1:sinddf,:),[1 2]);
Ksdf1 = unique(Ksdf1,'rows');
Ksff1 = sortrows(Ksff1(1:sindff,:),[1 2]);
Ksff1 = unique(Ksff1,'rows');
