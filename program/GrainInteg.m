% Integrate stress and strain
isw = 50;
switch isw%Task Switch
    
    case 1 %Get Material Properties
        
    case 3 %Get Stiffness, Force
        Kdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
        Kdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
        Kfd = sparse(Ksdf(:,2),Ksdf(:,1),0,nieq,neq);
        Kff = sparse(Ksff(:,1),Ksff(:,2),0,nieq,nieq);
        AssemQuant = 'AssemStifForc';
        Fd1 = Fext1;
        Fd3 = zeros(nieq,1);
    case 6 %Get Force
        hflgu = 1;
        h3flgu = 1;
        AssemQuant = 'AssemForc';
        Fd1 = Fext1;
        Fd3 = zeros(nieq,1);
    case 21 %Get Stiffness
        Kdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
        Kdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
        Kfd = sparse(Ksdf(:,2),Ksdf(:,1),0,nieq,neq);
        Kff = sparse(Ksff(:,1),Ksff(:,2),0,nieq,nieq);
        AssemQuant = 'AssemStif';
    case 50
        GrainSig = zeros(3+2,nummatCG);
        GrainEps = zeros(3,nummatCG);
        GrainVol = zeros(2,nummatCG);
        GrainDisp = zeros(2,nummatCG);
        AssemQuant = 'AssemIntegSE';
%         AssemQuant = 'AssemIntegSE_post';
end

for grain = 1:nummatCG
    

    iel = MatTypeTable(2,grain); %iel   = ie(nie-1,ma); same thing;
    mateprop = MateT(grain,:);
    grainG1 = find(RegionOnElement==grain);
    numelemg = length(grainG1);
    
    for el = 1:numelemg
        elem = grainG1(el);
        % Determine element size parameters
        nel = nnz(ix(elem,1:nen));
        nst = nen*ndf;
        
        ElemFlag = ix(elem,1:nen);
        actnode = find(ElemFlag>0);
        xl = zeros(ndm,nen);
        xl(1:ndm,actnode) = NodeTable(ElemFlag(actnode),1:ndm)';

        % Determine global to local equation numbers
        [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq);
        
        ul = zeros(ndf, nen);
        ul(ELDOFTa) = ModelDx(EGDOFTa)';
        ul(ELDOFTi) = gBC(EGDOFTi)';

        ElemRout
        
    evalin('caller',[AssemQuant ';']);
        
    end
end

GrainSig;
GrainEps;
GrainDisp;
