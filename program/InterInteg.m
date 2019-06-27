% Integrate stress and strain
isw = 51;
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
    case 51
        BoundSig = zeros(26,nummat-nummatCG);
        BoundEps = zeros(20,nummat-nummatCG);
        ElemS = zeros(26,1);
        ElemE = zeros(20,1);
        AssemQuant = 'AssemIntegTE';
end

for interI = 1:nummat-nummatCG
    
    ma = interI+nummatCG;
    matLR = RegionsOnInterface(interI,2:3);
    ellist = RegionsOnInterface(interI,5:6);
    iel = MatTypeTable(2,ma); %iel   = ie(nie-1,ma); same thing;
    mateprop = MateT(ma,:);
    sigLR = GrainSig(1:3,matLR);
    epsLR = GrainEps(1:3,matLR);
    volLR = GrainVol(1,matLR);
    disLR = GrainSig(4:5,matLR)./[volLR; volLR];
    
    for elem = ellist(1):ellist(2)
        
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

BoundSig = BoundSig./(ones(26,1)*BoundEps(1,:));