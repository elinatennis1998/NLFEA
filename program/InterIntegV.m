% Integrate stress and strain
isw = 53;
ndf = ndm;
psflags
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
        BoundSig = zeros(20,nummat-nummatCG);
        BoundEps = zeros(20,nummat-nummatCG);
        ElemS = zeros(20,1);
        ElemE = zeros(20,1);
        AssemQuant = 'AssemIntegTE';
    case 53
        BoundSig = zeros(1+18,nummat-nummatCG);
        BoundEps = zeros(3+18,nummat-nummatCG);
        ElemS = zeros(1+18,1);
        ElemE = zeros(3+18,1);
        AssemQuant = 'AssemIntegTE';
end

for interI = 1:nummat-nummatCG
    
    ma = interI+nummatCG;
    matLR = RegionsOnInterface(interI,2:3);
    ellist = RegionsOnInterface(interI,5:6);
    iel = MatTypeTable(2,ma); %iel   = ie(nie-1,ma); same thing;
    mateprop = MateT(ma,:);
    sigLR = GrainSig(:,matLR);
    epsLR = GrainEps(:,matLR);
    volLR = GrainVol(1,matLR);
    
    for elem = ellist(1):ellist(2)
        
        % Determine element size parameters
        nel = nnz(NodesOnElement(elem,1:nen));
        nst = nen*ndf;
        
        ElemFlag = NodesOnElement(elem,1:nen);
        actnode = find(ElemFlag>0);
        xl = zeros(ndm,nen);
        xl(1:ndm,actnode) = Coordinates(ElemFlag(actnode),1:ndm)';

        ElemRout
        
    evalin('caller',[AssemQuant ';']);
        
    end
end

BoundXmid = BoundEps(1:ndm,:)./(ones(ndm,1)*BoundSig(1,:));
% BoundSachs = BoundEps(4:21,:)./(ones(18,1)*BoundSig(1,:));
% BoundTay = BoundSig(2:19,:)./(ones(18,1)*BoundSig(1,:));