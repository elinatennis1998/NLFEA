% Integrate stress and strain
isw = 52;
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
    case 50
        GrainSig = zeros(3,nummatCG);
        GrainEps = zeros(3,nummatCG);
        GrainVol = zeros(2,nummatCG);
        AssemQuant = 'AssemIntegSE';
    case 52
        GrainSig = zeros(36,nummatCG);
        GrainEps = zeros(3,nummatCG);
        GrainVol = zeros(2,nummatCG);
        AssemQuant = 'AssemIntegSE';
end

for grain = 1:nummatCG
    

    iel = MatTypeTable(2,grain); %iel   = ie(nie-1,ma); same thing;
    mateprop = MateT(grain,:);
    grainG1 = find(RegionOnElement==grain);
    numelemg = length(grainG1);
    
    for el = 1:numelemg
        elem = grainG1(el);
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

GrainXmid = GrainEps(1:ndm,:)./(ones(ndm,1)*GrainVol(1,:));
% centroids of grains
Coordinates(NodeTypeNum(2):ndm+1:NodeTypeNum(3)-1,:) = GrainXmid';
% GrainXcent = Coorsinates

% Taylor K, Sachs K
DmatTaylor = reshape(sum(GrainSig(1:9,1:nummatCG),2),3,3)/sum(GrainVol(1,:));
DmatSachs = zeros(3,3);
for grain = 1:nummatCG
    
    Dinv = GrainVol(1,grain)*inv(reshape(GrainSig(1:9,grain),3,3));
    DmatSachs = DmatSachs + Dinv;
    
end
DmatSachs = inv(DmatSachs)*nummatCG;

%Saving variables for MRDG
save('GrainEps','GrainEps');
save('GrainVol','GrainVol');