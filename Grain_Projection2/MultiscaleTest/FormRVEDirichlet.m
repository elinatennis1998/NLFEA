% Tim Truster
% 06/30/2017
%
% Set up Dirichlet BC on boundary of RVE domain

NodeBC = zeros(ndm*(NodeTypeNum(4)-NodeTypeNum(3)),3);
numBC = 0;

for interI = MaterTypeNum(3):MaterTypeNum(4)-1
    
    jj = interI - MaterTypeNum(3);
    nodesI = NodeTypeNum(3)+(ndm+1)*jj:NodeTypeNum(3)-1+(ndm+1)*(jj+1);
    xmid_edge = BoundXmid(1:ndm,MaterTypeNum(3)-nummatCG+jj);
    Coordinates(nodesI(1),:) = xmid_edge';
    matLR = RegionsOnInterface(MaterTypeNum(3)-nummatCG+jj,2:3);
    if TFS == 3 || TFS == 4 %Sach strain correction
        DmatL = reshape(GrainSig(1:9,matLR(1)),3,3)/GrainVol(1,matLR(1));
        eG = factorTS*eRVE + (1-factorTS)*(inv(DmatL)*DmatSachs*eRVE')';
    else
        eG = eRVE;
    end
    xg = GrainXmid(1:ndm,matLR(1));
    xm = xmid_edge - xg;
    umid = uRVE + ...
          [eRVE(1)*xg(1)+0.5*eRVE(3)*xg(2) 0.5*eRVE(3)*xg(1)+eRVE(2)*xg(2)] + ...
          [eG(1)*xm(1)+0.5*eG(3)*xm(2) 0.5*eG(3)*xm(1)+eG(2)*xm(2)] + ...
          [0.5*wRVE*xmid_edge(2) -0.5*wRVE*xmid_edge(1)];
    if ndm == 2
    numBC = numBC + 1;
    NodeBC(numBC,:) = [nodesI(1) 1 umid(1)];
    numBC = numBC + 1;
    NodeBC(numBC,:) = [nodesI(1) 2 umid(2)];
    numBC = numBC + 1;
    NodeBC(numBC,:) = [nodesI(2) 1 eG(1)];
    numBC = numBC + 1;
    NodeBC(numBC,:) = [nodesI(2) 2 eG(2)];
    numBC = numBC + 1;
    NodeBC(numBC,:) = [nodesI(3) 1 eG(3)];
    numBC = numBC + 1;
    NodeBC(numBC,:) = [nodesI(3) 2 wRVE];
    else
        error('ndm = 3 not done yet')
    end
      
end