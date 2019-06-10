% Tim Truster
% 6/28/17
%
% Form constant strains in each grain

% egG = [[.02 -0.02*.25 0];
%        [.02 -0.02*.25 0];
%        [.02 -0.02*.25 0];
%        [.02 -0.02*.25 0];
%        [.02 -0.02*.25 0];
%        [.02 -0.02*.25 0];
%        [.02 -0.02*.25 0];
%        [.02 -0.02*.25 0];
%        [.02 -0.02*.25 0]];
egG = 9*[[0.001999126542199  -0.000520253302762 0];
       [0.001999126542199  -0.000520253302762 0];
       [0.001999126542199  -0.000520253302762 0];
       [0.001999126542199  -0.000520253302762 0];
       [0.004006987662409  -0.000837973577902 0];
       [0.001999126542199  -0.000520253302762 0];
       [0.001999126542199  -0.000520253302762 0];
       [0.001999126542199  -0.000520253302762 0];
       [0.001999126542199  -0.000520253302762 0]];
em = [eRVE(1) eRVE(3)/2; eRVE(3)/2 eRVE(2)];
ugG = GrainXmid'*em;
ugG = [ugG zeros(9,1)];
   
DispVecD = zeros(neq,1);
DispVecF = zeros(nieq,1);
Node_U_VC = zeros(numnp,2);
   
for grain = 1:9
    
    eg = [egG(grain,1) egG(grain,3)/2; egG(grain,3)/2 egG(grain,2)];
    wg = [0 ugG(grain,3)/2; -ugG(grain,3)/2 0];
    bg = [ugG(grain,1); ugG(grain,2)];
    mg = eg + wg;
    xint = GrainXmid(1:2,grain);
    
    for el = 1:numelemg
        elem = grainG(grain,el);
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
        
        % Set displacement
        us = zeros(ndf, nen);
        for node = 1:nel
            xn = xl(:,node);
            us(:,node) = mg*(xn-xint) + bg;
            Node_U_VC(ElemFlag(node),:) = us(:,node)';
        end
        
        usreshape = reshape(us,nen*ndm,1);
        
        ElemFTemp = usreshape(ELDOFTa);
        DispVecD(EGDOFTa) = ElemFTemp;
        ElemFTemp = usreshape(ELDOFTi);
        DispVecF(EGDOFTi) = ElemFTemp;
        
    end
end
    