% Tim Truster
% 11/27/2013
%
% Three-field Lagrange multiplier tying element
% Created for two-constituent growth model
% Current version assumes 8-node brick, dofs are arranged as:
% [node1xyz1 node2xyz1 node3xyz1 node4xyz1 node1xyz2 node2xyz2 node3xyz2...
%  node4xyz2 lamd1xyz1 lamd2xyz1 lamd3xyz1 lamd4xyz1 lamd1xyz2 lamd2xyz2...
%  lamd3xyz2 lamd4xyz2 psi_1xyz  psi_2xyz  psi_3xyz  psi_4xyz]
% First number is the node number (1-4, which are the counter-clockwise
% nodes from one face of a brick) and the second is the constituent.
%
% Sign conventions are from Brezzi_GAMM_3.pdf


%% Task Switch - Perform various element functions
switch isw 
%%    
    case 1 % Setup up elemental degrees of freedom
        
        if ndf > ndm
            
            for i = ndm+1:ndf
                lie(i,1) = 0;
            end
            
        end


%%
    case 3 % Stiffness and internal force vector
        
        ElemF = zeros(nst,1);
        
        lint = 4;
        nelm = 4;
        Nmat = zeros(ndm,ndm*nelm);
        elemk = zeros(ndm*nelm);
        Zero = elemk;
        
        der = 0;
        bf = 0;
        ib = 0;
        
        for je = 1:lint

            [Wgt,litr,lits] = intpntq(je,lint,ib);
            [shl,shld,shls,be] = shlq(litr,lits,4,4,der,bf);
                    
            t1 = xl(:,1:nelm)*shld(:,1);%[xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = xl(:,1:nelm)*shld(:,2);%[xs(:,2); 0];
            [tm2, tu2] = VecNormalize(t2);
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);

            c1 = Wgt*tm3;
            
            % Form B matrix
            for ie = 1:nelm
                
                for i = 1:ndm
                    Nmat(i,(ie-1)*ndf+i) = shl(ie);
                end
                 
            end
            
            % Evaluate solution fields
            u1  = ul(:,     1  :  nelm)*shl;
            u2  = ul(:,nelm+1  :2*nelm)*shl;
            l1  = ul(:,2*nelm+1:3*nelm)*shl;
            l2  = ul(:,3*nelm+1:4*nelm)*shl;
            psi = ul(:,4*nelm+1:5*nelm)*shl;
            
            ElemF = ElemF - c1*[-Nmat'*l1
                                -Nmat'*l2
                                Nmat'*(psi-u1)
                                Nmat'*(psi-u2)
                                Nmat'*(l1+l2)];
            elemk = elemk + c1*(Nmat'*Nmat);
            
        end %je
        
        ElemK = [Zero   Zero -elemk   Zero  Zero
                 Zero   Zero   Zero -elemk  Zero
               -elemk   Zero   Zero   Zero elemk
                 Zero -elemk   Zero   Zero elemk
                 Zero   Zero  elemk  elemk  Zero];
        
%%
    case 6 % Internal force vector
        
        ElemF = zeros(nst,1);
        
        lint = 4;
        nelm = 4;
        Nmat = zeros(ndm,ndm*nelm);
        
        der = 0;
        bf = 0;
        ib = 0;
        
        for je = 1:lint

            [Wgt,litr,lits] = intpntq(je,lint,ib);
            [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                    
            t1 = xl(:,1:nelm)*shld(:,1);%[xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = xl(:,1:nelm)*shld(:,2);%[xs(:,2); 0];
            [tm2, tu2] = VecNormalize(t2);
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);

            c1 = Wgt*tm3;
            
            % Form B matrix
            for ie = 1:nelm
                
                for i = 1:ndm
                    Nmat(i,(ie-1)*ndf+i) = shl(ie);
                end
                 
            end
            
            % Evaluate solution fields
            u1  = ul(:,     1  :  nelm)*shl;
            u2  = ul(:,nelm+1  :2*nelm)*shl;
            l1  = ul(:,2*nelm+1:3*nelm)*shl;
            l2  = ul(:,3*nelm+1:4*nelm)*shl;
            psi = ul(:,4*nelm+1:5*nelm)*shl;
            
            ElemF = ElemF - c1*[-Nmat'*l1
                                -Nmat'*l2
                                Nmat'*(psi-u1)
                                Nmat'*(psi-u2)
                                Nmat'*(l1+l2)];
            
        end %je

        
end %Task Switch
