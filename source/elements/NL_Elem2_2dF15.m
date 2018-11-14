function ElemF = NL_Elem2_2dF15(mateprop,bodyf,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,iprob)
        
% Purpose: Compute body force for element
% Called by: pbodyf.m
% Notes: The ordering of dofs in the external load vector are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        Body forces are treated as external loads, which can be scaled by
%        a proportional loading parameter.
%        ElemF is initialized as zero inside pbodyf.m and does not need to
%        be reinitialized here.
% Example: NL_Elem2_2dM.m
        
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            if exist('iprob','var') == 1 && iprob == 6
                lint = 91;
            else
                lint = 13;
            end
        elseif nel == 4
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;4;
            else
                lint = 4;
            end
        else
            if exist('iprob','var') == 1 && iprob == 6
                lint = 100;
            else
                lint = 9;
            end
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            if nel == 3
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0       
                    0         shl(1) 0         shl(2) 0         shl(3)];
            elseif nel == 4
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0        shl(4)  0       
                    0         shl(1) 0         shl(2) 0         shl(3) 0         shl(4)];
            elseif nel == 6
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0        shl(4)  0        shl(5)  0        shl(6)  0       
                    0         shl(1) 0         shl(2) 0         shl(3) 0         shl(4) 0         shl(5) 0         shl(6)];
            elseif nel == 9
%             % Form N matrix
            Nmat = [shl(1)  0        shl(2)  0        shl(3)  0        shl(4)  0        shl(5)  0        shl(6)  0        shl(7)  0        shl(8)  0        shl(9)  0       
                    0         shl(1) 0         shl(2) 0         shl(3) 0         shl(4) 0         shl(5) 0         shl(6) 0         shl(7) 0         shl(8) 0         shl(9)];
            end
                
            c1 = Wgt*Jdet*thick;
            
            if exist('iprob','var')
                if iprob == 6
                    mu = 40;
                    X = xl(1,1:nel)*shl;
                    Y = xl(2,1:nel)*shl;
                    fb = [ - (101*mu)/(101*X + 100) - (101*lam*(101*X + 100))/10000
                                                      0];
                else
                    fb = bodyf(1:2)';
                end
            else
                fb = bodyf(1:2)';
            end
            
            ElemF = ElemF + c1*Nmat'*fb;

        end %je
    ElemF;