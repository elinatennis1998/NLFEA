function [ElemK,ElemF,hr] = NL_Elem2_2dF03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,iprob)
        
% Purpose: Compute stiffness matrix and residual vector for element
% Called by: FormFE.m
% Notes: The ordering of dofs in the stiffness/force are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        The internal force vector should be the negative of its proper
%        value, e.g. -B'*sigma, due to the definition that the
%        Residual = F_external - F_internal.
%        While nst = ndf*nen in general, typically the assembly process
%        only requires the values for i=1:ndf*nel, where nel is the number
%        of nodes on the current element.
%        For DG elements, the stiffness of the left and right sides may be
%        computed separately and then combined into the total stiffness
%        ElemK before exiting.
% Example: L_Elem3_2dVMS.m
        
        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);

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
              [Qxy, shgs, Jdet] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
                
            if nel == 3
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1)];
            elseif nel == 4
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0        Qxy(5,1)  0        Qxy(6,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 0         Qxy(5,2) 0         Qxy(6,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1) Qxy(5,2)  Qxy(5,1) Qxy(6,2)  Qxy(6,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1) Qxy(5,2) -Qxy(5,1) Qxy(6,2) -Qxy(6,1)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0        Qxy(5,1)  0        Qxy(6,1)  0        Qxy(7,1)  0        Qxy(8,1)  0        Qxy(9,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 0         Qxy(5,2) 0         Qxy(6,2) 0         Qxy(7,2) 0         Qxy(8,2) 0         Qxy(9,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1) Qxy(5,2)  Qxy(5,1) Qxy(6,2)  Qxy(6,1) Qxy(7,2)  Qxy(7,1) Qxy(8,2)  Qxy(8,1) Qxy(9,2)  Qxy(9,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1) Qxy(5,2) -Qxy(5,1) Qxy(6,2) -Qxy(6,1) Qxy(7,2) -Qxy(7,1) Qxy(8,2) -Qxy(8,1) Qxy(9,2) -Qxy(9,1)];
            end
                
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
                
            c1 = Wgt*Jdet*thick;

            [sigma3, cmat] = SigmaCmat2(F,JxX,mateprop,lam);
            
            Smat = [sigma3(1) 0  sigma3(3)/2 sigma3(3)/2
                    0 sigma3(2)  sigma3(3)/2 -sigma3(3)/2
                    sigma3(3)/2  sigma3(3)/2 (sigma3(2)+sigma3(1))/4 (sigma3(2)-sigma3(1))/4
                    sigma3(3)/2 -sigma3(3)/2 (sigma3(2)-sigma3(1))/4 (sigma3(2)+sigma3(1))/4]+cmat;
            
            ElemF = ElemF - c1*(Bmat'*(sigma3));
            
            ElemK = ElemK + c1*(Bmat'*Smat*Bmat);

        end %je
    ElemK;