% 10/7/13
% Tim Truster
% UIUC
% Nonlinear large strain elasticity weak BC element.
% Based off of NL_Elem2_2dDG.m.
% Used to confirm quadratic convergence of general element routine.
% Debugged and verified 10/7/2013 to give quadratic convergence for
% patchtest5_2d.m with distorted loading, multiple load steps, ep = 20;

nitvms = 2;
if nitvms == 1 %VMS
pencoeff = 1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

h = 1/2;

switch isw %Task Switch
    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
%%
    case 3 %interface stiffness
        
        ElemKLL = zeros(nstL,nstL);
        ElemFL = zeros(nstL,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        bf1 = 0;
        bf2 = 0;
        
        PatchEL = matepropL(4);
        PatchvL = matepropL(5);
        muL = PatchEL/(2*(1+PatchvL));
        lamL = PatchvL*PatchEL/((1+PatchvL)*(1-2*PatchvL));
        I2 = eye(2);
        
        NmatL = zeros(2,nstL);
        BmatL = zeros(3,nstL);
        bnAdN1 = zeros(3,nstL);
        N1 = zeros(2,nstL);
        
        nelL = nel;
        ulL = ul;
        xlL = xl;
        ulresL = reshape(ulL(:,1:nelL),ndf*nelL,1);
        

        % Set jacobian for integration space
        drdr = 1;
        if nelL == 4 || nelL == 9
        eL1 = -1;
        roL = -1;
        else
        eL1 = 0;
        roL = 0;
        end
        
        % Perform integration of various matrices
        
        lint = 3;10;
        ib = 0;
        der = 0;
        bf = 0;
        
        if nitvms == 2
        % Nitsche
%         volL = getvol(xlL,nelL);
%         volR = getvol(xlR,nelR);
%         h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        ep = 20;100;
        end
        
        for ie = 1:lint
            
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
            
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nvect = [nLx 0 nLy
                     0 nLy nLx];
            nvec = [nLx; nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            for i = 1:nelL
                NmatL(1,(i-1)*ndf+1) = shlL(i);
                NmatL(2,(i-1)*ndf+2) = shlL(i);
                BmatL(Bcol1,(i-1)*ndf+1) = QxyL(i,col1);
                BmatL(Bcol2,(i-1)*ndf+2) = QxyL(i,col2);
            end
            
            [FL,JxXL,fiL] = kine2d(QxyL,ulL,nelL,0);
            PL = muL*(FL - fiL') + lamL*(JxXL^2-JxXL)*fiL';
                
            % Prescribed boundary condition
            gWBC = [0 0]';
            gvec = [1 0]'; % BC in x direction only
            BCmat = gvec*gvec'; % use eye(2) for BC in both x and y
%             bnAdN1 = gamL*nvect*cmatL*BmatL;
%             bnAdN2 = gamR*nvect*cmatR*BmatR;
            
            tvtr = PL*nvec; %average stress
            jumpu = gWBC - NmatL*ulresL; %displacement jump

            % Not needed for linear problems
            ElemFL = ElemFL - c1*( - NmatL'*(BCmat*tvtr + ep*BCmat*jumpu));%  + bnAdN1'*jumpu);% 
%             ElemFL = ElemFL - c1*( - NmatL'*(ep*BCmat*jumpu));%  + bnAdN1'*jumpu);% 

            % Symmetrizing term, "stresses" for weighting function
            WLi = 0;
            for Aw = 1:4
                WLi(1) = QxyL(Aw,1);
                WLi(2) = QxyL(Aw,2);
                for i = 1:2
                    for j = 1:2
                        termL = 0;
                        for I = 1:2
                            for J = 1:2
                                for k = 1:2
                                ALiIkJ = muL*(I2(i,k)*I2(I,J)+fiL(I,k)*fiL(J,i)) ...
                                       + lamL*((2*JxXL^2-JxXL)*fiL(I,i)*fiL(J,k) - (JxXL^2-JxXL)*fiL(I,k)*fiL(J,i));
                                termL = termL + c1*jumpu(j)*BCmat(j,k)*nvec(J)*ALiIkJ*WLi(I);
                                end
                            end
                        end
                        ElemFL(ndf*(Aw-1)+i) = ElemFL(ndf*(Aw-1)+i) - termL;
                    end
                end
            end
%             ElemFL = ElemFL - c1*( + bnAdN1'*jumpu);% );% 
%             ElemFR = ElemFR - c1*( + bnAdN2'*jumpu);% );% 

            % Non-symmetric term, enforcing traction continuity
            WLi = 0;
            dULj = 0;
            ElemKLL1 = zeros(nstL,nstL);
            for Aw = 1:4
                WLi = shlL(Aw);
                for Bu = 1:4
                    dULj(1) = QxyL(Bu,1);
                    dULj(2) = QxyL(Bu,2);
                    for i = 1:2
                        for j = 1:2
                            termLL = 0;
                            for I = 1:2
                                for J = 1:2
                                    for k = 1:2
                                    ALkIjJ = muL*(I2(k,j)*I2(I,J)+fiL(I,j)*fiL(J,k)) ...
                                           + lamL*((2*JxXL^2-JxXL)*fiL(I,k)*fiL(J,j) - (JxXL^2-JxXL)*fiL(I,j)*fiL(J,k));
                                    termLL = termLL - c1*WLi*BCmat(i,k)*nvec(I)*ALkIjJ*dULj(J);
                                    end
                                end
                            end
                            ElemKLL1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKLL1(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termLL;
                        end
                    end
                end
            end
%             ElemKLL = ElemKLL - c1*NmatL'*bnAdN1;
%             ElemKLR = ElemKLR - c1*NmatL'*bnAdN2;
%             ElemKRL = ElemKRL + c1*NmatR'*bnAdN1;
%             ElemKRR = ElemKRR + c1*NmatR'*bnAdN2;

            % Symmetrizing term, "stresses" for weighting function
            WLi = 0;
            dULj = 0;
            ElemKLL2 = zeros(nstL,nstL);
            for Aw = 1:4
                WLi(1) = QxyL(Aw,1);
                WLi(2) = QxyL(Aw,2);
                for Bu = 1:4
                    dULj = shlL(Bu);
                    for i = 1:2
                        for j = 1:2
                            termLL = 0;
                            for I = 1:2
                                for J = 1:2
                                    for k = 1:2
                                    ALiIkJ = muL*(I2(i,k)*I2(I,J)+fiL(I,k)*fiL(J,i)) ...
                                           + lamL*((2*JxXL^2-JxXL)*fiL(I,i)*fiL(J,k) - (JxXL^2-JxXL)*fiL(I,k)*fiL(J,i));
                                    termLL = termLL - c1*dULj*BCmat(j,k)*nvec(J)*ALiIkJ*WLi(I);
                                    end
                                end
                            end
                            ElemKLL2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKLL2(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termLL;
                        end
                    end
                end
            end
%             ElemKLL = ElemKLL - c1*bnAdN1'*NmatL;  % 
%             ElemKLR = ElemKLR + c1*bnAdN1'*NmatR;  % 
%             ElemKRL = ElemKRL - c1*bnAdN2'*NmatL;  %
%             ElemKRR = ElemKRR + c1*bnAdN2'*NmatR;  % 

            % Unexpected terms, involves jumpu and 6th order tensor
            WLi = 0;
            dULj = 0;
            ElemKLL3 = zeros(nstL,nstL);
            for Aw = 1:4
                WLi(1) = QxyL(Aw,1);
                WLi(2) = QxyL(Aw,2);
                for Bu = 1:4
                    dULj(1) = QxyL(Bu,1);
                    dULj(2) = QxyL(Bu,2);
                    for i = 1:2
                        for j = 1:2
                            termLL = 0;
                            for I = 1:2
                                for J = 1:2
                                    termDgLL = 0;
                                    for k = 1:2
                                        for K = 1:2
                                            for l = 1:2
DLiIjJlK = - muL*(fiL(I,l)*fiL(K,j)*fiL(J,i)+fiL(I,j)*fiL(J,l)*fiL(K,i)) ...
           + lamL*((4*JxXL^2-JxXL)*fiL(I,i)*fiL(J,j)*fiL(K,l) + (JxXL^2-JxXL)*(fiL(I,l)*fiL(K,j)*fiL(J,i) + fiL(I,j)*fiL(J,l)*fiL(K,i))...
           - (2*JxXL^2-JxXL)*(fiL(I,j)*fiL(J,i)*fiL(K,l) + fiL(I,l)*fiL(K,i)*fiL(J,j) + fiL(I,i)*fiL(J,l)*fiL(K,j)));
                                            termDgLL = termDgLL + DLiIjJlK*BCmat(k,l)*jumpu(k)*nvec(K);
                                            end
                                        end
                                    end
                                    termLL = termLL + c1*dULj(J)*termDgLL*WLi(I);
                                end
                            end
                            ElemKLL3(ndf*(Aw-1)+i,ndf*(Bu-1)+j) = ElemKLL3(ndf*(Aw-1)+i,ndf*(Bu-1)+j) + termLL;
                        end
                    end
                end
            end

            ElemKLL = ElemKLL + ElemKLL1 + ElemKLL2 + ElemKLL3 + c1*(NmatL'*ep*BCmat*NmatL);%
%             ElemKLL = ElemKLL + ElemKLL1 + c1*(NmatL'*ep*BCmat*NmatL);%
            
        end
% ElemKLL
            ElemK = ElemKLL;
            ElemF = ElemFL;
end
