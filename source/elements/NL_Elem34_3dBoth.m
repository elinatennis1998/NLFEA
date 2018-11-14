% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 07/14/2015
% UIUC

% Combo element for mixture theory
% Needs to have both solid elements defined together
% mateprop is a combo of each constituent, with
% [iel1 props1(1:13) iel2 props2(1:13)]

% NOTE: assumes only one growing constituent; otherwise, need to expand hr
% array.

% Default type of interactive force
if ~exist('InFoType','var')
    InFoType = 1;
end

switch isw %Task Switch
%%
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = 1*16;
%%
    case {3,6}
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat1 = zeros(9,nel*ndf/2);
        Nmat1 = zeros(3,nel*ndf/2);
        Bmat2 = zeros(9,nel*ndf/2);
        Nmat2 = zeros(3,nel*ndf/2);

        iel1 = mateprop(1);
        mateprop1 = mateprop(2:14);
        iel2 = mateprop(15);
        mateprop2 = mateprop(16:28);
        nel1 = nel/2;
        nel2 = nel/2;
        xl1 = xl(:,1:nel/2);
        ul1 = ul(:,1:nel/2);
        xl2 = xl(:,nel/2+1:nel);
        ul2 = ul(:,nel/2+1:nel);
        % For testing the calculations of dF and dN
%         rng('default')
%         ul1 = rand([3,8]);
%         ul2 = rand([3,8]);
%         P = rand(3);
        
        % Load Guass Integration Points

        if nel/2 == 4
            lint = 16;11;5;16;
        elseif nel/2 == 8
            lint = 8;
        elseif nel/2 == 10
            lint = 14;
%             lint = 27;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 0;

        
        if InFoType == 1
        % Compute Grad(F)
%         % Set default values for dF; needs to be pre-computed
%         dF1 = zeros(3,9);
%         dF2 = zeros(3,9);
        Fall1 = zeros(9,lint);
        Fall2 = zeros(9,lint);
        Ball1 = zeros(ndm*nel1,lint);
        Ball2 = zeros(ndm*nel2,lint);
        intermat = zeros(lint,4);
        xmid = mean(xl1,2);
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel1 == 4 || nel1 == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel1,nel1,der,bf);
                  [QXY1, shgs1, Jdet1] = shgtt(xl1(:,1:nel),nel1,shld,shls,nel1,bf,der,be);
                  [QXY2, shgs2, Jdet2] = shgtt(xl2(:,1:nel),nel2,shld,shls,nel2,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel1,nel1,der,bf);
                  [QXY1, shgs1, Jdet1] = shgb(xl1(:,1:nel1),nel1,shld,shls,nel1,bf,der,be);
                  [QXY2, shgs2, Jdet2] = shgb(xl2(:,1:nel2),nel2,shld,shls,nel2,bf,der,be);
                end
                
            F1 = ul1(1:3,1:nel1)*QXY1;%+eye(3);
            Fall1(1:9,ll) = reshape(F1,9,1);
            Ball1(:,ll) = reshape(QXY1',ndm*nel1,1);
            F2 = ul2(1:3,1:nel2)*QXY2;%+eye(3);
            Fall2(1:9,ll) = reshape(F2,9,1);
            Ball2(:,ll) = reshape(QXY2',ndm*nel2,1);
            
            intermat(ll,1:3) = xl1*shl - xmid;
%             intermat(ll,1:3) = ss(1:3)';
            intermat(ll,4) = 1.0;
            
        end
        
        % Perform least squares fit to a linear surface
        LHS1 = intermat\Fall1';
        LHS2 = intermat\Fall2';
        % Drop the constant coefficient and store
        dF1 = LHS1(1:3,1:9);
        dF2 = LHS2(1:3,1:9);
        % Note: F and P are reshaped the same way, such that dFvec'*Pvec =
        % dF11*P11 + dF12*P12 + ...
        % Entries are:
        %[ F1xx,x F1yx,x F1zx,x  F1xy,x F1yy,x F1zy,x  F1xz,x F1yz,x F1zz,x
        %[ F1xx,y F1yx,y F1zx,y  F1xy,y F1yy,y F1zy,y  F1xz,y F1yz,y F1zz,y
        %[ F1xx,z F1yx,z F1zx,z  F1xy,z F1yy,z F1zy,z  F1xz,z F1yz,z F1zz,z
        
        % Derivative terms for stiffness matrix
        LHS1 = intermat\Ball1';
        LHS2 = intermat\Ball2';
        % Drop the constant coefficient and store
        dN1 = LHS1(1:3,1:ndm*nel1);
        dN2 = LHS2(1:3,1:ndm*nel2);
        % Entries are (before reorder):
        %[ dN1x,x dN1y,x dN1z,x  dN2x,x dN2y,x dN2z,x ...
        %[ dN1x,y dN1y,y dN1z,y  dN2x,y dN2y,y dN2z,y ...
        %[ dN1x,z dN1y,z dN1z,z  dN2x,z dN2y,z dN2z,z ...
        % verify correct ordering
%         dN1(2,1:3:1+3*7)*ul1(1,1:8)', dF1(2,1) %F1xx,y
%         dN1(3,2:3:2+3*7)*ul1(1,1:8)', dF1(3,4) %F1xy,z
%         dN1(1,2:3:2+3*7)*ul1(2,1:8)', dF1(1,5) %F1yy,x
%         i = 3; j = 3; k = 2;
%         dN1(k,j:3:j+3*7)*ul1(i,1:8)', dF1(k,i+3*(j-1)) %F1yy,x
        % Reorder for better matrix multiply
%         for i = 1:nel1
%             dN1(1:3,3*i-2:3*i) = dN1(1:3,3*i-2:3*i)';
%         end
        dN1b = reshape(dN1,9,8); % This is a little faster
        dN1b = dN1b([1 4 7 2 5 8 3 6 9],:);
        dN1 = reshape(dN1b,3,24);
%         for i = 1:nel2
%             dN2(1:3,3*i-2:3*i) = dN2(1:3,3*i-2:3*i)';
%         end
        dN1b = reshape(dN2,9,8);
        dN1b = dN1b([1 4 7 2 5 8 3 6 9],:);
        dN2 = reshape(dN1b,3,24);
%         % verify correct ordering
%         i = 1; j = 3; k = 2;
%         dN1(j,k:3:k+3*7)*ul1(i,1:8)', dF1(k,i+3*(j-1)) %F1yy,x
%         Pvec = reshape(P,9,1);
%         dF1*Pvec
%         PN = P*dN1;
%         for i = 1:nel1
%             PN(1:3,3*i-2:3*i) = PN(1:3,3*i-2:3*i)';
%         end
%         PN*reshape(ul1,24,1)
        end
        
        
        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel1 == 4 || nel1 == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel1,nel1,der,bf);
                  [QXY1, shgs1, Jdet1] = shgtt(xl1(:,1:nel),nel1,shld,shls,nel1,bf,der,be);
                  [QXY2, shgs2, Jdet2] = shgtt(xl2(:,1:nel),nel2,shld,shls,nel2,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel1,nel1,der,bf);
                  [QXY1, shgs1, Jdet1] = shgb(xl1(:,1:nel1),nel1,shld,shls,nel1,bf,der,be);
                  [QXY2, shgs2, Jdet2] = shgb(xl2(:,1:nel2),nel2,shld,shls,nel2,bf,der,be);
                end
                
            for mm = 1:nel1
 Nmat1(:,3*mm-2:3*mm) = [shl(mm) 0         0         
                        0         shl(mm) 0         
                        0         0         shl(mm) ];
            end
            Nmat2 = Nmat1;
                
            [F1,JxX1,fi1] = kine3d(QXY1,ul1(:,1:nel1),nel1,0);
                
            c1 = Wgt*Jdet1;
        
            % Form B matrix; will use 1st PK stress for all materials
            for mm = 1:nel1
 Bmat1(:,3*mm-2:3*mm) = [QXY1(mm,1) 0         0         
                         0         QXY1(mm,1) 0         
                         0         0         QXY1(mm,1) 
                         QXY1(mm,2) 0         0         
                         0         QXY1(mm,2) 0
                         0         0         QXY1(mm,2) 
                         QXY1(mm,3) 0         0         
                         0         QXY1(mm,3) 0
                         0         0         QXY1(mm,3) ];
            end
            
            [F2,JxX2,fi2] = kine3d(QXY2,ul2(:,1:nel2),nel2,1);
                
            c2 = Wgt*Jdet2;
            
            % Form B matrix; will use 1st PK stress for all materials
            for mm = 1:nel2
 Bmat2(:,3*mm-2:3*mm) = [QXY2(mm,1) 0         0         
                         0         QXY2(mm,1) 0         
                         0         0         QXY2(mm,1) 
                         QXY2(mm,2) 0         0         
                         0         QXY2(mm,2) 0
                         0         0         QXY2(mm,2) 
                         QXY2(mm,3) 0         0         
                         0         QXY2(mm,3) 0
                         0         0         QXY2(mm,3) ];
            end
            
            if iel1 == 33
                lam1 = getlam(mateprop1);
                [sigma1, cmat1] = SigmaCmat3(F1,JxX1,mateprop1,lam1);
                P = tranr4(fi1,fi1); % transformation tensor F_iI*F_jJ
                P = [P' zeros(6,3); zeros(3,9)];
                Se = P*sigma1; % second PK stress
                Se_m = [Se(1) Se(4) Se(6); Se(4) Se(2) Se(5); Se(6) Se(5) Se(3)];
                Cmat_e = P*cmat1*P';
                Pe1 = F1*Se_m; % first PK stress % (5)
                Ae1 = CSFtoA(Cmat_e,Se_m,F1,ndm);
                Aeg_mat1 = reshape(Ae1,ndm*ndm,ndm*ndm);
                Pe_vec1 = reshape(Pe1,ndm*ndm,1);
            elseif iel1 == 26 
                i_var = hr(nha-1+ll);
                [Aeg1,Pe1,var] = cnst_den(F1,i_var,mateprop1,ndm); % density
                hr(nhb-1+ll) = var;

                Aeg_mat1 = reshape(Aeg1,ndm*ndm,ndm*ndm);
                Pe_vec1 = reshape(Pe1,ndm*ndm,1);
            elseif iel1 == 32
                i_var = hr(nha-1+ll);
                [var,Pe1,Aeg1] = updt_cyl(F1,i_var,mateprop1,ndm); % cylindrical
    %             [Aeg1,Pe1,var]=cnst_vol(F1,i_var,mateprop1,ndm); % isotropic
                hr(nhb-1+ll) = var;

                Aeg_mat1 = reshape(Aeg1,ndm*ndm,ndm*ndm);
                Pe_vec1 = reshape(Pe1,ndm*ndm,1);
            end
            
            if iel2 == 33
                lam2 = getlam(mateprop2);
                [sigma2, cmat2] = SigmaCmat3(F2,JxX2,mateprop2,lam2);
                P = tranr4(fi2,fi2); % transformation tensor F_iI*F_jJ
                P = [P' zeros(6,3); zeros(3,9)];
                Se = P*sigma2; % second PK stress
                Se_m = [Se(1) Se(4) Se(6); Se(4) Se(2) Se(5); Se(6) Se(5) Se(3)];
                Cmat_e = P*cmat2*P';
                Pe2 = F2*Se_m; % first PK stress % (5)
                Ae2 = CSFtoA(Cmat_e,Se_m,F2,ndm);
                Aeg_mat2 = reshape(Ae2,ndm*ndm,ndm*ndm);
                Pe_vec2 = reshape(Pe2,ndm*ndm,1);
            elseif iel2 == 26 
                i_var = hr(nha-1+ll+8);
                [Aeg2,Pe2,var] = cnst_den(F2,i_var,mateprop2,ndm); % density
                hr(nhb-1+ll+8) = var;

                Aeg_mat2 = reshape(Aeg2,ndm*ndm,ndm*ndm);
                Pe_vec2 = reshape(Pe2,ndm*ndm,1);
            elseif iel2 == 32
                i_var = hr(nha-1+ll+8);
                [var,Pe2,Aeg2] = updt_cyl(F2,i_var,mateprop2,ndm); % cylindrical
    %             [Aeg2,Pe2,var]=cnst_vol(F2,i_var,mateprop2,ndm); % isotropic
                hr(nhb-1+ll+8) = var;

                Aeg_mat2 = reshape(Aeg2,ndm*ndm,ndm*ndm);
                Pe_vec2 = reshape(Pe2,ndm*ndm,1);
            end
                
            % Interactive force
            
            InterForce4
            
            % Stress terms
            
            ElemF(1:nel1*ndf) = ElemF(1:nel1*ndf) - c1*(Bmat1'*(Pe_vec1));
            
            ElemK(1:nel1*ndf,1:nel1*ndf) = ElemK(1:nel1*ndf,1:nel1*ndf) + c1*(Bmat1'*Aeg_mat1*Bmat1);
            
            
            ElemF(nel1*ndf+1:nel*ndf) = ElemF(nel1*ndf+1:nel*ndf) - c2*(Bmat2'*(Pe_vec2));
            
            ElemK(nel1*ndf+1:nel*ndf,nel1*ndf+1:nel*ndf) = ElemK(nel1*ndf+1:nel*ndf,nel1*ndf+1:nel*ndf) + c2*(Bmat2'*Aeg_mat2*Bmat2);
            
            if ll == 1 && step > 0% Only do for the first integration point
            IForceList(1:3,elem,step) = IntFor;
            end 
            
        end %je
   ElemK;     
%%
    case -1

        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        nel = nel/2;
        SurfOrientEtoS
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        
        if nel == 4 || nel == 10
            lint = 13;
        else
            lint = 16;
        end
        der = 0;
        bf = 0;
        ib = 0;
        
        % Integration Loop
        for je = 1:lint

            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(je,lint,edge);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgtt(xl(:,1:nel),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(je,lint,5);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgb(xl(:,1:nel),nel,shld,shls,nel,bf,der,be);
            end
          
            %Evaluate tangent and normal vectors
            t1 = sx(:,2);
            [tm1, tu1] = VecNormalize(t1);
            t2 = sx(:,1);
            [tm2, tu2] = VecNormalize(t2);
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);

            %Force components are positive in positive coord. direction
            c1 = Wgt*tm3;
            for o=1:nel
                don = shl(o);
                F = don*traction';

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(2)*c1;

                ElemF(ndf*o-0) = ElemF(ndf*o-0) + F(3)*c1;

            end %o
            
        end %je
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
        nel = nel*2;
%%  
    case 22
        
        ElemM = zeros(nst);
        ElemF = zeros(nst,1);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        spvec0 = I1;
        
        % Load Guass Integration Points

        if nel == 4
            lint = 11; %minimum of 7 for all integrals in deformed state
            der = 0;
        elseif nel == 8
%             lint = 4;
            lint = 8;
            der = 1;
        elseif nel == 10
            lint = 14; %minimum of 13 for all integrals in deformed state
            der = 1;
        else
%             lint = 9;
            lint = 27;
            der = 0;
        end

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgtt(xl(:,1:nel)+ul(:,1:nel),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy, shgs, Jdet] = shgb(xl(:,1:nel)+ul(:,1:nel),nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 4 || nelP == 10
              [shlS,shld,shls] = shltt(ss,nelS,nel,0,0);
            else
              [shlS,shld,shls] = shlb(ss,nelS,nel,0,0);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul(:,1:nel),nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            sigma = SigmaCmat3(F,JxX,mateprop,lam);
            sigma = sigma/JxX;
            
            % Form B matrix
            Nmat = shlS';

            w = Wgt*Jdet;
            
            if stres <= 6 % stress components
                sigmas = sigma(stres);
            elseif stres >= 8
                if stres <= 10 % principal stresses
                sigma = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma);
                sigmas = psig(stres-7);
                else % hydrostatic stress
                sigmas = 1/3*sigma'*I1;
                end
            else % von Mises stress
                trs = sigma'*I1;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemF = ElemF + w*Nmat'*sigmas;
            
            ElemM = ElemM + w*(Nmat'*Nmat);

        end %je
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr*2);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        spvec0 = I1;
        

        iel1 = mateprop(1);
        mateprop1 = mateprop(2:14);
        iel2 = mateprop(15);
        mateprop2 = mateprop(16:28);
        nel1 = nel/2;
        nel2 = nel/2;
        xl1 = xl(:,1:nel/2);
        ul1 = ul(:,1:nel/2);
        xl2 = xl(:,nel/2+1:nel);
        ul2 = ul(:,nel/2+1:nel);
%         lam = getlam(mateprop);
%         kIF = mateprop(8);
        
        % Load Guass Integration Points

        if nel/2 == 4
            nint = 1;
        elseif nel/2 == 8
%             lint = 4;
            lint = 8;
            nint = 1;
        elseif nel/2 == 10
            lint = 11;
            nint = 4;
        else
            lint = 27;
            nint = 8;
        end
        der = 1;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel1 == 4 || nel1 == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel1,nel1,der,bf);
              [Qxy1, shgs1, Jdet1] = shgtt(xl1(:,1:nel)+ul1(:,1:nel),nel1,shld,shls,nel1,bf,der,be);
              [Qxy2, shgs2, Jdet2] = shgtt(xl2(:,1:nel)+ul2(:,1:nel),nel2,shld,shls,nel2,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel1,nel1,der,bf);
              [Qxy1, shgs1, Jdet1] = shgb(xl1(:,1:nel1)+ul1(:,1:nel1),nel1,shld,shls,nel1,bf,der,be);
              [Qxy2, shgs2, Jdet2] = shgb(xl2(:,1:nel2)+ul2(:,1:nel2),nel2,shld,shls,nel2,bf,der,be);
            end
            
            
            [fi1,JxX1,F1,QXY1] = kine3d(Qxy1,-ul1(:,1:nel1),nel1,1); %this is equivalent to ikine2d
            JxX1 = 1/JxX1; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet1 = Jdet1/JxX1;
                
            c1 = Wgt*Jdet1;
            
            if iel1 == 33
            for mm = 1:nel1   
 Bmat1(:,3*mm-2:3*mm) = [Qxy1(mm,1) 0         0         
                        0         Qxy1(mm,2) 0         
                        0         0         Qxy1(mm,3) 
                        Qxy1(mm,2) Qxy1(mm,1) 0         
                        0         Qxy1(mm,3) Qxy1(mm,2) 
                        Qxy1(mm,3) 0         Qxy1(mm,1) 
                        Qxy1(mm,2) -Qxy1(mm,1) 0         
                        0         Qxy1(mm,3) -Qxy1(mm,2) 
                        -Qxy1(mm,3) 0         Qxy1(mm,1) ];
            end
            elseif iel1 == 26 || iel1 == 32
            for mm = 1:nel1
 Bmat1(:,3*mm-2:3*mm) = [QXY1(mm,1) 0         0         
                         0         QXY1(mm,1) 0         
                         0         0         QXY1(mm,1) 
                         QXY1(mm,2) 0         0         
                         0         QXY1(mm,2) 0
                         0         0         QXY1(mm,2) 
                         QXY1(mm,3) 0         0         
                         0         QXY1(mm,3) 0
                         0         0         QXY1(mm,3) ];
            end
            end
            
            [fi2,JxX2,F2,QXY2] = kine3d(Qxy2,-ul2(:,1:nel2),nel2,1); %this is equivalent to ikine2d
            JxX2 = 1/JxX2; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet2 = Jdet2/JxX2;
                
            c2 = Wgt*Jdet2;
            
            if iel2 == 33
            for mm = 1:nel2
 Bmat2(:,3*mm-2:3*mm) = [Qxy2(mm,1) 0         0         
                        0         Qxy2(mm,2) 0         
                        0         0         Qxy2(mm,3) 
                        Qxy2(mm,2) Qxy2(mm,1) 0         
                        0         Qxy2(mm,3) Qxy2(mm,2) 
                        Qxy2(mm,3) 0         Qxy2(mm,1) 
                        Qxy2(mm,2) -Qxy2(mm,1) 0         
                        0         Qxy2(mm,3) -Qxy2(mm,2) 
                        -Qxy2(mm,3) 0         Qxy2(mm,1) ];
            end
            elseif iel2 == 26 || iel2 == 32
            for mm = 1:nel2
 Bmat2(:,3*mm-2:3*mm) = [QXY2(mm,1) 0         0         
                         0         QXY2(mm,1) 0         
                         0         0         QXY2(mm,1) 
                         QXY2(mm,2) 0         0         
                         0         QXY2(mm,2) 0
                         0         0         QXY2(mm,2) 
                         QXY2(mm,3) 0         0         
                         0         QXY2(mm,3) 0
                         0         0         QXY2(mm,3) ];
            end
            end
            
            if iel1 == 33
                lam1 = getlam(mateprop1);
                [sigma1, cmat1] = SigmaCmat3(F1,JxX1,mateprop1,lam1);
                sigma1 = sigma1/JxX1;
                densit1 = mateprop(29);
            elseif iel1 == 26 
                i_var = hr(nha-1+ll);
                [Aeg1,Pe1,var] = cnst_den(F1,i_var,mateprop1,ndm); % density
                sigma1 = Pe1*F1'/JxX1;
                sigma1 = [sigma1(1,1); sigma1(2,2); sigma1(3,3); sigma1(1,2); sigma1(2,3); sigma1(3,1); zeros(3,1)];
                densit1 = (1 + var)*mateprop1(3);
            elseif iel1 == 32
                i_var = hr(nha-1+ll);
                [var,Pe1,Aeg1] = updt_cyl(F1,i_var,mateprop1,ndm); % cylindrical
    %             [Aeg1,Pe1,var]=cnst_vol(F1,i_var,mateprop1,ndm); % isotropic
                mo  = mateprop1(10:12); %fiber orientation vector
                [~, mo] = VecNormalize(mo); % ensure mo is a unit vector
                mo = mo';
                delta  = eye(ndm); 
                theta  = var(1) + 1;
    %             Fg_inv = (1/theta)*delta;
                Fg_inv = (1/theta)*delta + (1-1/theta)*(mo*mo');
                Fe = F1 * Fg_inv;
                Fe_inv = inv(Fe);
                Je = det(Fe);
                sigma1 = Pe1*Fe'/Je;
                sigma1 = [sigma1(1,1); sigma1(2,2); sigma1(3,3); sigma1(1,2); sigma1(2,3); sigma1(3,1); zeros(3,1)];
            end
            
            if iel2 == 33
                lam2 = getlam(mateprop2);
                [sigma2, cmat2] = SigmaCmat3(F2,JxX2,mateprop2,lam2);
                sigma2 = sigma2/JxX2;
                densit2 = mateprop(30);
            elseif iel2 == 26 
                i_var = hr(nha-1+ll+8);
                [Aeg2,Pe2,var] = cnst_den(F2,i_var,mateprop2,ndm); % density
                sigma2 = Pe2*F2'/JxX2;
                sigma2 = [sigma2(1,1); sigma2(2,2); sigma2(3,3); sigma2(1,2); sigma2(2,3); sigma2(3,1); zeros(3,1)];
                densit2 = (1 + var)*mateprop2(3);
            elseif iel2 == 32
                i_var = hr(nha-1+ll+8);
                [var,Pe2,Aeg2] = updt_cyl(F2,i_var,mateprop2,ndm); % cylindrical
                mo  = mateprop2(10:12); %fiber orientation vector
                [~, mo] = VecNormalize(mo); % ensure mo is a unit vector
                mo = mo';
                delta  = eye(ndm); 
                theta  = var(1) + 1;
    %             Fg_inv = (1/theta)*delta;
                Fg_inv = (1/theta)*delta + (1-1/theta)*(mo*mo');
                Fe = F2 * Fg_inv;
                Fe_inv = inv(Fe);
                Je = det(Fe);
                sigma2 = Pe2*Fe'/Je;
                sigma2 = [sigma2(1,1); sigma2(2,2); sigma2(3,3); sigma2(1,2); sigma2(2,3); sigma2(3,1); zeros(3,1)];
            end
            
            for stres = 1:npstr
            
            if stres <= 6 % stress components
                sigmas = sigma1(stres);
            elseif stres >= 7
                if stres <= 9 % Interactive force
                sigmas = IForceList(stres-6,elem,step);
                elseif stres == 10 % hydrostatic stress
                sigmas = 1/3*sigma1'*I1;
                elseif stres == 11 % internal variable
                    if iel1 == 26
                        sigmas = densit1;
                    elseif iel1 == 32
                        sigmas = theta;
                    else
                        sigmas = 0;
                    end
                end
            end
            
            ElemS2(ll,stres) = sigmas;
            
            if stres <= 6 % stress components
                sigmas = sigma2(stres);
            elseif stres >= 7
                if stres <= 9 % Interactive force
                sigmas = IForceList(stres-6,elem,step);
                elseif stres == 10 % hydrostatic stress
                sigmas = 1/3*sigma2'*I1;
                elseif stres == 11 % internal variable
                    if iel2 == 26
                        sigmas = densit2;
                    elseif iel2 == 32
                        sigmas = theta;
                    else
                        sigmas = 0;
                    end
                end
            end
            
            ElemS2(ll,stres+11) = sigmas;
            
            end

        end %je
        
        % interpolate stress at nodes
        if nel1 == 4
            plist = [1 0 0 0
                     0 1 0 0
                     0 0 0 1];
        elseif nel1 == 8
            plist = [-1 1 1 -1 -1 1 1 -1
                     -1 -1 1 1 -1 -1 1 1
                     -1 -1 -1 -1 1 1 1 1];
        elseif nel1 == 10
            plist = [ 1.927050983124842 -0.309016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947 -0.309016994374947  0.809016994374947  0.809016994374947 -0.309016994374947 -0.309016994374947
                     -0.309016994374947  1.927050983124842 -0.309016994374947 -0.309016994374947  0.809016994374947  0.809016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947 -0.309016994374947
                     -0.309016994374947 -0.309016994374947 -0.309016994374947  1.927050983124842 -0.309016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947  0.809016994374947  0.809016994374947];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 0 -sqr3 sqr3 0 0 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0
                     -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0 0 0 0 0];
        end
        
        for ll = 1:nel1
            
            r = plist(1,ll);
            s = plist(2,ll);
            t = plist(3,ll);
            shpS = sshp3d(r,s,t,nint);
            
            for stres = 1:npstr
                
                sigmas = ElemS2(1:nint,stres)'*shpS;
                ElemS(ll,stres) = sigmas;
                
                sigmas = ElemS2(1:nint,stres+11)'*shpS;
                ElemS(ll+nel1,stres) = sigmas;
                
            end
            
        end
%         
%         %Integration Loop
%         Vol = 0;
%         for ll = 1:lint
% 
%             %Evaluate first derivatives of basis functions at int. point
%             if nel == 4 || nel == 10
%               [Wgt,ss] =  int3d_t(ll,lint,ib);
%               [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
%               [Qxy,shgs,Jdet,be] = shgtt(xl(:,1:nel)+ul(:,1:nel),nel,shld,shls,nel,bf,der,be);
%             else
%               [Wgt,ss] =  intpntb(ll,lint,ib);
%               [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
%               [Qxy,shgs,Jdet,be] = shgb(xl(:,1:nel)+ul(:,1:nel),nel,shld,shls,nel,bf,der,be);
%             end
%             
%             [fi,JxX,F] = kine3d(Qxy,-ul(:,1:nel),nel,0); %this is equivalent to ikine2d
%             JxX = 1/JxX; %this is equivalent to ikine2d
%     %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
%             Jdet = Jdet/JxX;
% 
%             w = Wgt*Jdet*thick;
%             
%             Vol = Vol + w;
% 
%         end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
        
        
end