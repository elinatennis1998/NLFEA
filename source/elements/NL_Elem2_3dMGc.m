% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 10/01/2013
% cylindrical growth model
% first tested using new implementation of isotropic volumetric growth
%
% Equation numbers (#) refer to Computational Modelling of Isotropic Multiplicative Growth
% Himpel, Kuhl, Menzel, Steinmann, CMES, vol.8, no.2, pp.119-134, 2005
%
% Stress calculations in isw=25 are up-to-date

switch isw %Task Switch
%%
    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = 1*nen;
%%
    case 3
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(9,3*nel);
        Bmat2 = zeros(9,3*nel);

        emod = mateprop(1);   nue = mateprop(2);
%         mo  = mateprop(10:12); %fiber orientation vector
%         [~, mo] = VecNormalize(mo); % ensure mo is a unit vector
%         mo = mo';
        xlm = getlam([1 2 0 emod nue]);
        
        % Load Guass Integration Points

        if nel == 4
            lint = 16;11;5;16;
        elseif nel == 8
            lint = 8;
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
              [QXY, shgs, Jdet] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
              [QXY, shgs, Jdet] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            [F,JxX,fi] = kine3d(QXY,ul,nel,0);
                
            c1 = Wgt*Jdet;
            
            % Material update to find F=Fe*Fg, isotropic growth stretch
            % ratio theta=var+1
            i_var = hr(nha-1+ll);
            [var,fac1,fac2]=updt_vol2(F,i_var,mateprop,ndm);
            hr(nhb-1+ll) = var;
            
            delta  = eye(ndm); 
            theta  = var + 1;
            Fg_inv = (1/theta)*delta;% + (1-1/theta)*(mo*mo'); % (6)
            Fe = F * Fg_inv;
            Ce = Fe' * Fe; % (7)
            Fe_inv = inv(Fe);
            Ce_inv = inv(Ce);
            Je = det(Fe);

            % compute Se, Cijkl_e - using same material model as Kuhl
            [taue, cmat] = SigmaCmat3(Fe,Je,[1 2 0 emod nue],xlm); % taue is the Kirchhoff stress tensor (i.e. no Jacobian)
            P = tranr4(Fe_inv,Fe_inv); % transformation tensor F_iI*F_jJ
            P = [P' zeros(6,3); zeros(3,9)];
            Se = P*taue; % second PK stress
            Se_m = [Se(1) Se(4) Se(6); Se(4) Se(2) Se(5); Se(6) Se(5) Se(3)];
            Cmat_e = P*cmat*P';
            Pe = Fe*Se_m; % first PK stress

            sigma2 = taue;

            % Nonsymmetric terms:
            Ce_v = [Ce(1,1); Ce(2,2); Ce(3,3); 2*Ce(1,2); 2*Ce(2,3); 2*Ce(3,1); zeros(3,1)];
            CeCe = Cmat_e*Ce_v; % (55)
            intermed = CeCe + Se;
            intermed_m = [intermed(1) intermed(4) intermed(6); intermed(4) intermed(2) intermed(5); intermed(6) intermed(5) intermed(3)];
            
            AeFe = Fe*intermed_m; % this formula uses the acoustic to reference moduli transform:
            % A_iIjJ = S_IJ*g_ij + F_iK*C_KILJ*F_jL ->
            % A_iIjJ*F_jJ = F_jJ*S_JI*g_ij + F_iK*C_KILJ*F_jL*F_jJ ->
            % A_iIjJ*F_jJ = F_iJ*S_JI + F_iK*C_KILJ*C_LJ ->
            % A_iIjJ*F_jJ = F_iJ*(S_JI + C_JILK*C_LK)
            ten1 = fac1 *  AeFe;
            ten2 = fac2 * (Pe+AeFe);
            ten1v = [ten1(1,1); ten1(1,2); ten1(1,3); ten1(2,1); ten1(2,2); ten1(2,3); ten1(3,1); ten1(3,2); ten1(3,3)];
            ten2v = [ten2(1,1); ten2(1,2); ten2(1,3); ten2(2,1); ten2(2,2); ten2(2,3); ten2(3,1); ten2(3,2); ten2(3,3)];
            
            % equivalent method using only reference config quantities;
            % verified for random Fe+delta
            ten1V = fac1 *  intermed;
            ten2V = fac2 *  (intermed + Se);
            P = tranr4(Fe,Fe);
            P = [P' zeros(6,3); zeros(3,9)];
            ten1V = P*ten1V;
            ten2V = P*ten2V;
            
            Smat = ...
[    sigma2(1),        0,        0,           sigma2(4)/2,                 0,           sigma2(6)/2,           sigma2(4)/2,                 0,          -sigma2(6)/2
         0,    sigma2(2),        0,           sigma2(4)/2,           sigma2(5)/2,                 0,          -sigma2(4)/2,           sigma2(5)/2,                 0
         0,        0,    sigma2(3),                 0,           sigma2(5)/2,           sigma2(6)/2,                 0,          -sigma2(5)/2,           sigma2(6)/2
   sigma2(4)/2,  sigma2(4)/2,        0, sigma2(1)/4 + sigma2(2)/4,           sigma2(6)/4,           sigma2(5)/4, sigma2(2)/4 - sigma2(1)/4,           sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2,  sigma2(5)/2,           sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,           sigma2(4)/4,          -sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,           sigma2(4)/4
   sigma2(6)/2,        0,  sigma2(6)/2,           sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4,           sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4
   sigma2(4)/2, -sigma2(4)/2,        0, sigma2(2)/4 - sigma2(1)/4,          -sigma2(6)/4,           sigma2(5)/4, sigma2(1)/4 + sigma2(2)/4,          -sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2, -sigma2(5)/2,           sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,          -sigma2(4)/4,          -sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,          -sigma2(4)/4
  -sigma2(6)/2,        0,  sigma2(6)/2,          -sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4,          -sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4];
            
            Smat = Smat/theta + cmat/theta; %not sure why to divide by theta
            
            Qxy = QXY*Fe_inv; %map shape functions from intermediate to spatial
                
            for mm = 1:nel    
 Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                        0         Qxy(mm,2) 0         
                        0         0         Qxy(mm,3) 
                        Qxy(mm,2) Qxy(mm,1) 0         
                        0         Qxy(mm,3) Qxy(mm,2) 
                        Qxy(mm,3) 0         Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1) 0         
                        0         Qxy(mm,3) -Qxy(mm,2) 
                        -Qxy(mm,3) 0         Qxy(mm,1) ];
 Bmat2(:,3*mm-2:3*mm) = [QXY(mm,1) 0         0         
                         QXY(mm,2) 0         0         
                         QXY(mm,3) 0         0         
                         0         QXY(mm,1) 0         
                         0         QXY(mm,2) 0
                         0         QXY(mm,3) 0
                         0         0         QXY(mm,1) 
                         0         0         QXY(mm,2) 
                         0         0         QXY(mm,3) ];
            end
            
            ElemF = ElemF - c1*(Bmat'*(sigma2));
            
%             ElemK = ElemK + c1*(Bmat'*Smat*Bmat) + c1/theta*(Bmat2'*ten1v*ten2v'*Bmat2);
            ElemK = ElemK + c1*(Bmat'*Smat*Bmat) + c1/theta*(Bmat'*ten1V*ten2V'*Bmat);
            
        end %je
   ElemK;     
%%
    case -1

        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
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
              [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(je,lint,5);
              [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgb(xl,nel,shld,shls,nel,bf,der,be);
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
%%
    case 6
        
        ElemF = zeros(nst,1);
        Bmat = zeros(9,3*nel);

        lam = getlam(mateprop);
        
        % Load Guass Integration Points

        if nel == 4
            lint = 11;5;16;
        elseif nel == 8
            lint = 8;
        elseif nel == 10
            lint = 14;
%             lint = 27;
        else
            lint = 27;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 4 || nel == 10
                  [Wgt,ss] =  int3d_t(ll,lint,ib);
                  [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgtt(xl+ul,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl+ul,nel,shld,shls,nel,bf,der,be);
                end
                

            for mm = 1:nel    
 Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                        0         Qxy(mm,2) 0         
                        0         0         Qxy(mm,3) 
                        Qxy(mm,2) Qxy(mm,1) 0         
                        0         Qxy(mm,3) Qxy(mm,2) 
                        Qxy(mm,3) 0         Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1) 0         
                        0         Qxy(mm,3) -Qxy(mm,2) 
                        -Qxy(mm,3) 0         Qxy(mm,1) ];
            end
                
                
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
                
            c1 = Wgt*Jdet;
            
            sigma2 = SigmaCmat3(F,JxX,mateprop,lam);
            
            ElemF = ElemF - c1*(Bmat'*(sigma2));
            
        end %je
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
              [shl,shld,shls,be] = shltt(ss,nel,nen,der,bf);
              [Qxy, shgs, Jdet] = shgtt(xl+ul,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
              [Qxy, shgs, Jdet] = shgb(xl+ul,nel,shld,shls,nel,bf,der,be);
            end

            if nelP == 4 || nelP == 10
              [shlS,shld,shls] = shltt(ss,nelS,nel,0,0);
            else
              [shlS,shld,shls] = shlb(ss,nelS,nel,0,0);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
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
        ElemS2 = zeros(nel,npstr);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        spvec0 = I1;
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 4
            lint = 1;
            nint = 1;
        elseif nel == 8
%             lint = 4;
            lint = 8;
            nint = 1;
        elseif nel == 10
            lint = 11;
            nint = 4;
        else
            lint = 27;
            nint = 8;
        end
        
        der = 0;
        bf = 0;
        ib = 0;

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,nint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,nint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            end
            
            if nelP == 4 || nelP == 10
              shlp = shltt(ss,nelP,nel,0,0);
            else
              shlp = shlb(ss,nelP,nel,0,0);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;            
            
            i_var = hr(nha-1+ll);
            [var,fac1,fac2]=updt_vol2(F,i_var,mateprop,ndm);
            hr(nhb-1+ll) = var;
            
            delta  = eye(ndm); 
            theta  = var(1) + 1;
            Fg_inv = (1/theta)*delta;% + (1-1/theta)*(mo*mo');
            Fe = F * Fg_inv;
            Ce = Fe' * Fe;
            Fe_inv = inv(Fe);
            Ce_inv = inv(Ce);
            Je = det(Fe);

            % compute Se, Cijkl_e
            sigma = SigmaCmat3(Fe,Je,[1 2 0 emod nue],xlm); % taue is the Kirchhoff stress tensor (i.e. no Jacobian)
            sigma = sigma/Je;
            
            for stres = 1:npstr
            
            if stres <= 6 % stress components
                sigmas = sigma(stres);
            elseif stres >= 8
                if stres <= 10 % principal stresses
                sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma2);
                sigmas = psig(stres-7);
                else % hydrostatic stress
                sigmas = 1/3*sigma'*I1;
                end
            else % von Mises stress
                trs = sigma'*I1;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS2(ll,stres) = sigmas;
            
            end

        end %je
        
        % interpolate stress at nodes
        if nel == 4
            plist = [1 0 0 0
                     0 1 0 0
                     0 0 0 1];
        elseif nel == 8
            plist = [-1 1 1 -1 -1 1 1 -1
                     -1 -1 1 1 -1 -1 1 1
                     -1 -1 -1 -1 1 1 1 1];
        elseif nel == 10
            plist = [ 1.927050983124842 -0.309016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947 -0.309016994374947  0.809016994374947  0.809016994374947 -0.309016994374947 -0.309016994374947
                     -0.309016994374947  1.927050983124842 -0.309016994374947 -0.309016994374947  0.809016994374947  0.809016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947 -0.309016994374947
                     -0.309016994374947 -0.309016994374947 -0.309016994374947  1.927050983124842 -0.309016994374947 -0.309016994374947 -0.309016994374947  0.809016994374947  0.809016994374947  0.809016994374947];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 0 -sqr3 sqr3 0 0 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0
                     -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0 0 0 0 0];
        end
        
        for ll = 1:nelS
            
            r = plist(1,ll);
            s = plist(2,ll);
            t = plist(3,ll);
            shpS = sshp3d(r,s,t,nint);
            
            for stres = 1:npstr
                
                sigmas = ElemS2(1:nint,stres)'*shpS;
                ElemS(ll,stres) = sigmas;
                
            end
            
        end
        
        %Integration Loop
        Vol = 0;
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            w = Wgt*Jdet*thick;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
end