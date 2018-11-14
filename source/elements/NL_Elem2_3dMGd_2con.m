% Tim Truster
% 07/15/2013
% UIUC

% Density Growth Element
% Displacement and Density unknowns
% Based off of IJNME58-Kuhl paper on bone growth
% Gives same results for ex_patch3 as ex_patch2 with smooth loading, namely
% a constant density field. Convergence is quadratic.

switch isw %Task Switch
%%
    case 1
        
        if ndf > 4
            
            for i = 5:ndf
                lie(i,1) = 0;
            end
            
        end
%%
    case 3
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(10,4*nel);
        Nmat = zeros(4,4*nel);

        lam = getlam(mateprop);
        PatchE = mateprop(4);
        Patchv = mateprop(5);
        mu = PatchE/(2*(1+Patchv));
        rho0 = mateprop(6);   psi0 = mateprop(7);   
        expm = mateprop(8);   expn = mateprop(9);   dt   = mateprop(10);
        c_rho = mateprop(11); % Interactive force constant
        kIF = mateprop(12); % Interactive force constant
        
        % Load Guass Integration Points

        if nel == 4
            lint = 4;11;5;16;
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
                  [Qxy, shgs, Jdet] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
                end
                
            for mm = 1:nel    
 Bmat(:,4*mm-3:4*mm) = [Qxy(mm,1) 0         0         0
                        0         Qxy(mm,2) 0         0
                        0         0         Qxy(mm,3) 0
                        Qxy(mm,2) Qxy(mm,1) 0         0
                        0         Qxy(mm,3) Qxy(mm,2) 0
                        Qxy(mm,3) 0         Qxy(mm,1) 0
                        Qxy(mm,2) -Qxy(mm,1) 0         0
                        0         Qxy(mm,3) -Qxy(mm,2) 0
                        -Qxy(mm,3) 0         Qxy(mm,1) 0
                        0         0         0          shl(mm)];
 Nmat(:,4*mm-3:4*mm) = [shl(mm) 0         0         0
                        0         shl(mm) 0         0
                        0         0         shl(mm) 0
                        0         0         0          0];
            end
                
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
                
            c1 = Wgt*Jdet;
            
            % Density calcs
            C  = F'*F;
            I1 = trace(C);
            psi0_neo = lam/2 * log(JxX)^2  +  mu/2 *(I1 - ndm - 2*log(JxX));
            var = ul(4,:)*shl; % relative density is the unknown
            var_n = ul_n(4,:)*shl;
            rho = (1 + var) * rho0;
            rho_k0 = (1 + var_n) * rho0;
            facs= (rho/rho0)^expn;
            dres=  rho0*(1/dt    - c_rho*(expn-expm) * (rho/rho0)^(expn-expm) / rho *psi0_neo); % 1/dt - dR0/drho , also need to include drho/dvar
            Rres = (rho-rho_k0)/dt - c_rho*((rho/rho0)^(expn-expm)*psi0_neo-psi0); % residual of mass-balance
            
            % Displacement calcs
            [sigma2, cmat] = SigmaCmat3(F,JxX,mateprop,lam);
            
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
            
            Smat = Smat + cmat; % dP/dF
            svec = -c_rho*(rho/rho0)^(expn-expm)*sigma2; % dr0/dF
            svec2 = expn/rho*facs*sigma2*rho0; % dP/drho*drho/dvar
            Dmat = [facs*Smat svec2
                    svec' dres];
                
            % Interactive force
            
            InterForce3
            
            ElemF = ElemF - c1*(Bmat'*[facs*sigma2; Rres] + Nmat'*[IntFor; 0]);
            
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat + JxX*Nmat'*kIF*Nmat);
            
        end %je
   ElemK;   
%%
    case -1
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS

        ElemF = zeros(nst,1);
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

                ElemF(ndf*o-3) = ElemF(ndf*o-3) + F(1)*c1;

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(2)*c1;

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(3)*c1;

            end %o
            
        end %je
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
        
%%
    case 6
        
        ElemF = zeros(nst,1);
        Bmat = zeros(10,4*nel);
        Nmat = zeros(4,4*nel);

        lam = getlam(mateprop);
        PatchE = mateprop(4);
        Patchv = mateprop(5);
        mu = PatchE/(2*(1+Patchv));
        rho0 = mateprop(6);   psi0 = mateprop(7);   
        expm = mateprop(8);   expn = mateprop(9);   dt   = mateprop(10);
        c_rho = mateprop(11); % Interactive force constant
        kIF = mateprop(12); % Interactive force constant
        
        % Load Guass Integration Points

        if nel == 4
            lint = 4;11;5;16;
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
                  [Qxy, shgs, Jdet] = shgtt(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,ss] =  intpntb(ll,lint,ib);
                  [shl,shld,shls,be] = shlb(ss,nel,nen,der,bf);
                  [Qxy, shgs, Jdet] = shgb(xl+ul(1:3,:),nel,shld,shls,nel,bf,der,be);
                end
                
            for mm = 1:nel    
 Bmat(:,4*mm-3:4*mm) = [Qxy(mm,1) 0         0         0
                        0         Qxy(mm,2) 0         0
                        0         0         Qxy(mm,3) 0
                        Qxy(mm,2) Qxy(mm,1) 0         0
                        0         Qxy(mm,3) Qxy(mm,2) 0
                        Qxy(mm,3) 0         Qxy(mm,1) 0
                        Qxy(mm,2) -Qxy(mm,1) 0         0
                        0         Qxy(mm,3) -Qxy(mm,2) 0
                        -Qxy(mm,3) 0         Qxy(mm,1) 0
                        0         0         0          shl(mm)];
 Nmat(:,4*mm-3:4*mm) = [shl(mm) 0         0         0
                        0         shl(mm) 0         0
                        0         0         shl(mm) 0
                        0         0         0          0];
            end
                
            [fi,JxX,F] = kine3d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
                
            c1 = Wgt*Jdet;
            
            % Density calcs
            C  = F'*F;
            I1 = trace(C);
            psi0_neo = lam/2 * log(JxX)^2  +  mu/2 *(I1 - ndm - 2*log(JxX));
            var = ul(4,:)*shl; % relative density is the unknown
            var_n = ul_n(4,:)*shl;
            rho = (1 + var) * rho0;
            rho_k0 = (1 + var_n) * rho0;
            facs= (rho/rho0)^expn;
            dres=  rho0*(1/dt    - c_rho*(expn-expm) * (rho/rho0)^(expn-expm) / rho *psi0_neo); % 1/dt - dR0/drho , also need to include drho/dvar
            Rres = (rho-rho_k0)/dt - c_rho*((rho/rho0)^(expn-expm)*psi0_neo-psi0); % residual of mass-balance
            
            % Displacement calcs
            [sigma2, cmat] = SigmaCmat3(F,JxX,mateprop,lam);
                
            % Interactive force
            
            InterForce1
            
            ElemF = ElemF - c1*(Bmat'*[facs*sigma2; Rres] + Nmat'*[IntFor; 0]);
            
        end %je
   ElemF;   
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

        ElemS = zeros(nel,numstr+1);
        ElemS2 = zeros(nel,numstr);

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
            
            sigma = SigmaCmat3(F,JxX,mateprop,lam);
            sigma = sigma/JxX;
            
            for stres = 1:numstr
            
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
            
            for stres = 1:numstr
                
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
        ElemS(i,numstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
end