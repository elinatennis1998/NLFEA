% Tim Truster
% 05/06/2014
%
% Hypoelasto-plasticity element as a precursor for crystal plasticity
% Implementation uses Green-Naghdi stress rate, based off of the WARP3D
% implementation, version 17.5.0.
%
% Verified 05/08(9)/2014 against WARP3D. Using files WARP3D_Hypo3, I had to
% modify material properties slightly and watch the precision of values,
% but I would get good agreement. There is still the question of R*ddt*R'
% to answer as well.
%
% 5/9/14 morning: finished comparing with WARP3D and found all
% equivalences, using WARP3D_Hypo3. I found one bug in
% rstgp1(line=112)\gtmat1(line=117) that should be changed in order to
% realign with the theory. HOWEVER, changing it made the NR convergence
% slower for plasticity.
%
% 5/9/14 afternoon: talked with Mark Messner and fixed my engin. vs tensor.
% strain factors of 2. Now I agree with the WARP3D version.

switch isw %Task Switch
%%

    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        nh1 = (11+6+3)*nen;
        
%%
    case {3,6}
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(4,2*nel);
        Bmat2 = Bmat;

        PatchE = mateprop(4);
        Patchv = mateprop(5);
        beta = mateprop(6);
        tan_e = mateprop(7);
        yld_n1 = mateprop(8);
        ym_n1 = PatchE;
        nu_n1 = Patchv;
        ym_n = PatchE;
        nu_n = Patchv;
        hprime_n1 = tan_e*ym_n1/(ym_n1 - tan_e); % manual page 188 section 3.5.1
        
        % Load Guass Integration Points

        if nel == 3
            lint = 1;
        elseif nel == 4
            lint = 4;
        elseif nel == 6
            lint = 7;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 0;
        
        uldres = reshape(uld,nst,1);

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [QXY, shgs, Jdet,bub,xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [QXY, shgs, Jdet,bub,xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            % Step 1: Deformation gradient (1.143) and (1.144)
            [F,JxX,fi,Qxy] = kine2d(QXY,ul(:,1:nel),nel,1); %F_n+1
            Jdet = Jdet*JxX;
            ulm = 1/2*(ul_n + ul);
            [F2,JxX2,fi2,Qxy2] = kine2d(QXY,ulm(:,1:nel),nel,1); %F_n+1/2
            
            F = [F zeros(2,1); zeros(1,2) 1];
            F2 = [F2 zeros(2,1); zeros(1,2) 1];
            
            % Step 2: Compute polar decomposition (1.145) and (1.146)
            [Rmat, Umat] = poldec(F);
%             [Rmat, Umat_i] = rtcmp1(F);
            [R2, U2] = poldec(F2);
            
            % Form B matrix
            for mm = 1:nel  
% shape functions
 Nmat(:,2*mm-1:2*mm) = [shl(mm,1)     0          
                           0        shl(mm,1)    ];
% derivatives w.r.t. x_n+1^i
 Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0                  
                        0         Qxy(mm,2)  
                        Qxy(mm,2) Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1)];
% derivatives w.r.t. x_n+1/2^i
 Bmat2(:,2*mm-1:2*mm) = [Qxy2(mm,1) 0                 
                        0         Qxy2(mm,2) 
                        Qxy2(mm,2) Qxy2(mm,1)          
                        Qxy2(mm,2) -Qxy2(mm,1)];
            end
            % B-bar modification would follow right HERE
            
            
%% Branching to take care of WARP3D implementation: first solve
% is with linear K and no stress updates; then an update of the
% stresses and residual vector is performed
            if initia == 1 && step ~= 1% Follow lnstiff.f version
                
            
                % Load stresses and rotation from last converged step
                i = nh1-1+(ll-1)*20; %pointer to first history index
                cgn1 = hr(i+12:i+20);
                
                
%                 if step == 1 % set initial stiffness as in lnstff.f, line 646
                
%                 else % use tangent from last load step, and extrapolate stresses
                    % using recstr_cep_uddt_for_block idea: 
                    % stress @ n+1 = stress @ n + [Dt]* deps
                    
                    % In Warp3d, the order of operations is as follows:
                    % 1. extrapolate displacements (in NR_Loop3)
                    % 2. compute tangent stiffness with stifup.f using:
                    %    a. old displacements (Rmat3)
                    %    b. previous step TRUE Cauchy stress
                    %    c. for mm10, the cep from previous step
                    % 3. compute F_int using recstr_cep_uddt_for_block,
                    %    which is:
                    %    a. previous unrotated Cauchy stress
                    %    b. cep stored during stifup.f call above
                    %    c. old displacements (Rmat3)
                    
                    % 2. Reload cep and urcs from last converged step
                    % Can't use cnst10 because Jdet and Wgt get added on
                    % for now, FORCE it to be elastic; add this feature
                    % back later when I want plastic Warp3d elements
                    cep = cnst1(0,nu_n1,ym_n1,0,0,beta,0,Jdet,Wgt,0);
                    full_cep = cep/(Wgt*Jdet);
                    urcs_blk_n = cgn1(1:6);            % compute unrotated material tangent tensor
            % This is the ONLY thing that changes for each material model            
% c
% c              set material states and get [cep]. transform
% c              unrotated (material) -> spatial using qn1.
% c              get temperature dependent modulus and nu if needed.
% c
            
                    % Recompute TRUE Cauchy stress using rotations from last step
                    [F3,JxX3,fi3,Qxy3] = kine2d(QXY,ul_n(:,1:nel),nel,1); %F_n+1/2

                    F3 = [F3 zeros(2,1); zeros(1,2) 1];
            
                    % Step 2: Compute polar decomposition (1.145) and (1.146)
                    [Rmat3, Umat3] = poldec(F3);
                    qn3 = getrm1(Rmat3,2);
                    sigma3 = qn3*urcs_blk_n; % verified
            
% c
% c                       convert [Dt] from unrotated cauchy to cauchy
% c                       at current deformed configuration for geometric
% c                       nonlinear analysis. no computations
% c                       for cohesive or deformation plasticity. for UMAT with
% c                       hyperelastic formulations which use [F] to get strains, the
% c                       [Dt] stored in WARP3D is really for Cauchy stress - not
% c                       unrotated Cauchy stress. The code below skips the
% c                       rotation but may include the [Q] modification as
% c                       requested in user input.
% c
                    cep = ctran1(cep,qn3,sigma3,1,Jdet,Wgt); % verified
                    
                    cep = cep([1 2 4],[1 2 4]);

                    % As far as I can tell from the WARP manual and code (kg1), the
                    % geometric stress term looks exactly like the one I normally
                    % use, so that's why I copy it here.
                    Smat = Wgt*Jdet*[sigma3(1) 0  sigma3(3)/2 sigma3(3)/2
                            0 sigma3(2)  sigma3(3)/2 -sigma3(3)/2
                            sigma3(3)/2  sigma3(3)/2 (sigma3(2)+sigma3(1))/4 (sigma3(2)-sigma3(1))/4
                            sigma3(3)/2 -sigma3(3)/2 (sigma3(2)-sigma3(1))/4 (sigma3(2)+sigma3(1))/4];

                    Smat = Smat + [cep zeros(3,1); zeros(1,4)];
                    
                    
                    % 3. Extrapolate the stresses using cep and deps
                    del_eps = Bmat2*uldres; % (1.147)
                    del_eps = [del_eps(1) del_eps(2) 0 del_eps(3) 0 0]';
                    qn2 = getrm1(R2,1);
                    deps = qn2*del_eps(1:6);
                    
                    % execute according to recstr_cep_uddt_for_block
                    urcs_blk_n1 = urcs_blk_n + full_cep*deps;
                    
                    sigma = qn3*urcs_blk_n1;
                    sigma2 = [sigma; 0; 0; 0]; % verified
                    
                    sigma2 = sigma2([1 2 4 7]);
            
                
%                 end
                
            
            else
%% do full material update using mm01 using rknstr

                
            % Step 3: Spatial displacement gradient
            del_eps = Bmat2*uldres; % (1.147)
            del_eps = [del_eps(1) del_eps(2) 0 del_eps(3) 0 0]';
%             % Tim's version using conversion to tensors
%             Dten = [del_eps(1) del_eps(4)/2 del_eps(6)/2 
%                     del_eps(4)/2 del_eps(2) del_eps(5)/2
%                     del_eps(6)/2 del_eps(5)/2 del_eps(3)]; % (1.148)
%                 
%             % Step 4: Unrotated configuration
%             dten = R2'*Dten*R2; % (1.149)
%             
%             deps2 = [dten(1,1); dten(2,2); dten(3,3); 2*dten(1,2); 2*dten(2,3); 2*dten(3,1)];
            qn2 = getrm1(R2,1); % Mark Messner clarified: Tim forgot about engineering vs tensorial strain
            deps = qn2*del_eps(1:6);
            

%%
            % Step 5: Update the unrotated Cauchy stress
            % This is the ONLY thing that changes for each material model
            i = nh1-1+(ll-1)*20; %pointer to first history index
            history = hr(i+1:i+11);
            cgn = hr(i+12:i+20);
            [cgn1,history1,rtse,yield] = mm01(ym_n1,nu_n1,beta,hprime_n1,yld_n1,cgn,deps,history,ym_n,nu_n);
            i = nh2-1+(ll-1)*20; %pointer to first history index
            hr(i+1:i+11) = history1;
            hr(i+12:i+20) = cgn1;
            
            
%%            
%             % Step 6: Transform to spatial Cauchy stress
%             % Tim's version using conversion to tensors; agrees EXACTLY
%             % with WARP3D version
%             tten = [cgn1(1) cgn1(4) cgn1(6) 
%                     cgn1(4) cgn1(2) cgn1(5) 
%                     cgn1(6) cgn1(5) cgn1(3)];
%             
%             sigma3 = Rmat*tten*Rmat'; % (1.151)
%             sigma4 = [sigma3(1,1); sigma3(2,2); sigma3(3,3); sigma3(1,2); sigma3(2,3); sigma3(3,1); 0; 0; 0];
            qn1 = getrm1(Rmat,2);
            sigma = qn1*cgn1(1:6);
            sigma2 = [sigma; 0; 0; 0];
            
            
%%
            % compute unrotated material tangent tensor
            % This is the ONLY thing that changes for each material model            
% c
% c              set material states and get [cep]. transform
% c              unrotated (material) -> spatial using qn1.
% c              get temperature dependent modulus and nu if needed.
% c
            cep = cnst1(rtse,nu_n1,ym_n1,history1(2),history1(5),beta,history1(1),Jdet,Wgt,yield);
            
            
%%
% c
% c                       convert [Dt] from unrotated cauchy to cauchy
% c                       at current deformed configuration for geometric
% c                       nonlinear analysis. no computations
% c                       for cohesive or deformation plasticity. for UMAT with
% c                       hyperelastic formulations which use [F] to get strains, the
% c                       [Dt] stored in WARP3D is really for Cauchy stress - not
% c                       unrotated Cauchy stress. The code below skips the
% c                       rotation but may include the [Q] modification as
% c                       requested in user input.
% c
%             qn1 = getrm1(Rmat,2); % duplicate call to re-form qn1
            cep = ctran1(cep,qn1,sigma2,1,Jdet,Wgt);
            
            sigma2 = sigma2([1 2 4 7]);
            cep = cep([1 2 4],[1 2 4]);
                
            % As far as I can tell from the WARP manual and code (kg1), the
            % geometric stress term looks exactly like the one I normally
            % use, so that's why I copy it here.
            Smat = Wgt*Jdet*[sigma2(1) 0  sigma2(3)/2 sigma2(3)/2
                    0 sigma2(2)  sigma2(3)/2 -sigma2(3)/2
                    sigma2(3)/2  sigma2(3)/2 (sigma2(2)+sigma2(1))/4 (sigma2(2)-sigma2(1))/4
                    sigma2(3)/2 -sigma2(3)/2 (sigma2(2)-sigma2(1))/4 (sigma2(2)+sigma2(1))/4];

            Smat = Smat + [cep zeros(3,1); zeros(1,4)];
            
            end

            c1 = Wgt*Jdet;        

            % In rknifv.f, the shape function derivatives are computed
            % using ce_n1 which according to drive_eps_.../dupstr_blocked
            % is the current coordinates
            ElemF(1:ndf*nel) = ElemF(1:ndf*nel) - c1*Bmat'*(sigma2);
            
            ElemK(1:ndf*nel,1:ndf*nel) = ElemK(1:ndf*nel,1:ndf*nel) + (Bmat'*Smat*Bmat);
            
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

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            elseif nel == 4 || nel == 9
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

            if nel == 3 || nel == 6
              [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
              [QXY, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
              [QXY, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
                Traction = traction;
            
            c1 = Wgt*tm3*drdr*thick;
            for o=1:nel
                don = shl(o);
                F = don*Traction';

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(1)*c1;

                ElemF(ndf*o-0)   = ElemF(ndf*o-0)   + F(2)*c1;

            end %o

        end %ie
        
        ElemF;
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
%%        
    case 25 %Stress Projection2

        PatchE = mateprop(4);
        Patchv = mateprop(5);
        beta = mateprop(6);
        tan_e = mateprop(7);
        yld_n1 = mateprop(8);
        ym_n1 = PatchE;
        nu_n1 = Patchv;
        ym_n = PatchE;
        nu_n = Patchv;
        hprime_n1 = tan_e*ym_n1/(ym_n1 - tan_e); % manual page 188 section 3.5.1

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0];
        spvec0 = I1;
        
        thick = 1;
        
        %%%%% NOTE: These will not be correct for history dependent
        %%%%% materials since we are pulling from the wrong Gauss point
        
        % Load Guass Integration Points

        if nel == 3
            lint = 1;
            nint = 1;
        elseif nel == 4
%             lint = 4;
            lint = 4;
            nint = 4;
        elseif nel == 6
            lint = 7;
            nint = 4;
        else
            lint = 9;
            nint = 8;
        end
        
        der = 0;
        bf = 0;
        ib = 0;
        
        uldres = reshape(uld,nst,1);

        %Stress Loop
        for ll = 1:nint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [QXY, shgs, Jdet,bub,xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [QXY, shgs, Jdet,bub,xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            % Step 1: Deformation gradient (1.143) and (1.144)
            [F,JxX,fi,Qxy] = kine2d(QXY,ul(:,1:nel),nel,1); %F_n+1
            Jdet = Jdet*JxX;
            ulm = 1/2*(ul_n + ul);
            [F2,JxX2,fi2,Qxy2] = kine2d(QXY,ulm(:,1:nel),nel,1); %F_n+1/2
            
            F = [F zeros(2,1); zeros(1,2) 1];
            F2 = [F2 zeros(2,1); zeros(1,2) 1];
            
            % Step 2: Compute polar decomposition (1.145) and (1.146)
            [Rmat, Umat] = poldec(F);
%             [Rmat, Umat_i] = rtcmp1(F);
            [R2, U2] = poldec(F2);
            
            % Form B matrix
            for mm = 1:nel  
% shape functions
 Nmat(:,2*mm-1:2*mm) = [shl(mm,1)     0          
                           0        shl(mm,1)    ];
% derivatives w.r.t. x_n+1^i
 Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0                  
                        0         Qxy(mm,2)  
                        Qxy(mm,2) Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1)];
% derivatives w.r.t. x_n+1/2^i
 Bmat2(:,2*mm-1:2*mm) = [Qxy2(mm,1) 0                 
                        0         Qxy2(mm,2) 
                        Qxy2(mm,2) Qxy2(mm,1)          
                        Qxy2(mm,2) -Qxy2(mm,1)];
            end
            % B-bar modification would follow right HERE
            
            % Step 3: Spatial displacement gradient
            del_eps = Bmat2*uldres; % (1.147)
            del_eps = [del_eps(1) del_eps(2) 0 del_eps(3) 0 0]';
%             % Tim's version using conversion to tensors
%             Dten = [del_eps(1) del_eps(4)/2 del_eps(6)/2 
%                     del_eps(4)/2 del_eps(2) del_eps(5)/2
%                     del_eps(6)/2 del_eps(5)/2 del_eps(3)]; % (1.148)
%                 
%             % Step 4: Unrotated configuration
%             dten = R2'*Dten*R2; % (1.149)
%             
%             deps2 = [dten(1,1); dten(2,2); dten(3,3); 2*dten(1,2); 2*dten(2,3); 2*dten(3,1)];
            qn2 = getrm1(R2,1); % Mark Messner clarified: Tim forgot about engineering vs tensorial strain
            deps = qn2*del_eps(1:6);
            

%%
            % Step 5: Update the unrotated Cauchy stress
            % This is the ONLY thing that changes for each material model
            i = nh1-1+(ll-1)*20; %pointer to first history index
            history = hr(i+1:i+11);
            cgn = hr(i+12:i+20);
            [cgn1,history1,rtse,yield] = mm01(ym_n1,nu_n1,beta,hprime_n1,yld_n1,cgn,deps,history,ym_n,nu_n);
            i = nh2-1+(ll-1)*20; %pointer to first history index
            hr(i+1:i+11) = history1;
            hr(i+12:i+20) = cgn1;
            
            
%%            
%             % Step 6: Transform to spatial Cauchy stress
%             % Tim's version using conversion to tensors; agrees EXACTLY
%             % with WARP3D version
%             tten = [cgn1(1) cgn1(4) cgn1(6) 
%                     cgn1(4) cgn1(2) cgn1(5) 
%                     cgn1(6) cgn1(5) cgn1(3)];
%             
%             sigma3 = Rmat*tten*Rmat'; % (1.151)
%             sigma4 = [sigma3(1,1); sigma3(2,2); sigma3(3,3); sigma3(1,2); sigma3(2,3); sigma3(3,1); 0; 0; 0];
            qn1 = getrm1(Rmat,2);
            sigma = qn1*cgn1(1:6);
            
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
                sigmas = sqrt(3/2*(dsig'*diag([1 1 1 2 2 2])*dsig));
            end
            
            ElemS2(ll,stres) = sigmas;
            
            end

        end %je
        
        % interpolate stress at nodes
        if nel == 3
            plist = [1 0 0 0
                     0 1 0 0
                     0 0 0 1];
        elseif nel == 4
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3
                     -sqr3 -sqr3 sqr3 sqr3];
%             plist = [-1 1 1 -1 -1 1 1 -1
%                      -1 -1 1 1 -1 -1 1 1
%                      -1 -1 -1 -1 1 1 1 1];
        elseif nel == 6
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
            shpS = sshp2d(r,s,nint);
            
            for stres = 1:npstr
                
                sigmas = ElemS2(1:nint,stres)'*shpS;
                ElemS(ll,stres) = sigmas;
                
            end
            
        end
        
        %Integration Loop
        Vol = 0;
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            [F,JxX,fi] = kine2d(QXY,ul,nel,0); %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            w = Wgt*Jdet*thick;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemS(i,npstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
%%        
    case 26 % Element Stress

        ElemS = zeros(1,nestr);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0];
        spvec0 = I1;
        
        %%%%% NOTE: These will not be correct for history dependent
        %%%%% materials since we are pulling from the wrong Gauss point
        
        der = 0;
        bf = 0;
        ib = 0;


            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [shl,shld,shls,be] = shlt(1/3,1/3,nel,nel,der,bf);
              [QXY, shgs, Jdet,bub,xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [shl,shld,shls,be] = shlq(0,0,nel,nel,der,bf);
              [QXY, shgs, Jdet,bub,xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            % Step 1: Deformation gradient (1.143) and (1.144)
            [F,JxX,fi,Qxy] = kine2d(QXY,ul(:,1:nel),nel,1); %F_n+1
            Jdet = Jdet*JxX;
            ulm = 1/2*(ul_n + ul);
            [F2,JxX2,fi2,Qxy2] = kine2d(QXY,ulm(:,1:nel),nel,1); %F_n+1/2
            
            F = [F zeros(2,1); zeros(1,2) 1];
            F2 = [F2 zeros(2,1); zeros(1,2) 1];
            
            % Step 2: Compute polar decomposition (1.145) and (1.146)
            [Rmat, Umat] = poldec(F);
%             [Rmat, Umat_i] = rtcmp1(F);
            [R2, U2] = poldec(F2);
            
            % Form B matrix
            for mm = 1:nel  
% shape functions
 Nmat(:,2*mm-1:2*mm) = [shl(mm,1)     0          
                           0        shl(mm,1)    ];
% derivatives w.r.t. x_n+1^i
 Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0                  
                        0         Qxy(mm,2)  
                        Qxy(mm,2) Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1)];
% derivatives w.r.t. x_n+1/2^i
 Bmat2(:,2*mm-1:2*mm) = [Qxy2(mm,1) 0                 
                        0         Qxy2(mm,2) 
                        Qxy2(mm,2) Qxy2(mm,1)          
                        Qxy2(mm,2) -Qxy2(mm,1)];
            end
            % B-bar modification would follow right HERE
            
            % Step 3: Spatial displacement gradient
            del_eps = Bmat2*uldres; % (1.147)
            del_eps = [del_eps(1) del_eps(2) 0 del_eps(3) 0 0]';
%             % Tim's version using conversion to tensors
%             Dten = [del_eps(1) del_eps(4)/2 del_eps(6)/2 
%                     del_eps(4)/2 del_eps(2) del_eps(5)/2
%                     del_eps(6)/2 del_eps(5)/2 del_eps(3)]; % (1.148)
%                 
%             % Step 4: Unrotated configuration
%             dten = R2'*Dten*R2; % (1.149)
%             
%             deps2 = [dten(1,1); dten(2,2); dten(3,3); 2*dten(1,2); 2*dten(2,3); 2*dten(3,1)];
            qn2 = getrm1(R2,1); % Mark Messner clarified: Tim forgot about engineering vs tensorial strain
            deps = qn2*del_eps(1:6);
            

%%
            % Step 5: Update the unrotated Cauchy stress
            % This is the ONLY thing that changes for each material model
            i = nh1-1+(ll-1)*20; %pointer to first history index
            history = hr(i+1:i+11);
            cgn = hr(i+12:i+20);
            [cgn1,history1,rtse,yield] = mm01(ym_n1,nu_n1,beta,hprime_n1,yld_n1,cgn,deps,history,ym_n,nu_n);
            i = nh2-1+(ll-1)*20; %pointer to first history index
            hr(i+1:i+11) = history1;
            hr(i+12:i+20) = cgn1;
            
            
%%            
%             % Step 6: Transform to spatial Cauchy stress
%             % Tim's version using conversion to tensors; agrees EXACTLY
%             % with WARP3D version
%             tten = [cgn1(1) cgn1(4) cgn1(6) 
%                     cgn1(4) cgn1(2) cgn1(5) 
%                     cgn1(6) cgn1(5) cgn1(3)];
%             
%             sigma3 = Rmat*tten*Rmat'; % (1.151)
%             sigma4 = [sigma3(1,1); sigma3(2,2); sigma3(3,3); sigma3(1,2); sigma3(2,3); sigma3(3,1); 0; 0; 0];
            qn1 = getrm1(Rmat,2);
            sigma = qn1*cgn1(1:6);
            sigma2 = [sigma; 0; 0; 0];
            
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
                sigmas = sqrt(3/2*(dsig'*diag([1 1 1 2 2 2])*dsig));
            end
            
            ElemS(1,stres) = sigmas;
            
            end

    case 24
        
        ElemP = zeros(12,nel);
        
        Bmat = zeros(9,3*nel);
        Bmat2 = Bmat;

        PatchE = mateprop(4);
        Patchv = mateprop(5);
        beta = mateprop(6);
        tan_e = mateprop(7);
        yld_n1 = mateprop(8);
        ym_n1 = PatchE;
        nu_n1 = Patchv;
        ym_n = PatchE;
        nu_n = Patchv;
        hprime_n1 = tan_e*ym_n1/(ym_n1 - tan_e);
        
        % Load Guass Integration Points

        if nel == 3
            lint = 1;
        elseif nel == 4
            lint = 4;
        elseif nel == 6
            lint = 7;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 0;
        
        uldres = reshape(uld,nst,1);

        %Integration Loop
        for ll = 1:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              [Wgt,litr,lits] =  intpntt(ll,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [QXY, shgs, Jdet,bub,xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,litr,lits] =  intpntq(ll,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [QXY, shgs, Jdet,bub,xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            % Step 1: Deformation gradient (1.143) and (1.144)
            [F,JxX,fi,Qxy] = kine2d(QXY,ul(:,1:nel),nel,1); %F_n+1
            Jdet = Jdet*JxX;
            ulm = 1/2*(ul_n + ul);
            [F2,JxX2,fi2,Qxy2] = kine2d(QXY,ulm(:,1:nel),nel,1); %F_n+1/2
            
            F = [F zeros(2,1); zeros(1,2) 1];
            F2 = [F2 zeros(2,1); zeros(1,2) 1];
            
            % Step 2: Compute polar decomposition (1.145) and (1.146)
            [Rmat, Umat] = poldec(F);
%             [Rmat, Umat_i] = rtcmp1(F);
            [R2, U2] = poldec(F2);
            
            % Form B matrix
            for mm = 1:nel  
% shape functions
 Nmat(:,2*mm-1:2*mm) = [shl(mm,1)     0          
                           0        shl(mm,1)    ];
% derivatives w.r.t. x_n+1^i
 Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0                  
                        0         Qxy(mm,2)  
                        Qxy(mm,2) Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1)];
% derivatives w.r.t. x_n+1/2^i
 Bmat2(:,2*mm-1:2*mm) = [Qxy2(mm,1) 0                 
                        0         Qxy2(mm,2) 
                        Qxy2(mm,2) Qxy2(mm,1)          
                        Qxy2(mm,2) -Qxy2(mm,1)];
            end
            % B-bar modification would follow right HERE
            
            % Step 3: Spatial displacement gradient
            del_eps = Bmat2*uldres; % (1.147)
%             % Tim's version using conversion to tensors
%             Dten = [del_eps(1) del_eps(4)/2 del_eps(6)/2 
%                     del_eps(4)/2 del_eps(2) del_eps(5)/2
%                     del_eps(6)/2 del_eps(5)/2 del_eps(3)]; % (1.148)
%                 
%             % Step 4: Unrotated configuration
%             dten = R2'*Dten*R2; % (1.149)
%             
%             deps2 = [dten(1,1); dten(2,2); dten(3,3); 2*dten(1,2); 2*dten(2,3); 2*dten(3,1)];
            qn2 = getrm1(R2,1); % Mark Messner clarified: Tim forgot about engineering vs tensorial strain
            deps = qn2*del_eps(1:6);
            

%%
            % Step 5: Update the unrotated Cauchy stress
            % This is the ONLY thing that changes for each material model
            i = nh1-1+(ll-1)*20; %pointer to first history index
            history = hr(i+1:i+11);
            cgn = hr(i+12:i+20);
            [cgn1,history1,rtse,yield] = mm01(ym_n1,nu_n1,beta,hprime_n1,yld_n1,cgn,deps,history,ym_n,nu_n);
            i = nh2-1+(ll-1)*20; %pointer to first history index
            hr(i+1:i+11) = history1;
            hr(i+12:i+20) = cgn1;
            
            
%%            
%             % Step 6: Transform to spatial Cauchy stress
%             % Tim's version using conversion to tensors; agrees EXACTLY
%             % with WARP3D version
%             tten = [cgn1(1) cgn1(4) cgn1(6) 
%                     cgn1(4) cgn1(2) cgn1(5) 
%                     cgn1(6) cgn1(5) cgn1(3)];
%             
%             sigma3 = Rmat*tten*Rmat'; % (1.151)
%             sigma4 = [sigma3(1,1); sigma3(2,2); sigma3(3,3); sigma3(1,2); sigma3(2,3); sigma3(3,1); 0; 0; 0];
            qn1 = getrm1(Rmat,2);
            sigma = qn1*cgn1(1:6);
            sigma2 = [sigma; 0; 0; 0];
            
            
%%
            % compute unrotated material tangent tensor
            % This is the ONLY thing that changes for each material model            
% c
% c              set material states and get [cep]. transform
% c              unrotated (material) -> spatial using qn1.
% c              get temperature dependent modulus and nu if needed.
% c
            cep = cnst1(rtse,nu_n1,ym_n1,history1(2),history1(5),beta,history1(1),Jdet,Wgt,yield);
            
            
%%
% c
% c                       convert [Dt] from unrotated cauchy to cauchy
% c                       at current deformed configuration for geometric
% c                       nonlinear analysis. no computations
% c                       for cohesive or deformation plasticity. for UMAT with
% c                       hyperelastic formulations which use [F] to get strains, the
% c                       [Dt] stored in WARP3D is really for Cauchy stress - not
% c                       unrotated Cauchy stress. The code below skips the
% c                       rotation but may include the [Q] modification as
% c                       requested in user input.
% c
%             qn1 = getrm1(Rmat,2); % duplicate call to re-form qn1
            cep = ctran1(cep,qn1,sigma2,1,Jdet,Wgt);
                
                % output stuff
                
                ElemP(1,ll) = sigma2(1);
                ElemP(2,ll) = sigma2(2);
                ElemP(3,ll) = sigma2(3);
                ElemP(4,ll) = sigma2(4);
                ElemP(5,ll) = history1(5);
                ElemP(6,ll) = history1(6);
                ElemP(7,ll) = history1(7);
                ElemP(8,ll) = history1(8);
                ElemP(9,ll) = history1(9);
                ElemP(10,ll) = cep(1,1)/(Jdet*Wgt);
                ElemP(11,ll) = cep(2,2)/(Jdet*Wgt);
                ElemP(12,ll) = cep(4,4)/(Jdet*Wgt);
            
        end %je
            
    case 40 % Initialize history terms
        
        if nel == 3
            lint = 1;
        elseif nel == 4
            lint = 4;
        elseif nel == 6
            lint = 7;
        else
            lint = 9;
        end
        
        zero = 0;
        PatchE = mateprop(4);
        sigma_o = mateprop(8);
        tan_e = mateprop(7);
        ym_n1 = PatchE;
        hprime = tan_e*ym_n1/(ym_n1 - tan_e);
%     sigy = mateprop(8);
%     K = mateprop(9);
%     H = mateprop(10);
%     eta = mateprop(11);

        root3 = sqrt(3);
            kn           = sigma_o / root3;
        
        % Loop over integration points
        for l = 1:lint
                
            % Store history variables
            i = nh1-1+(l-1)*20; %pointer to first history index
            hr(i+1) = zero;
            hr(i+2) = kn;
            hr(i+3) = zero;
            hr(i+4) = 3;
            hr(i+5) = hprime;
            hr(i+6:i+20) = zero;
%             history(i+7) = zero
%             history(i+8) = zero
%             history(i+9) = zero
%             history(i+10) = zero
%             history(i+11) = zero

        end %je
        
end