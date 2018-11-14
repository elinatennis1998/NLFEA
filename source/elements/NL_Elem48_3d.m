% Tim Truster
% 06/01/2014
%
% Crystal plasticity element based on Mark's first implementation
%
% 06/24/2014 - revisions to set/allocate material properties
% 07/31/2014 - Subcycling added to this element and to NR_Loop3
%              Verified to give identical results for MarkCPTest6Voche, Lie
%              and Mark versions (ReMeth), numsubcycle=0 and numReiter=0;
% 04/29/2016 - New element, only for bricks, has 8 integration points but
%              only uses the stress at the first point; gradient model is
%              disabled
[max_slip_sys,max_uhard] = maxparamsCP;
vec1 = zeros(max_uhard,1);
vec2 = zeros(max_uhard,1);
if ~(exist('tstep','var')) % set dt for CP material routine to tstep from input file
    dtWARP = 1.0e+10;
else
    dtWARP = tstep;
end

switch isw %Task Switch
%%

    case 1
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end
        
        hist_sz = 76+12+max_slip_sys+(30+max_slip_sys+3*max_uhard);
        nh1CP = hist_sz;
        nh1 = hist_sz*8;
        % Sets number of Lie variables for interpolation, later
        numLie = 18; % tensors Re and Ue
        istv = 11; % number of stresses per node for post-processing
        iste = 11*nen+38+6+5+2+3+9+9; % number of stresses per element (total number, including all integration points)
        
%%
    case {3,6}
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(9,3*nel);
        Bmat2 = Bmat;
        
%         % Use this to test whether adaptive time stepping works
%         if step == 1 && iter == 1 && (numincrem1 == 1 || numincrem2 == 1)
%             failFEelem = 1;
%             return
%         end
        
        % Load Guass Integration Points

        lint = 8;
        ngp = lint;
        der = 0;
        bf = 0;
        ib = 0;
        
        uldres = reshape(uld,nst,1);
          
        
%% Set up local work and history data in format expected by mm10
%                     case ( 10 )
        local_work = setup_mm10_rknstr( 1, mateprop.cp_prop, mateprop.cp_other, elem, dtWARP); % dt = 1.0e+10;
%      from line 410,rknstr in WARP3D;
        ncrystals = 1;
        gp_temp_inc = zeros(lint,1);
        % technically, this should only be for one Gauss point but for all
        % elements in the block; but, only the first number ever gets
        % accessed.
        gp_temps = mateprop.cp_allTemps*ones(lint,1); % set a constant temperature for all Gauss points, from input file.
        stepW = step;
        if stepW == 0
            stepW = 1;
        end
        local_work.step = stepW;

        % Copy history data for all integration points into local_work
        for ll = 1:1
            
            i = nh1-1+(ll-1)*nh1CP; %pointer to first history index
            history_n = hr(i+1:i+hist_sz);
            
            local_work.elem_hist(1,1:hist_sz,ll) = history_n;
            
            
        end
        

%% Integration Loop
        for ll = 1:1

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet,bbb,xs] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet,bbb,xs] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            % Step 1: Deformation gradient (1.143) and (1.144)
            [F,JxX,fi,Qxy] = kine3d(QXY,ul(:,1:nel),nel,1); %F_n+1
            Jdet = Jdet*JxX;
            ulm = 1/2*(ul_n + ul);
            [F2,JxX2,fi2,Qxy2] = kine3d(QXY,ulm(:,1:nel),nel,1); %F_n+1/2
            
            % Step 2: Compute polar decomposition (1.145) and (1.146)
            [Rmat, Umat] = poldec(F);
%             [Rmat, Umat_i] = rtcmp1(F);
            [R2, U2] = poldec(F2);
            
            % Form B matrix
            for mm = 1:nel  
% shape functions
 Nmat(:,3*mm-2:3*mm) = [shl(mm,1)     0          0
                           0        shl(mm,1)     0
                           0            0       shl(mm,1) ];
% derivatives w.r.t. x_n+1^i
 Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                        0         Qxy(mm,2) 0         
                        0         0         Qxy(mm,3) 
                        Qxy(mm,2) Qxy(mm,1) 0         
                        0         Qxy(mm,3) Qxy(mm,2) 
                        Qxy(mm,3) 0         Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1) 0         
                        0         Qxy(mm,3) -Qxy(mm,2) 
                        -Qxy(mm,3) 0         Qxy(mm,1) ];
% derivatives w.r.t. x_n+1/2^i
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
            % B-bar modification would follow right HERE
            
            
%% Branching to take care of WARP3D implementation: first solve
% is with linear K and no stress updates; then an update of the
% stresses and residual vector is performed
            if initia == 1 % Follow lnstiff.f version
                
            
                % Load stresses and rotation from last converged step
                i = nh1-1+(ll-1)*nh1CP; %pointer to first history index
                history_n = hr(i+1:i+hist_sz);
                
                co = 76+12+max_slip_sys;
%                 cn = 76+max_slip_sys+(1)*(25+max_uhard)-1;
                
                cgn1 = history_n(co:co+5);
                
                
                if step == 1 % set initial stiffness as in lnstff.f, line 646
                    
                qn1 = getrm1(Rmat,2);
                sigma = qn1*cgn1(1:6);
                sigma2 = [sigma; 0; 0; 0]; % get stress from last time step
                
                local_work = local_work_lin_ek;
                local_work.cp_stiff(1,:,:) = mateprop.cp_prop.elast_stiff;
                
%                 if (imatprp(107,matnum) .eq. 1) then
                angles(1) = mateprop.cp_other.angles(1);
                angles(2) = mateprop.cp_other.angles(2);
                angles(3) = mateprop.cp_other.angles(3);
%                   elseif (imatprp(107,matnum) .eq. 2) then
%                         osn = data_offset(elnum)
%                         angles(1:3) = angle_input(osn,ci,1:3)
%                   else
%                         write (out,9502)
%                         call die_gracefully
%                   end if
%                   aci = imatprp(102,matnum);
%                   ati = imatprp(103,matnum);
% c                 Call a helper to get the crystal -> reference rotation
%                   if (ati == 1) %then
                atype = 'degrees';
%                   elseif (ati .eq. 2) then
%                         atype = "radians"
%                   else
%                         write(out,9503)
%                         call die_gracefully
%                   end if

%                   if (aci == 1) %then
                aconv='kocks';
%                   else
%                         write(out,9504)
%                         call die_gracefully
%                   end if
                local_work.cp_g_rot(1,1:3,1:3,1) = ...
                    mm10_rotation_matrix(angles, aconv, atype);
                
                local_work.det_jac_block = Jdet;
                local_work.weights = Wgt;
                local_work = lnstff10(1, 1, local_work);
                cep = squeeze(local_work.cep);
                
                Smat = [cep zeros(6,3); zeros(3,9)];
                
                else % use tangent from last load step, and extrapolate stresses
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
                    cep = reshape(squeeze(local_work.elem_hist(1,1:36,ll)),6,6);
                    full_cep = cep;
                    cep = Wgt*Jdet*cep;
                    urcs_blk_n = cgn1(1:6);
            
                    % Recompute TRUE Cauchy stress using rotations from last step
                    [F3,JxX3,fi3,Qxy3] = kine3d(QXY,ul_n(:,1:nel),nel,1); %F_n+1/2
            
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
                    
            Smat = Wgt*Jdet*...
[    sigma3(1),        0,        0,           sigma3(4)/2,                 0,           sigma3(6)/2,           sigma3(4)/2,                 0,          -sigma3(6)/2
         0,    sigma3(2),        0,           sigma3(4)/2,           sigma3(5)/2,                 0,          -sigma3(4)/2,           sigma3(5)/2,                 0
         0,        0,    sigma3(3),                 0,           sigma3(5)/2,           sigma3(6)/2,                 0,          -sigma3(5)/2,           sigma3(6)/2
   sigma3(4)/2,  sigma3(4)/2,        0, sigma3(1)/4 + sigma3(2)/4,           sigma3(6)/4,           sigma3(5)/4, sigma3(2)/4 - sigma3(1)/4,           sigma3(6)/4,          -sigma3(5)/4
         0,  sigma3(5)/2,  sigma3(5)/2,           sigma3(6)/4, sigma3(2)/4 + sigma3(3)/4,           sigma3(4)/4,          -sigma3(6)/4, sigma3(3)/4 - sigma3(2)/4,           sigma3(4)/4
   sigma3(6)/2,        0,  sigma3(6)/2,           sigma3(5)/4,           sigma3(4)/4, sigma3(1)/4 + sigma3(3)/4,           sigma3(5)/4,          -sigma3(4)/4, sigma3(1)/4 - sigma3(3)/4
   sigma3(4)/2, -sigma3(4)/2,        0, sigma3(2)/4 - sigma3(1)/4,          -sigma3(6)/4,           sigma3(5)/4, sigma3(1)/4 + sigma3(2)/4,          -sigma3(6)/4,          -sigma3(5)/4
         0,  sigma3(5)/2, -sigma3(5)/2,           sigma3(6)/4, sigma3(3)/4 - sigma3(2)/4,          -sigma3(4)/4,          -sigma3(6)/4, sigma3(2)/4 + sigma3(3)/4,          -sigma3(4)/4
  -sigma3(6)/2,        0,  sigma3(6)/2,          -sigma3(5)/4,           sigma3(4)/4, sigma3(1)/4 - sigma3(3)/4,          -sigma3(5)/4,          -sigma3(4)/4, sigma3(1)/4 + sigma3(3)/4];

                    Smat = Smat + [cep zeros(6,3); zeros(3,9)];
                    
                    
                    % 3. Extrapolate the stresses using cep and deps
                    del_eps = Bmat2*uldres; % (1.147)
                    qn2 = getrm1(R2,1);
                    deps = qn2*del_eps(1:6);
                    
                    % execute according to recstr_cep_uddt_for_block
                    urcs_blk_n1 = urcs_blk_n + full_cep*deps;
                    
                    sigma = qn3*urcs_blk_n1;
                    sigma2 = [sigma; 0; 0; 0]; % verified
            
                
                end
                
            
            else
%% do full material update using mm10 using rknstr
                
        
%         % Use this to test whether adaptive time stepping works
%         if step == 1 && iter == 0 && (numincrem1 == 1 || numincrem2 == 1)
%             failFEelem = 1;
%             return
%         end
            
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
            
%             i = nh1-1+(ll-1)*nh1CP; %pointer to first history index
%             history_n = hr(i+1:i+132);
%             ncrystals = 1;
%             gp_temp_inc = 0;
%             gp_temps = 0;
%             hist_sz = 2000;
%             stepW = step;
%             if stepW == 0
%                 stepW = 1;
%             end
%             local_work.step = stepW;
%             history_n = reshape(history_n,1,132);
            local_work.jac(1,:,:) = xs;
            uddt = reshape(deps,1,6);
            local_work.rot_blk_n1(1, 1:9, ll) = reshape(Rmat,1,9);
            
             [local_work.elem_hist1(1,1:hist_sz,ll),local_work.elem_hist(1,1:hist_sz,ll), ...
            local_work, gp_temps, gp_temp_inc] = mm10(ll, 1,...
            ncrystals, hist_sz, local_work.elem_hist(1,1:hist_sz,ll), local_work, uddt,...
            gp_temps, gp_temp_inc, 1,subcycle);
        
            if local_work.material_cut_step
                failFEelem = 1;
                return
            end
        
            i = nh2-1+(ll-1)*nh1CP; %pointer to first history index
            hr(i+1:i+hist_sz) = local_work.elem_hist1(1,1:hist_sz,ll);
            cgn1 = squeeze(local_work.urcs_blk_n1(1,1:6,ll))';
            
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
                local_work2 = local_work_tan_ek;
                local_work2.det_jac_block(ll) = Jdet;
                local_work2.weights(ll) = Wgt;
            local_work2 = cnst10(ll, 1, local_work.elem_hist1(1,1:hist_sz,ll), local_work2);
            cep = squeeze(local_work2.cep(1,:,:));
            
            
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
            
                
            % As far as I can tell from the WARP manual and code (kg1), the
            % geometric stress term looks exactly like the one I normally
            % use, so that's why I copy it here.
            Smat = Wgt*Jdet*...
[    sigma2(1),        0,        0,           sigma2(4)/2,                 0,           sigma2(6)/2,           sigma2(4)/2,                 0,          -sigma2(6)/2
         0,    sigma2(2),        0,           sigma2(4)/2,           sigma2(5)/2,                 0,          -sigma2(4)/2,           sigma2(5)/2,                 0
         0,        0,    sigma2(3),                 0,           sigma2(5)/2,           sigma2(6)/2,                 0,          -sigma2(5)/2,           sigma2(6)/2
   sigma2(4)/2,  sigma2(4)/2,        0, sigma2(1)/4 + sigma2(2)/4,           sigma2(6)/4,           sigma2(5)/4, sigma2(2)/4 - sigma2(1)/4,           sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2,  sigma2(5)/2,           sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,           sigma2(4)/4,          -sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,           sigma2(4)/4
   sigma2(6)/2,        0,  sigma2(6)/2,           sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4,           sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4
   sigma2(4)/2, -sigma2(4)/2,        0, sigma2(2)/4 - sigma2(1)/4,          -sigma2(6)/4,           sigma2(5)/4, sigma2(1)/4 + sigma2(2)/4,          -sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2, -sigma2(5)/2,           sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,          -sigma2(4)/4,          -sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,          -sigma2(4)/4
  -sigma2(6)/2,        0,  sigma2(6)/2,          -sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4,          -sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4];

            Smat = Smat + [cep zeros(6,3); zeros(3,9)];
            
            
            
            end

            c1 = Wgt*Jdet;        

            % In rknifv.f, the shape function derivatives are computed
            % using ce_n1 which according to drive_eps_.../dupstr_blocked
            % is the current coordinates
            ElemF(1:ndf*nel) = ElemF(1:ndf*nel) - c1*Bmat'*(sigma2);
            
            ElemK(1:ndf*nel,1:ndf*nel) = ElemK(1:ndf*nel,1:ndf*nel) + (Bmat'*Smat*Bmat);
            
        end %je
        ElemK;
        
        %% Handle rest of integration points
        for ll = 2:lint

            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet,bbb,xs] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet,bbb,xs] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            % Step 1: Deformation gradient (1.143) and (1.144)
            [F,JxX,fi,Qxy] = kine3d(QXY,ul(:,1:nel),nel,1); %F_n+1
            Jdet = Jdet*JxX;
            ulm = 1/2*(ul_n + ul);
            [F2,JxX2,fi2,Qxy2] = kine3d(QXY,ulm(:,1:nel),nel,1); %F_n+1/2
            
            % Step 2: Compute polar decomposition (1.145) and (1.146)
            [Rmat, Umat] = poldec(F);
%             [Rmat, Umat_i] = rtcmp1(F);
            [R2, U2] = poldec(F2);
            
            % Form B matrix
            for mm = 1:nel  
% shape functions
 Nmat(:,3*mm-2:3*mm) = [shl(mm,1)     0          0
                           0        shl(mm,1)     0
                           0            0       shl(mm,1) ];
% derivatives w.r.t. x_n+1^i
 Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                        0         Qxy(mm,2) 0         
                        0         0         Qxy(mm,3) 
                        Qxy(mm,2) Qxy(mm,1) 0         
                        0         Qxy(mm,3) Qxy(mm,2) 
                        Qxy(mm,3) 0         Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1) 0         
                        0         Qxy(mm,3) -Qxy(mm,2) 
                        -Qxy(mm,3) 0         Qxy(mm,1) ];
% derivatives w.r.t. x_n+1/2^i
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
            % B-bar modification would follow right HERE
            

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

            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(je,lint,edge);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(je,lint,5);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY,shgs,Jdet,be,sx] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
          
            %Evaluate tangent and normal vectors
            t1 = sx(:,2);
            [tm1, tu1] = VecNormalize(t1);
            t2 = sx(:,1);
            [tm2, tu2] = VecNormalize(t2);
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            
%             if exist('iprob','var')
%             else
                Traction = traction;
%             end

            %Force components are positive in positive coord. direction
            c1 = Wgt*tm3;
            for o=1:nel
                don = shl(o);
                F = don*Traction';  %traction for the one without body force

                ElemF(ndf*o-2) = ElemF(ndf*o-2) + F(1)*c1;

                ElemF(ndf*o-1) = ElemF(ndf*o-1) + F(2)*c1;

                ElemF(ndf*o-0) = ElemF(ndf*o-0) + F(3)*c1;

            end %o
            
        end %je
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);
        
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);
        
        I1 = [1; 1; 1; 0; 0; 0];
        spvec0 = I1;
        
        thick = 1;
        
        %%%%% NOTE: These will not be correct for history dependent
        %%%%% materials since we are pulling from the wrong Gauss point
        
        % Load Guass Integration Points

        if nel == 4
            lint = 1;
            nint = 1;
        elseif nel == 8
%             lint = 4;
            lint = 8;
            nint = 8;
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
              [QXY, shgs, Jdet,bbb,xs] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,nint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet,bbb,xs] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            % Step 1: Deformation gradient (1.143) and (1.144)
            [F,JxX,fi,Qxy] = kine3d(QXY,ul(:,1:nel),nel,1); %F_n+1
            Jdet = Jdet*JxX;
            ulm = 1/2*(ul_n + ul);
            [F2,JxX2,fi2,Qxy2] = kine3d(QXY,ulm(:,1:nel),nel,1); %F_n+1/2
            
            % Step 2: Compute polar decomposition (1.145) and (1.146)
            [Rmat, Umat] = poldec(F);
%             [Rmat, Umat_i] = rtcmp1(F);
            [R2, U2] = poldec(F2);
            
            % Form B matrix
            for mm = 1:nel  
% shape functions
 Nmat(:,3*mm-2:3*mm) = [shl(mm,1)     0          0
                           0        shl(mm,1)     0
                           0            0       shl(mm,1) ];
% derivatives w.r.t. x_n+1^i
 Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                        0         Qxy(mm,2) 0         
                        0         0         Qxy(mm,3) 
                        Qxy(mm,2) Qxy(mm,1) 0         
                        0         Qxy(mm,3) Qxy(mm,2) 
                        Qxy(mm,3) 0         Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1) 0         
                        0         Qxy(mm,3) -Qxy(mm,2) 
                        -Qxy(mm,3) 0         Qxy(mm,1) ];
% derivatives w.r.t. x_n+1/2^i
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
            % B-bar modification would follow right HERE
            
            
            %% Load stresses and rotation from current converged step
            i = nh2-1+(ll-1)*nh1CP; %pointer to first history index
            history_n = hr(i+1:i+132);

            co = 76+12+max_slip_sys;
%             cn = 76+max_slip_sys+(1)*(25+max_uhard)-1;
            
            cgn1 = history_n(co:co+5);
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
        if nel == 4
            plist = [1 0 0 0
                     0 1 0 0
                     0 0 0 1];
        elseif nel == 8
            sqr3 = sqrt(3);
            plist = [-sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 0 -sqr3 sqr3 0 0 0
                     -sqr3 -sqr3 sqr3 sqr3 -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0 sqr3 0 -sqr3 -sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0
                     -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 -sqr3 -sqr3 -sqr3 -sqr3 sqr3 sqr3 sqr3 sqr3 0 0 0 0 -sqr3 sqr3 0 0 0 0 0];
%             plist = [-1 1 1 -1 -1 1 1 -1
%                      -1 -1 1 1 -1 -1 1 1
%                      -1 -1 -1 -1 1 1 1 1];
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
%%        
    case 26 % Element Stress
        
        % Mimics Mark's results output in WARP, file oust_elem.f line 132
        % First come stresses at each GP, then the averaged CP output.

        ElemS = zeros(1,nestr);
        
        % Load Guass Integration Points

        if nel == 4
            lint = 1;
        elseif nel == 8
            lint = 8;
        elseif nel == 10
            lint = 14;
        else
            lint = 27;
        end
        
        I1 = [1; 1; 1; 0; 0; 0];
        spvec0 = I1;
        str_ind = 0;
        
        der = 0;
        bf = 0;
        ib = 0;

        
        for ll = 1:lint
            
            %Evaluate first derivatives of basis functions at int. point
            if nel == 4 || nel == 10
              [Wgt,ss] =  int3d_t(ll,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet,bbb,xs] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet,bbb,xs] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
                
            % Step 1: Deformation gradient (1.143) and (1.144)
            [F,JxX,fi,Qxy] = kine3d(QXY,ul(:,1:nel),nel,1); %F_n+1
            Jdet = Jdet*JxX;
            ulm = 1/2*(ul_n + ul);
            [F2,JxX2,fi2,Qxy2] = kine3d(QXY,ulm(:,1:nel),nel,1); %F_n+1/2
            
            % Step 2: Compute polar decomposition (1.145) and (1.146)
            [Rmat, Umat] = poldec(F);
%             [Rmat, Umat_i] = rtcmp1(F);
            [R2, U2] = poldec(F2);
            
            % Form B matrix
            for mm = 1:nel  
% shape functions
 Nmat(:,3*mm-2:3*mm) = [shl(mm,1)     0          0
                           0        shl(mm,1)     0
                           0            0       shl(mm,1) ];
% derivatives w.r.t. x_n+1^i
 Bmat(:,3*mm-2:3*mm) = [Qxy(mm,1) 0         0         
                        0         Qxy(mm,2) 0         
                        0         0         Qxy(mm,3) 
                        Qxy(mm,2) Qxy(mm,1) 0         
                        0         Qxy(mm,3) Qxy(mm,2) 
                        Qxy(mm,3) 0         Qxy(mm,1) 
                        Qxy(mm,2) -Qxy(mm,1) 0         
                        0         Qxy(mm,3) -Qxy(mm,2) 
                        -Qxy(mm,3) 0         Qxy(mm,1) ];
% derivatives w.r.t. x_n+1/2^i
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
            % B-bar modification would follow right HERE
            
            
            %% Load stresses and rotation from current converged step
            i = nh2-1+(ll-1)*nh1CP; %pointer to first history index
            history_n = hr(i+1:i+hist_sz);

            co = 76+12+max_slip_sys;
%             cn = 76+max_slip_sys+(1)*(25+max_uhard)-1;
            
            cgn1 = history_n(co:co+5);
            qn1 = getrm1(Rmat,2);
            sigma = qn1*cgn1(1:6);
            
            for stres = 1:11
                
                str_ind = str_ind + 1;
            
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
            
            ElemS(str_ind) = sigmas;
            
            end
            
            % Pad a zero
            str_ind = str_ind + 1;
            
        end
        
        
%% Now do CP results averaged over the element

        CPresults = zeros(38+6+5+2+3,1);

        local_work = setup_mm10_rknstr( 1, mateprop.cp_prop, mateprop.cp_other, elem, dtWARP);
        ncrystals = 1;
        gp_temp_inc = zeros(lint,1);
        gp_temps = mateprop.cp_allTemps*ones(lint,1); % technically not used here

        % Copy history data for all integration points into local_work
        for ll = 1:lint
            
            i = nh2-1+(ll-1)*nh1CP; %pointer to first history index
            history_n = hr(i+1:i+hist_sz);
            local_work.elem_hist(1,1:hist_sz,ll) = history_n;
            
        end
        
%             matnum = iprops(38, felem)
%             if (imatprp(104,matnum) .eq. 1) then
              cnum = 1;
%             elseif (imatprp(104,matnum) .eq. 2) then
%               osn = data_offset(elnum)
%               cnum = crystal_input(osn,1)
%             else
%               write(*,*) "Invalid crystal" 
%               call die_gracefully
%             end if
% c
            nslip = mateprop.cp_prop(cnum).nslip;
% c           Why did we have to do this as a ptr?
%             history_dump = reshape(history_blocks(i).ptr,
%      &           (/history_blk_list(i),int_points,elblks(0,i)/))
% c
% c           Result 1: nslip as a double
            CPresults(1) = nslip;
% c
% c           Results 2 - 25: the slip totals, padded with zeros as required
% c                 and averaged over Gauss points
%             CPresults(2:25) = 0.0;
            s = 76+12;
            e = 76+12+nslip-1;
            CPresults(2:(2+nslip-1)) = ...
               sum(local_work.elem_hist(1,s:e,1:lint),3)/lint;
% c           Results 26-34: Nye tensor, averaged over gauss points
            s = 37;
            e = 63;
            gradFe = reshape(sum(local_work.elem_hist(1,s:e,1:lint),3)/ ...
                  lint, 3,3,3);
            nye = zeros(3,3);
            for d = 1:3
              for z = 1:3
                for b = 1:3
                  for w = 1:3
                    l_c = (z-b)*(b-w)*(w-z)/2;
                    nye(d,z) = nye(d,z) - l_c*gradFe(d,b,w);
                  end %do
                end %do
              end %do
            end %do
            CPresults(26:34) = reshape(nye, 9,1);

% c           Results 35-37: 12 tau_tilde of first crystal, avged
            s = 76+12+max_slip_sys+30+max_slip_sys;
            e = 76+12+max_slip_sys+30+max_slip_sys+11;
            CPresults(35:46) = sum(local_work.elem_hist(1,s:e,1:lint),3)/ ...
                  lint;
% c           Results 38-40: euler angles of the first crystal, averaged over
% c           GPs
            s = 76+12+max_slip_sys+6;
            e = 76+12+max_slip_sys+8;
            CPresults(47:49) = sum(local_work.elem_hist(1,s:e,1:lint),3)/ ...
                  lint;
% c           Results 41-46: elastic strain of the first crystal, averaged over
% c           GPs
            s = 76+12+max_slip_sys+24;
            e = 76+12+max_slip_sys+29;
            CPresults(50:55) = sum(local_work.elem_hist(1,s:e,1:lint),3)/ ...
                  lint;
% c           Results 47-51: first 5 uhard variables
            s = 76+12+max_slip_sys+30+max_slip_sys+max_uhard;
            e = 76+12+max_slip_sys+30+max_slip_sys+max_uhard+4;
            CPresults(56:60) = sum(local_work.elem_hist(1,s:e,1:lint),3)/ ...
                  lint;
% c           Results 52-54: first 3 tt_rates 
            s = 76+12+max_slip_sys+30+max_slip_sys+2*max_uhard;
            e = 76+12+max_slip_sys+30+max_slip_sys+2*max_uhard+11;
            CPresults(61:72) = sum(local_work.elem_hist(1,s:e,1:lint),3)/ ...
                  lint;

%           if (oubin) then
%             write(fileno) enum, p_el_type, (real(results(k)),k=1,width)
%           end if
% 
%           if (ouasc) then
% c           There's a method to my madness, casting as a real avoids
% c           3 digit exponents
%             write(fileno,'(2I8,/,(6E13.6))') enum, p_el_type, 
%      &            (real(results(k)),k=1,width)
%           end if
            ElemS(str_ind+1:str_ind+38+6+5+2+3+9+9) = CPresults;

    case 24
        
            
    case 40 % Initialize history terms
        
        if nel == 4
            lint = 1;
        elseif nel == 8
            lint = 8;
        elseif nel == 10
            lint = 14;
        else
            lint = 27;
        end
        
        local_work = setup_mm10_rknstr( 1, mateprop.cp_prop, mateprop.cp_other, elem, dtWARP); % dt = 1.0e+10;
%         
        % Loop over integration points
        for ll = 1:lint
                
            % Store history variables
            i = nh1-1+(ll-1)*nh1CP; %pointer to first history index
            history_n = hr(i+1:i+hist_sz);
            history_n = mm10_initialize_history(ll, 1,1,history_n', local_work);
            hr(i+1:i+hist_sz) = history_n;

        end %je
        
end