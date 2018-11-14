% Tim Truster
% 08/14/2014
%
% Crystal plasticity element based on Mark's first implementation
% 2D version for plane-strain analyses
%
%
% 04/19/2015 - Omar's try to modify the routine to include subcycling,
% tested against MarkCP2D and give same results
[max_slip_sys,max_uhard] = maxparamsCP;
vec1 = zeros(max_uhard,1);
vec2 = zeros(max_uhard,1);
if ~(exist('tstep','var')) % set dt for CP material routine to tstep from input file
    dtWARP = 1.0e+10;
else
    dtWARP = tstep;
end
if ~(exist('initcycle','var'))
    initcycle = 1;2; % Way to initialize Rt, gradFe values for subcycling;
end
% initcycle=1 is the old way, where the previous step value (n) is used,
% while initcycle=2 uses the converged value of Rt from the previous
% iteration as an initial guess for the new internal variables at step
% n+1. This seems to work pretty well.
if numsubcycles == 0 % Iterations for implicit exponential integration of Rp
    numReiter = 0; % Explicit
else
    numReiter = 20; % Implicit
end
RpTol = 1e-12; % Rp is orthogonal, so its norm should be well conditioned and not need a ratio
if ~(exist('ReMeth','var'))
    ReMeth = 1;0; % 1 for Lie, 0 for Mark;
end

switch isw %Task Switch
    %%
    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        hist_sz = 76+12+max_slip_sys+(30+max_slip_sys+3*max_uhard);
        nh1CP = hist_sz;
        if nen == 6
        nh1 = hist_sz*7;
        nh3CP = 9+27+27+27; % Re, actual grad to use, Mark's, smooothed
        nh3 = nh3CP*7;
        else
        nh1 = hist_sz*nen;
        nh3CP = 9+27+27+27; % Re, actual grad to use, Mark's, smooothed
        nh3 = nh3CP*nen;
        end
        % Sets number of Lie variables for interpolation, later
        numLie = 18; % tensors Re and Ue
        istv = 11; % number of stresses per node for post-processing
        iste = 11*nen+38+6+5+2+3; % number of stresses per element (total number, including all integration points)
        if nen == 10
            iste = iste + 1*11;
        end
        
        %%
    case {3,6}
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Bmat = zeros(4,2*nel);
        Bmat2 = Bmat;
        
        %         % Use this to test whether adaptive time stepping works
        %         if step == 1 && iter == 1 && (numincrem1 == 1 || numincrem2 == 1)
        %             failFEelem = 1;
        %             return
        %         end
        
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
        ngp = lint;
        der = 0;
        bf = 0;
        ib = 0;
        
        uldres = reshape(uld,nst,1);
        if step == 30
            step;
        end
        
        
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
        if initcycle == 2
            subcycle2 = (initia==0) + 1; % only copy over or use the values from the previous time step for the initial prediction iteration; after that, use values from last iteration.
        else
            subcycle2 = subcycle; % always use the values from the previous iteration as the initial guess.
        end
        local_work.step = stepW;
        
        % Copy history data for all integration points into local_work
        for ll = 1:lint
            
            i = nh1-1+(ll-1)*nh1CP; %pointer to first history index
            history_n = hr(i+1:i+hist_sz);
            
            if subcycle2 == 1 % Copy history of Re, gradRe into the holding spot
                
                co = 76+12+max_slip_sys;
                Rp = history_n(co-1+10:co-1+18);
                gradFe = history_n(37:63);
                backstress = history_n(76:76+11);
                j = nh3-1+(ll-1)*nh3CP; %pointer to first history index
                initialRegrad = [Rp; gradFe];
                hr(j+1:j+36) = initialRegrad;
            end
            
            % Extract Re, gradRe from holding spot and into Mark's arrays
            % This way, I didn't have to modify his mm10 routines very much
            j = nh3-1+(ll-1)*nh3CP; %pointer to first history index
            Regrad = hr(j+1:j+36);
            co = 76+12+max_slip_sys;
            history_n(co-1+10:co-1+18) = Regrad(1:9); % Use latest value for Rp_n+1, which would be coming from the last subcycle
            history_n(37:63) = Regrad(10:36);

            local_work.elem_hist(1,1:hist_sz,ll) = history_n;
            
            
        end
        
        
        %% Integration Loop
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
            if initia == 1 % Follow lnstiff.f version
                
                
                % Load stresses and rotation from last converged step
                i = nh1-1+(ll-1)*nh1CP; %pointer to first history index
                history_n = hr(i+1:i+hist_sz);
                
                co = 76+12+max_slip_sys;
                %                 cn = 76+12+max_slip_sys+(1)*(25+max_uhard)-1;
                
                cgn1 = history_n(co:co+5);
                
                
                if step == 1 % set initial stiffness as in lnstff.f, line 646
                    
                    qn1 = getrm1(Rmat,2);
                    sigma = qn1*cgn1(1:6);
                    sigma2 = [sigma; 0; 0; 0]; % get stress from last time step
                    sigma2 = sigma2([1 2 4 7]);
                    
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
                    cep = cep([1 2 4],[1 2 4]);
                    
                    Smat = [cep zeros(3,1); zeros(1,4)];
                    
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
                local_work.jac(1,:,:) = [xs zeros(2,1); [zeros(1,2) 1]];
                uddt = reshape(deps,1,6);
                local_work.rot_blk_n1(1, 1:9, ll) = reshape(Rmat,1,9);
                
                if subcycle > 1
                    % Copy values of Cauchy stress, tau_tilde from previous
                    % subcycle in order to compute the norm of the difference
                    % with the current subcycle
                    co = 76+12+max_slip_sys;
                    cn = 76+12+max_slip_sys+(30+max_slip_sys+3*max_uhard)-1;
                    i = nh2-1+(ll-1)*nh1CP; %pointer to first history index
                    oldHist = hr(i+co:i+cn);
                    cgn_old = oldHist(1:6);
                    tt_old = oldHist(30+max_slip_sys+1:30+max_slip_sys+local_work.c_props(1,1).num_hard);
                    qn1 = getrm1(Rmat,2);
                    sigma_old = qn1*cgn_old(1:6);
                end
                
                [local_work.elem_hist1(1,1:hist_sz,ll),local_work.elem_hist(1,1:hist_sz,ll), ...
                    local_work, gp_temps, gp_temp_inc] = mm10(ll, 1,...
                    ncrystals, hist_sz, local_work.elem_hist(1,1:hist_sz,ll), local_work, uddt,...
                    gp_temps, gp_temp_inc, 1,subcycle);
                
                if local_work.material_cut_step
                    failFEelem = 1;
                    return
                end
                
                
                %% Implicit exponential update to rotation matrix
                % Get Rp_n value from history array, from previous step n
                i = nh1-1+(ll-1)*nh1CP; %pointer to first history index
                history_n = hr(i+1:i+hist_sz);
                co = 76+12+max_slip_sys;
                Rp_n = reshape(history_n(co-1+10:co-1+18),3,3);
                % Move current guess for Rp_n+1 into np1 vector
                % technically using the one computed inside mm10, which won't
                % be right anymore, but this is an initial guess.
                Rp_n1 = reshape(local_work.elem_hist1(1,co-1+10:co-1+18,ll),3,3);
                % Re-create data structures that occur inside mm10 at the time
                % when the routine converged to a stress state
                np1 = mm10_setup_np1...
                    (local_work.rot_blk_n1(1,1:9,ll), ...
                    uddt(1,1:6), local_work.dt, gp_temps(1), ...
                    local_work.step,1+local_work.felem, ll);
                np1.dg = sqrt(2.0/3.0*(dot(np1.D(1:3),np1.D(1:3))...
                    +0.5*dot(np1.D(4:6),np1.D(4:6))));
                np1.stress = local_work.urcs_blk_n1(1,1:6,ll);
              np1.tau_tilde = local_work.elem_hist1(1,co-1+30+max_slip_sys+1:co-1+30+max_slip_sys+local_work.c_props(1,1).num_hard,ll);
                % above lines are from mm10, mm10_setup, mm10_store_gp, and mm10_store_cryhist
                np = crystal_state;
                np1.Rp = Rp_n1;
                np.Rp = Rp_n;
                % Load ms and qs arrays
                cc_props= mm10_init_cc_props(local_work.c_props(1,1)...
                    ,local_work.angle_type(1),...
                    local_work.angle_convention(1));
                Rp_n1_old = Rp_n1;
                % Loop to satisfy the nonlinear equation:
                % Rp_n+1 = exp(A(Rp_n+1))*Rp_n
                for impReiter = 1:numReiter
                    
                    % Bring forward the lattice arrays to most current
                    % configuration, as they are involved in computing the A
                    % tensor in the exponential.
                    % Technically the stress and hardening variables are too, but
                    % they are taken care of in mm10, where the new Rp_n+1 is used.
                    RE = mm10_RT2RVE(transpose(np1.Rp));
                    RW = mm10_RT2RVW(transpose(np1.Rp));
                    RWC = ...
                        mm10_RT2RVW(np1.R*transpose(np1.Rp));
                    for i=1:cc_props.nslip
                        np1.ms(1:6,i) = RE*cc_props.ms(1:6,i);
                        np1.qs(1:3,i) = RW*cc_props.qs(1:3,i);
                        np1.qc(1:3,i) = RWC*cc_props.qs(1:3,i);
                    end
 [cc_props, np1, np, np1.stress, np1.tau_tilde, vec1,vec2] = mm10_formvecs(cc_props,np1,np,np1.stress,... 
 np1.tau_tilde,vec1,vec2);
              [cc_props, np1, np,vec1,vec2] = mm10_update_rotation(cc_props, np1, np,vec1,vec2);
                    Rp_n1 = np1.Rp; % Just for seeing the value during debug
                    if norm(Rp_n1-Rp_n1_old) < RpTol
                        %                   impReiter
                        break
                    else
                        Rp_n1_old = Rp_n1;
                    end
                end
                
                % Copy new Rp_n+1 from np1 into hr3 and hr2
                Rp_n1 = np1.Rp;
                local_work.elem_hist1(1,co-1+10:co-1+18,ll) = reshape(Rp_n1,9,1);
                if impReiter > 0 % hack to make sure explicit version still works
                    j = nh3-1+(ll-1)*nh3CP; %pointer to first history index
                    hr(j+1:j+9) = reshape(Rp_n1,9,1);
                end
                %             if ll == 1 && step > 2
                %                 Rp_n1
                %             end
                
                
                %             % copy step n history only if first step
                %             if stepW == 1
                %             i = nh1-1+(ll-1)*nh1CP; %pointer to first history index
                %             hr(i+1:i+132) = local_work.elem_hist(1,1:132,ll);
                %             end
                
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
                
                
                co = 76+12+max_slip_sys;
            tautilde = local_work.elem_hist1(1,co-1+30+max_slip_sys+1:co-1+30+max_slip_sys+local_work.c_props(1,1).num_hard,ll);
                if subcycle > 1 % Compute subcycle changes to tau_tilde and stress
                    SubNorms(2) = SubNorms(2) + norm(sigma-sigma_old);
                    SubNorms(3) = SubNorms(3) + norm(tautilde-tt_old);
                else % Store the initial norm of the stresses for scaling tolerances
                    SubNorms(2) = SubNorms(2) + norm(sigma);
                    SubNorms(3) = SubNorms(3) + norm(tautilde);
                end
                
                
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
                local_work2 = cnst10(ll, 1, local_work.elem_hist1(1,1:132,ll), local_work2);
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
        
        
        %% Compute gradFeInv according to Mark's method
        if ~initia
            
            % Don't do FeGrad for ORNL model
            if local_work.c_props.h_type ~= 4 || local_work.c_props.h_type == 4
                
                local_work.span = 1;
                % c  In rknstr, right after rstgp1, line 451
                % c           For CP model, calculate the gradient of the elastic
                % c           rotations at the element level by linear curve fit.
                % c
                % c           For linear models this will just be based on the plastic rotations
                % c           which may or may not be a realistic assumption
                % c
                %       if (local_work%mat_type .eq. 10) then
                for i=1:local_work.span
                    if (local_work.ncrystals(i) > 1) %then
                        local_work.elem_hist1(i,37:63,1:ngp) = 0.0;
                    else
                        %              rs = 76+max_slip_sys+9;
                        %              re = 76+max_slip_sys+17;
                        %              Rps = squeeze(local_work.elem_hist(i,rs:re,1:ngp));
                        Rps = zeros(9,lint);
                        for ll = 1:lint
                            j = nh3-1+(ll-1)*nh3CP; %pointer to first history index
                            Rps(:,ll) = hr(j+1:j+9);
                        end
                        local_work.elem_hist1(i,37:63,1:ngp) = mm10_calc_grads2(ngp, 1, 1, 1, ...
                            squeeze(local_work.rot_blk_n1(i,1:9,1:ngp)), ... % R tensor from polar decomposition F=RU (namely R = F*Uinv)
                            squeeze(local_work.jac(i,1:3,1:3)), ... % element jacobian at last integration point, w.r.t. reference coords, according to gtmat1.f
                            Rps,nel); %Rps, stored in mm10_store_cryhist line 25 from np1.Rp
                    end %if
                end %do
                %       end if
                % Copy gradFeInv back into hr history array
                for ll = 1:lint
                    i = nh2-1+(ll-1)*nh1CP; %pointer to first history index
                    hr(i+1:i+hist_sz) = local_work.elem_hist1(1,1:hist_sz,ll);
                    % Also put gradFe into hr3 part for use in next subcycle
                    i = nh3-1+(ll-1)*nh3CP; %pointer to first history index
                    hr(i+37:i+36+27) = local_work.elem_hist1(1,37:63,ll);
                end
                
            end
            
        end
        
        %%
    case -1
        
        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS
        
        % Determine bounds of integration
        
        if nel == 4 || nel == 9
            
            dr = 2;
            ro = -1;
            
            % Upper Limit
            if nodeA == ElemFlag(2)
                eR2 = 1;
            else %
                eR2 = 0;
            end
            % Lower Limit
            if nodeB == ElemFlag(1)
                eR1 = -1;
            else %nodeA == ElemFlag(5)
                eR1 = 0;
            end
            
        elseif nel == 3 || nel == 6
            
            dr = 1;
            ro = 0;
            
            % Upper Limit
            if nodeA == ElemFlag(2)
                eR2 = 1;
            else %nodeA == ElemFlagR(5)
                eR2 = 1/2;
            end
            % Lower Limit
            if nodeB == ElemFlag(1)
                eR1 = 0;
            else %nodeA == ElemFlag(5)
                eR1 = 1/2;
            end
            
        end
        
        % Set jacobian for integration space
        drdr = (eR2 - eR1)/dr;
        
        lint = 4;
        
        ideriv = 0;
        der = 0;
        bf = 0;
        thick = 1;
        
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
        
        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);
        
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
            
            
            %% Load stresses and rotation from current converged step
            i = nh2-1+(ll-1)*nh1CP; %pointer to first history index
            history_n = hr(i+1:i+hist_sz);
            
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
        
        % Mimics Mark's results output in WARP, file oust_elem.f line 132
        % First come stresses at each GP, then the averaged CP output.
        
        ElemS = zeros(1,nestr);
        
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
        
        I1 = [1; 1; 1; 0; 0; 0];
        spvec0 = I1;
        str_ind = 0;
        
        der = 0;
        bf = 0;
        ib = 0;
        
        
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
        
        CPresults = zeros(38+6+5,1);
        
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
        
        % c           Results 35: tau_tilde of first crystal, avged
        s = 76+12+max_slip_sys+30+max_slip_sys;
        CPresults(35) = sum(local_work.elem_hist(1,s,1:lint))/ ...
            lint;
        % c           Results 36-38: euler angles of the first crystal, averaged over
        % c           GPs
        s = 76+12+max_slip_sys+6;
        e = 76+12+max_slip_sys+8;
            CPresults(36:38) = sum(local_work.elem_hist(1,s:e,1:lint),3)/ ...
                  lint;
        % c           Results 39-44: elastic strain of the first crystal, averaged over
        % c           GPs
        s = 76+12+max_slip_sys+24;
        e = 76+12+max_slip_sys+29;
        CPresults(39:44) = sum(local_work.elem_hist(1,s:e,1:lint),3)/ ...
              lint;
        % c           Results 45-49: first 5 uhard variables
        s = 76+12+max_slip_sys+30+max_slip_sys+max_uhard;
        e = 76+12+max_slip_sys+30+max_slip_sys+max_uhard+4;
        CPresults(45:49) = sum(local_work.elem_hist(1,s:e,1:lint),3)/ ...
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
        ElemS(str_ind+1:str_ind+38+6+5) = CPresults;
        
    case 24
        
        
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
        
    case 100 % Nye tensor field post-processing, subcycle
        
        if nel == 3
            lint = 1;
        elseif nel == 4
            lint = 4;
        elseif nel == 6
            lint = 7;
        else
            lint = 9;
        end
        
        ElemL = zeros(nel,numLie+1);
        ElemL2 = zeros(numLie,lint);
        
        local_work = setup_mm10_rknstr( 1, mateprop.cp_prop, mateprop.cp_other, elem, dtWARP); % dt = 1.0e+10;
        ncrystals = 1;
        gp_temp_inc = zeros(lint,1);
        gp_temps = mateprop.cp_allTemps*ones(lint,1); % technically not used here
        local_work.step = step;
        
        
        % Loop over nodes, convert Lie group to Lie algebra
        for ll = 1:lint
            
            i = nh3-1+(ll-1)*nh3CP; %pointer to first history index
            history_n = hr(i+1:i+36);
            
            
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
            
            local_work.jac(1,:,:) = [xs zeros(2,1); [zeros(1,2) 1]];
            local_work.rot_blk_n1(1, 1:9, ll) = reshape(Rmat,1,9);
            
            rs = 76+12+max_slip_sys+9;
            re = 76+12+max_slip_sys+17;
            Rps = history_n(1:9);
            rot_blk = squeeze(local_work.rot_blk_n1(1,1:9,ll));
            
            if (local_work.geo_non_flg)
                jacinv = squeeze(local_work.jac(1,:,:));
                jacinv = inv(jacinv);
                Rt_grp = reshape(Rps(1:9),3,3)*...
                    reshape(rot_blk(1:9),3,3)';
            else
                Rt_grp = reshape(Rps(1:9), 3,3);
            end
            
            % convert rotation
            Rt_alg = real(logm(Rt_grp));
            Rt_alg = 1/2*(Rt_alg - Rt_alg'); % Force to be skew
            ElemL2(1:9,ll) = reshape(Rt_alg,1,9);
            
        end
        
        
        %% Simple extrapolation of integration points to nodes and averaging
        if LieLump == 0
            
            % Extrapolation parametric locations for nodes
            if nel == 4
                lint = 4;
                nelS = 4;
                sqr3 = sqrt(3);
                plist = [-sqr3 sqr3 sqr3 -sqr3
                    -sqr3 -sqr3 sqr3 sqr3];
                %          n1    n2   n3   n4
                % Map history ordering to node ordering
                %             ElemL2 = ElemL2(:,[1 2 4 3]);
                %         elseif nel == 6
                %             lint = 7;
                %             nint = 3;
                %             plist = [-1/3 5/3 -1/3 2/3 2/3 -1/3
                %                      -1/3 -1/3 5/3 -1/3 2/3 2/3];
                %         else
                %             lint = 9;
                %             nint = 4;
                %             plist = [-sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 -sqr3 0
                %                      -sqr3 -sqr3 sqr3 sqr3 -sqr3 0 sqr3 0 0];
            end
            
            der = 0;
            bf = 0;
            
            % Extrapolate Lie algebra values from integration points to nodes
            for ll = 1:nelS
                
                ss = plist(:,ll);
                
                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                    shl = shlt(ss(1),ss(2),nel,nel,der,bf);
                else
                    shl = shlq(ss(1),ss(2),nel,nel,der,bf);
                end
                
                % Extrapolation is a product of shape functions and integration
                % point values
                Extra_vals = ElemL2*shl;
                ElemL(ll,1:numLie) = Extra_vals';
                
            end
            
            %Integration Loop for area
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
                
                w = Wgt*Jdet;
                
                Vol = Vol + w;
                
            end %je
            
            for i = 1:nel
                ElemL(i,numLie+1) = 1; % use this for simple average Vol; % use this for weighted average
            end
            
            % Multiply areas onto the nodal points
            for k = 1:nel
                ElemL(k,1:numLie) = ElemL(k,1:numLie)*ElemL(k,numLie+1);
            end
            
        else % create lumped mass matrix using HRZ method, notation from Mota paper
            
            der = 0;
            bf = 0;
            
            ElemM = zeros(nel,nel);
            
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
                
                w = Wgt*Jdet;
                
                % Sigma_beta = integral(shape_function*integ_pt_value)
                ElemL(1:nel,1:numLie) = ElemL(1:nel,1:numLie) + w*shl*ElemL2(1:numLie,ll)';
                
                % H_ab = consistent mass entry
                ElemM = ElemM + w*(shl*shl');
                
            end %je
            
            % HRZ lumping technique
            total_mass = sum(sum(ElemM));
            s_sum = sum(diag(ElemM));
            ElemL(1:nel,numLie+1) = diag(ElemM)*total_mass/s_sum;
            
        end %LieLump
        
        
    case 101 % Convert Lie algebra to Lie group, subcycle
        
        ErrorsR = zeros(3,3);
        
        % Loop over nodes, interpolate z(X), and compute Lie algebra Z = exp(z)
        for ll = 1:nel
            
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
            
            
            % Compute dZiI/dXJ = dexp(ziI)/dzkK*dzkK/dXJ
            grads = zeros(3,3,3);
            ziI = reshape(ElemL(1:nel,1:9)'*shl,3,3);
            dzkKdXJ = ElemL(1:nel,1:9)'*QXY; % inverse of Jacobian
            for L = 1:2
                dzkKdXL = reshape(dzkKdXJ(1:9,L),3,3);
                [~,dZiIdXL] = dexpmt(ziI,dzkKdXL);
                grads(1:3,1:3,L) = dZiIdXL;
            end
            
            
            i = nh3-1+(ll-1)*nh3CP; %pointer to first history index
            history_n = hr(i+1:i+nh3CP);
            gradFes = reshape(grads,27,1);
            gradFesMark = history_n(37:63); % From Mark's averaging routine
            hr(i+64:i+63+27) = gradFes; % put in a safe place
            hr(i+37:i+36+27) = gradFesMark; % copy Mark's too
            if ReMeth == 1
                hr(i+10:i+36) = gradFes; % overwrite Mark's values
                i = nh2-1+(ll-1)*nh1CP; %pointer to first history index
                hr(i+37:i+63) = gradFes; % overwrite Mark's values
            else
                hr(i+10:i+36) = gradFesMark; % use Mark's values
                i = nh2-1+(ll-1)*nh1CP; %pointer to first history index
                hr(i+37:i+63) = gradFesMark; % overwrite Mark's values
            end
            
            %             if ll == 1 && step > 2
            %                 gradFes(1:3)
            %             end
            
            % Compare Re
            
            % Step 1: Deformation gradient (1.143) and (1.144)
            [F,JxX,fi,Qxy] = kine2d(QXY,ul(:,1:nel),nel,1); %F_n+1
            Jdet = Jdet*JxX;
            ulm = 1/2*(ul_n + ul);
            [F2,JxX2,fi2,Qxy2] = kine2d(QXY,ulm(:,1:nel),nel,1); %F_n+1/2
            F = [F zeros(2,1); zeros(1,2) 1];
            
            % Step 2: Compute polar decomposition (1.145) and (1.146)
            [Rmat, Umat] = poldec(F);
            
            Rt_int = expm(ziI);
            rs = 76+12+max_slip_sys+9;
            re = 76+12+max_slip_sys+17;
            Rps = history_n(1:9);
            rot_blk = Rmat;
            
            %              if (local_work.geo_non_flg)
            Rt_grp = reshape(Rps(1:9),3,3)*rot_blk';
            %              else
            %                   Rt_grp = reshape(Rps(1:9), 3,3);
            %              end
            if ll == 1
                %                 Rt_int,Rt_grp
            end
            ErrorsR = ErrorsR +  Rt_int-Rt_grp;
            if ll == 1
                %                 ErrorsR
            end
            
        end
        %         % Need to think about and update this portion
        %
        %         Rt_alg = ElemL(1:9);
        %         Rt_alg = reshape(Rt_alg,3,3);
        %         Rt_grp = expm(Rt_alg);
        %         ElemL(1:9) = reshape(Rt_grp,1,9);
        
end