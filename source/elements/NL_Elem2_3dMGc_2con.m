% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 10/01/2013
% cylindrical growth model
% second version, using updt_cyl for the material model
% Verified for isotropic growth against the other implementations
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
        
        nh1 = 1*8;%nen;
%%
    case 3
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        Nmat = zeros(3,nst);
        Bmat2 = zeros(9,nst);

        emod = mateprop(1);   nue = mateprop(2);
%         mo  = mateprop(10:12); %fiber orientation vector
%         [~, mo] = VecNormalize(mo); % ensure mo is a unit vector
%         mo = mo';
        kIF = mateprop(13);
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
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgtt(xl(:,1:nel),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgb(xl(:,1:nel),nel,shld,shls,nel,bf,der,be);
            end
                
            [F,JxX,fi] = kine3d(QXY,ul(:,1:nel),nel,0);
                
            c1 = Wgt*Jdet;
            
            % Material update to find F=Fe*Fg, isotropic growth stretch
            % ratio theta=var+1
            i_var = hr(nha-1+ll);
            [var,Pe,Aeg] = updt_cyl(F,i_var,mateprop,ndm); % cylindrical
%             [Aeg,Pe,var]=cnst_vol(F,i_var,mateprop,ndm); % isotropic
            hr(nhb-1+ll) = var;
            
            Aeg_mat = reshape(Aeg,ndm*ndm,ndm*ndm);
            Pe_vec = reshape(Pe,ndm*ndm,1);
                
            for mm = 1:nel    
 Bmat2(:,3*mm-2:3*mm) = [QXY(mm,1) 0         0         
                         0         QXY(mm,1) 0         
                         0         0         QXY(mm,1) 
                         QXY(mm,2) 0         0         
                         0         QXY(mm,2) 0
                         0         0         QXY(mm,2) 
                         QXY(mm,3) 0         0         
                         0         QXY(mm,3) 0
                         0         0         QXY(mm,3) ];
 Nmat(:,3*mm-2:3*mm) = [shl(mm) 0         0        
                        0         shl(mm) 0        
                        0         0         shl(mm)];
            end
                
            % Interactive force
            
            InterForce3
            % NOTE: kIF may get reset within this routine, since
            % interaction is between a pair, not a single material
            
            % Note Also: this means that to really compute the interactive
            % force, I should do something like DG and have bigger elements
            
            ElemF = ElemF - c1*(Bmat2'*(Pe_vec) + Nmat'*IntFor);
            
            ElemK = ElemK + c1*(Bmat2'*Aeg_mat*Bmat2 + JxX*Nmat'*kIF/tstep*Nmat);
           
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
              [QXY,shgs,Jdet,be,sx] = shgtt(xlv,nel,shld,shls,nel,bf,der,be);
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
%%
    case 6
        
        ElemF = zeros(nst,1);
        Nmat = zeros(3,3*nel);
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
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgtt(xl(:,1:nel),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [QXY, shgs, Jdet] = shgb(xl(:,1:nel),nel,shld,shls,nel,bf,der,be);
            end
                
            [F,JxX,fi] = kine3d(QXY,ul(:,1:nel),nel,0);
                
            c1 = Wgt*Jdet;
            
            % Material update to find F=Fe*Fg, isotropic growth stretch
            % ratio theta=var+1
            i_var = hr(nha-1+ll);
            [var,Pe,Aeg] = updt_cyl(F,i_var,mateprop,ndm); % cylindrical
%             [Aeg,Pe,var]=cnst_vol(F,i_var,mateprop,ndm); % isotropic
            hr(nhb-1+ll) = var;
            
            Pe_vec = reshape(Pe,ndm*ndm,1);
                
            for mm = 1:nel    
 Bmat2(:,3*mm-2:3*mm) = [QXY(mm,1) 0         0         
                         0         QXY(mm,1) 0         
                         0         0         QXY(mm,1) 
                         QXY(mm,2) 0         0         
                         0         QXY(mm,2) 0
                         0         0         QXY(mm,2) 
                         QXY(mm,3) 0         0         
                         0         QXY(mm,3) 0
                         0         0         QXY(mm,3) ];
 Nmat(:,3*mm-2:3*mm) = [shl(mm) 0         0        
                        0         shl(mm) 0        
                        0         0         shl(mm)];
            end
                
            % Interactive force
            
            InterForce3
            % NOTE: kIF may get reset within this routine, since
            % interaction is between a pair, not a single material
            
            % Note Also: this means that to really compute the interactive
            % force, I should do something like DG and have bigger elements
            
            ElemF(1:ndf*nel) = ElemF(1:ndf*nel) - c1*(Bmat2'*(Pe_vec) + Nmat'*IntFor);
            if ll == 1 % Only do for the first integration point
            IForceList(1:3,elem,step) = IntFor;
            end            
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

        ElemS = zeros(nel,numstr+1);
        ElemS2 = zeros(nel,numstr);

        emod = mateprop(1);   nue = mateprop(2);
        mo  = mateprop(10:12); %fiber orientation vector
        [~, mo] = VecNormalize(mo); % ensure mo is a unit vector
        mo = mo';
        xlm = getlam([1 2 0 emod nue]);
        
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
              [Qxy,shgs,Jdet,be] = shgtt(xl(:,1:nel)+ul(1:3,1:nel),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,nint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl(:,1:nel)+ul(1:3,1:nel),nel,shld,shls,nel,bf,der,be);
            end
            
            if nelP == 4 || nelP == 10
              shlp = shltt(ss,nelP,nel,0,0);
            else
              shlp = shlb(ss,nelP,nel,0,0);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul(1:3,1:nel),nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;            
            
            % Material update to find F=Fe*Fg, isotropic growth stretch
            % ratio theta=var+1
            i_var = hr(nha-1+ll);
            [var,Pe,Aeg] = updt_cyl(F,i_var,mateprop,ndm); % cylindrical
%             [Aeg,Pe,var]=cnst_vol(F,i_var,mateprop,ndm); % isotropic
            hr(nhb-1+ll) = var;
            
            delta  = eye(ndm); 
            theta  = var(1) + 1;
%             Fg_inv = (1/theta)*delta;
            Fg_inv = (1/theta)*delta + (1-1/theta)*(mo*mo');
            Fe = F * Fg_inv;
            Fe_inv = inv(Fe);
            Je = det(Fe);
            sigma = Pe*Fe'/Je;
            sigma = [sigma(1,1); sigma(2,2); sigma(3,3); sigma(1,2); sigma(2,3); sigma(3,1); zeros(3,1)];
            
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
              [Qxy,shgs,Jdet,be] = shgtt(xl(:,1:nel)+ul(1:3,1:nel),nel,shld,shls,nel,bf,der,be);
            else
              [Wgt,ss] =  intpntb(ll,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgb(xl(:,1:nel)+ul(1:3,1:nel),nel,shld,shls,nel,bf,der,be);
            end
            
            [fi,JxX,F] = kine3d(Qxy,-ul(1:3,1:nel),nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;

            w = Wgt*Jdet*thick;
            
            Vol = Vol + w;

        end %je

        for i = 1:nel
        ElemS(i,numstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
        
    case 40 % Initialize Intermediate Configuration
        
        
        if exist('numconst','var') && numconst==3 && elem > numel_mix/3 && elem <= 2*numel_mix/3 %ma==2
            % For 3 constituent model, start the middle constituent at a
            % grown state.
            
            lint = 8;

            % Loop over integration points
            for l = 1:lint

                % Store history variables
                hr(nhb-1+l) = 0.0030;
                hr(nha-1+l) = 0.0030;

            end %je
            
        end
end