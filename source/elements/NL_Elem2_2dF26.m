function ElemS = NL_Elem2_2dF26(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,nestr)
        
% Purpose: Compute element stresses, e.g. at integration points
% Called by: FormFE.m
% Notes: Number of stresses set as nestr in pmatin; should be e.g. number
%        of integration points times number of stresses per point
% Example: NL_Elem2_2dM.m

        ElemS = zeros(1,nestr);

        lam = getlam(mateprop);
        
        I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
        stresID = [1 2 4 0 1 3];
        spvec0 = I1;
        
        % Load Guass Integration Points

        
        der = 0;
        bf = 0;
        ib = 0;


            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
              ss = 1/3*[1 1];
              [shl,shld,shls,be] = shlt(ss(1),ss(2),nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgt(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            else
              ss = [0 0];
              [shl,shld,shls,be] = shlq(ss(1),ss(2),nel,nel,der,bf);
              [Qxy,shgs,Jdet,be] = shgq(xl+ul(1:2,:),nel,shld,shls,nel,bf,der,be);
            end
            
            if nelP == 3 || nelP == 6
              shlp = shlt(ss(1),ss(2),nelP,nel,0,0);
            else
              shlp = shlq(ss(1),ss(2),nelP,nel,0,0);
            end
            
            [fi,JxX,F] = kine2d(Qxy,-ul,nel,0); %this is equivalent to ikine2d
            JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            Jdet = Jdet/JxX;
            
            sigma = SigmaCmat3([F zeros(2,1); zeros(1,2) 1],JxX,mateprop,lam);
            sigma = sigma/JxX;
            
            for stres = 1:nestr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                sigma2 = [sigma(1) sigma(4) sigma(6); sigma(4) sigma(2) sigma(5); sigma(6) sigma(5) sigma(3)];
                psig = eig(sigma2);
                sigmas = psig(stresID(stres));
                else % hydrostatic stress
                sigmas = 1/3*sigma'*I1;
                end
            else % von Mises stress
                trs = sigma'*I1;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS(1,stres) = sigmas;
            
            end