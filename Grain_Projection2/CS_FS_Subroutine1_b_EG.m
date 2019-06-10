% Two dimenisional plane stress/strain element
% Modified subroutine forsolid elements with flags for coarse and flag
% scale
% Eline Geut
% Created 3/7/2019. Last modified 3/7/2019 
PatchE = mateprop(1);
Patchv = mateprop(2);
thick = mateprop(3);
GrainA = mateprop(4);
if ~exist('PSPS','var')
PSPS = 'n'; %plane stresS/straiN flag
end
if ~exist('iprob','var')
iprob = 0;
end

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

switch isw %Task Switch
    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        istv = 7; % number of stresses per node for post-processing
        iste = 7; % number of stresses per element
%%
    case 3
 
        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        
        if PSPS == 'n'
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        else
        Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                                  Patchv  1      0
                                  0      0      (1-Patchv)/2];
        end

        der = 0;
        bf = 0;
        ib = 0; 
        ulres = reshape(ul,ndf*nen,1);
        
        if mateprop(5) == 1
            flag_CC = 1;
            flag_FF = 0;
        elseif mateprop(5) == 0
            flag_CC = 0;
            flag_FF = 1;
        else 
            disp('error, wrong flag value');
        end 

        %Coarse Scale Stiffness Martix 
        if flag_CC == 1
            K_CC = zeros(ndf*nel);
            ElemF_CC = zeros(ndf*nel,1);
            
            % Load Gauss Points for quadrature
            lint = 1;
            Bmat_CC = [0 0 1 0 0 0
                0 0 0 1 0 0
                0 0 0 0 1 0];
            elemA = GrainA/numelemg;
            c1_CC = elemA*thick;
            
            sigma_CC = Dmat*(Bmat_CC*ulres(1:nel*ndf));
            ElemF_CC = ElemF_CC - c1_CC*Bmat_CC'*sigma_CC;
            K_CC = K_CC + c1_CC*Bmat_CC'*Dmat*Bmat_CC;
            
        elseif flag_FF == 1
            %Fine Scale Stiffness Matrix
            
            ElemK = zeros(ndf*nel);
            ElemF = zeros(ndf*nel,1);
            Nmat = zeros(2,ndf*nel);
            Bmat = zeros(3,ndf*nel);
            
            % Load Gauss Points for quadrature

        K_FF = zeros(ndf*nel);
        ElemF_FF = zeros(ndf*nel,1);
        Nmat_FF = zeros(2,ndf*nel);
        Bmat_FF = zeros(3,ndf*nel);

        % Load Gauss Points for quadrature
        if nel == 3
            lint = lintt3;%13;
        elseif nel == 4
            lint = lintq4;
        elseif nel == 6
            lint = lintt6;%13;
        elseif nel == 9
            lint = lintq9;
        end
        
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(je,lint,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [Wgt,litr,lits] = intpntq(je,lint,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
                    
            c1_FF = Wgt*Jdet*thick;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat_FF(1,(ie-1)*ndf+1) = shl(ie);
              Nmat_FF(2,(ie-1)*ndf+2) = shl(ie);
                
              Bmat_FF(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat_FF(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
                 
            end
            
            xint = xl(:,1:nel)*shl;
            cvec = zeros(3,1);
            sigma_FF = Dmat*(Bmat_FF*ulres(1:nel*ndf) - cvec);
            ElemF_FF = ElemF_FF - c1_FF*Bmat_FF'*sigma_FF;
            
            K_FF = K_FF + c1_FF*Bmat_FF'*Dmat*Bmat_FF;

        end 
        end %je
        
        %Coarse-Fine stiffness matrices 
        %Still needs to be anjusted and finilized. Mainly c1 constant
        %determination needs to be worked on
        
%         K_CF = c1_CC.*Bmat_CC'*Dmat*Bmat_FF; %CF element stiffness matrix
%         K_FC = c1_FF.*Bmat_FF'*Dmat*Bmat_CC; %FC element stiffness matrix
%         
%         ElemF = [ElemF_CC; ElemF_FF];
%         %Assembly of the General Stiffness Matrix
%         ElemK = [K_CC K_CF; K_FC K_FF];
        
        %Modification of 4/27/19 to separate Meso and Micro 
        
        if flag_CC == 1
            ElemK = K_CC;
            ElemF = ElemF_CC;
        elseif flag_FF == 1
            ElemK = K_FF;
            ElemF = ElemF_FF;
        else
            disp('error, wrong flag value');
        end
     ElemK;         
  
%%       
 case 26
        %% General Values
        
        stresID = [1 2 3 0 1 3];
        %         ElemS = zeros(nestr,1);
        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        if PSPS == 'n'
            Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        else
            Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                Patchv  1      0
                0      0      (1-Patchv)/2];
        end
    
        if mateprop(5) == 1
            flag_CC = 1;
            flag_FF = 0;
        elseif mateprop(5) == 0
            flag_CC = 0;
            flag_FF = 1;
        else
            disp('error, wrong flag value');
        end
        
        %Fine Scale
        
        if flag_FF == 1
            
            Nmat_FF = zeros(2,ndf*nel);
            Bmat_FF = zeros(3,ndf*nel);
            I1 = [1; 1; 0];
            
            % Load Guass Integration Points
            
            if nel == 3
                lint = 1;
                nint = 1;
            elseif nel == 4
                %             lint = 4;
                lint = 4;
                nint = 1;
            elseif nel == 6
                lint = 7;
                nint = 3;
            else
                lint = 9;
                nint = 4;
            end
            
            der = 0;
            bf = 0;
            ib = 0;
            
            ulres = reshape(ul,nst,1);
            
            %Stress Loop
            ll = 1;
            %         for ll = 1:nint
            
            %Evaluate first derivatives of basis functions at int. point
            if nel == 3 || nel == 6
                [w,litr,lits] =  intpntt(ll,ll,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg,shgs,Jdet,be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
            else
                [w,litr,lits] =  intpntq(ll,ll,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg,shgs,Jdet,be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            % Form B matrix
            for ie = 1:nel
                
                Bmat_FF(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
                Bmat_FF(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
                
            end
            
            epsil_FF = Bmat_FF*ulres(1:ndf*nel);
            sigma_FF = Dmat*epsil_FF;
            
            for stres = 1:npstr
                
                if stres <= 3 % stress components
                    sigmas_FF = sigma_FF(stresID(stres));
                elseif stres >= 5
                    if stres <= 6 % principal stresses
                        if PSPS == 'n'
                            sigz_FF = lam*(epsil_FF(1)+epsil_FF(2));
                        else
                            sigz_FF = 0;
                        end
                        sigma2_FF = [sigma_FF(1) sigma_FF(3) 0; sigma_FF(3) sigma_FF(2) 0; 0 0 sigz_FF];
                        psig_FF = eig(sigma2_FF);
                        sigmas_FF = psig_FF(stresID(stres));
                    else % hydrostatic stress
                        if PSPS == 'n'
                            sigz_FF = lam*(epsil_FF(1)+epsil_FF(2));
                        else
                            sigz_FF = 0;
                        end
                        sigmas_FF = 1/3*(sigma_FF'*I1 + sigz_FF);
                    end
                else % von Mises stress
                    if PSPS == 'n'
                        sigz_FF = lam*(epsil_FF(1)+epsil_FF(2));
                    else
                        sigz_FF = 0;
                    end
                    trs_FF = sigma_FF'*I1 + sigz_FF;
                    dsig_FF = sigma_FF - 1/3*trs_FF*I1;
                    sigmas_FF = sqrt(3/2*(dsig_FF'*dsig_FF));
                end
                
                ElemS(stres) = sigmas_FF;
                
            end
            ElemS;
            
        elseif flag_CC == 1
            % Case 26 Coarse Scale
            ElemS_CC = zeros(nestr,1);
            Bmat_CC = [0 0 1 0 0 0
                0 0 0 1 0 0
                0 0 0 0 1 0];
            I1 = [1; 1; 0];
            
            % Load Guass Integration Points
            
            lint = 1;
            nint = 1;
            
            der = 0;
            bf = 0;
            ib = 0;
            
            %Stress Loop
            ll = 1;
            
            epsil_CC = Bmat_CC*ulres(1:ndf*nel);
            sigma_CC = Dmat*epsil_CC;
            
            for stres = 1:npstr
                
                if stres <= 3 % stress components
                    sigmas_CC = sigma_CC(stresID(stres));
                elseif stres >= 5
                    if stres <= 6 % principal stresses
                        if PSPS == 'n'
                            sigz_CC = lam*(epsil_CC(1)+epsil_CC(2));
                        else
                            sigz_CC = 0;
                        end
                        sigma2_CC = [sigma_CC(1) sigma_CC(3) 0; sigma_CC(3) sigma_CC(2) 0; 0 0 sigz_CC];
                        psig_CC = eig(sigma2_CC);
                        sigmas_CC = psig_CC(stresID(stres));
                    else % hydrostatic stress
                        if PSPS == 'n'
                            sigz_CC = lam*(epsil_CC(1)+epsil_CC(2));
                        else
                            sigz_CC = 0;
                        end
                        sigmas_CC = 1/3*(sigma_CC'*I1 + sigz_CC);
                    end
                else % von Mises stress
                    if PSPS == 'n'
                        sigz_CC = lam*(epsil_CC(1)+epsil_CC(2));
                    else
                        sigz_CC = 0;
                    end
                    trs_CC = sigma_CC'*I1 + sigz_CC;
                    dsig_CC = sigma_CC - 1/3*trs_CC*I1;
                    sigmas_CC = sqrt(3/2*(dsig_CC'*dsig_CC));
                end
                
                ElemS(stres) = sigmas_CC;
            end
        else
            disp('error, not yet implemented');
        end
        ElemS;
        
%% Case 50. Fine Scale
   case 50

       if flag_FF == 1
           
           ElemS = zeros(3+2,1);
           ElemE = zeros(3,1);
           ElemV = zeros(2,1);
           Nmat = zeros(2,ndf*nel);
           Bmat_FF = zeros(3,ndf*nel);
           
           % Load Gauss Points for quadrature
           if nel == 3
               lint = lintt3;%13;
           elseif nel == 4
               lint = lintq4;
           elseif nel == 6
               lint = lintt6;%13;
           elseif nel == 9
               lint = lintq9;
           end
           
           lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
           mu = PatchE/(2*(1+Patchv));
           
           cvec = zeros(3,1); %temperature vector, see Hughes book
           if PSPS == 'n'
               Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
           else
               Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                   Patchv  1      0
                   0      0      (1-Patchv)/2];
           end
           
           der = 0;
           bf = 0;
           ib = 0;
           
           ulres = reshape(ul,ndf*nen,1);
           
           for je = 1:lint
               
               if nel == 3 || nel == 6
                   [Wgt,litr,lits] = intpntt(je,lint,ib);
                   [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                   [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
               elseif nel == 4 || nel == 9
                   [Wgt,litr,lits] = intpntq(je,lint,ib);
                   [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                   [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
               end
               
               c1 = Wgt*Jdet*thick;
               
               % Form B matrix
               for ie = 1:nel
                   
                   Nmat(1,(ie-1)*ndf+1) = shl(ie);
                   Nmat(2,(ie-1)*ndf+2) = shl(ie);
                   
                   Bmat_FF(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
                   Bmat_FF(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
                   
               end
               
               xint = xl(:,1:nel)*shl;
               
               displac = Nmat*ulres(1:nel*ndf);
               epsil_FF = (Bmat_FF*ulres(1:nel*ndf) - cvec);
               sigma = Dmat*epsil;
               
               ElemE = ElemE + c1*epsil_FF;
               ElemS = ElemS + c1*[sigma; displac];
               ElemV = ElemV + c1*[1; sigma'*epsil_FF];
               
           end %je
       end
 %%   
    case 52
        
%         if flag_FF == 1
            ElemS = zeros(9,1);
            ElemE = zeros(3,1);
            ElemV = zeros(2,1);
            Nmat = zeros(2,ndf*nel);
            
            % Load Gauss Points for quadrature
            if nel == 3
                lint = lintt3;%13;
            elseif nel == 4
                lint = lintq4;
            elseif nel == 6
                lint = lintt6;%13;
            elseif nel == 9
                lint = lintq9;
            end
            
            lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
            mu = PatchE/(2*(1+Patchv));
            
            cvec = zeros(3,1); %temperature vector, see Hughes book
            if PSPS == 'n'
                Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
            else
                Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                    Patchv  1      0
                    0      0      (1-Patchv)/2];
            end
            Dvec = reshape(Dmat,9,1);
            
            der = 0;
            bf = 0;
            ib = 0;
            
            for je = 1:lint
                
                if nel == 3 || nel == 6
                    [Wgt,litr,lits] = intpntt(je,lint,ib);
                    [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                    [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
                elseif nel == 4 || nel == 9
                    [Wgt,litr,lits] = intpntq(je,lint,ib);
                    [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                    [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
                end
                
                c1 = Wgt*Jdet*thick;
                
                % Form B matrix
                for ie = 1:nel
                    
                    Nmat(1,(ie-1)*ndf+1) = shl(ie);
                    Nmat(2,(ie-1)*ndf+2) = shl(ie);
                    
                end
                
                xint = xl(:,1:nel)*shl;
                
                ElemS(1:9) = ElemS(1:9) + c1*Dvec;
                ElemE(1:2) = ElemE(1:2) + c1*xint;
                ElemV = ElemV + c1*[1; 0];
                
            end %je
%         end
end %Task Switch
