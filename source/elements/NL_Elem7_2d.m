% function ElemK = NL_Elem2_2d(ElemXYZ, nel, ndf, R, S, T, pr, ps, PatchE, Patchv)
%
% Tim Truster
% 7/2009
% UIUC
%
%   nel:             = number of nodes on current element
%
%   ndf:             = max number of DOF per node
%
%   ndm:             = space dimension of mesh
%
%   nst:             = size of element arrays (normally ndf*nel)
%
%
%   wf:              = weight factor; adds additional dimension to
%                      coordinate arrays (0 = Lagrangian elements, 
%                      1 = NURBS elements)
%
%   ElemXYZ:          = Local array containing (x,y) coordinates of nodes
%                      forming the element passed to the ComElemStif
%                      function
%
%                       X       Y
%                       1.54    -2.367
%
%   pr,ps,pt=0:      = degree of polynomial for elements in
%                      uvw-directions
%
%   mr,ms,mt:        = number of konts on element in uvw-directions
%
%   nr,ns,nt:        = number of nodes on element in uvw-directions
%
%   R,S,T:           = Element knot vectors
%

% CEE576 Fall 2012, Energy Conserving version
% Uses material model from the midterm that year

% Set Material Properties

aE = mateprop(2);
bE = mateprop(3);
rho0 = mateprop(4);

switch isw %Task Switch

    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(nst);
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);
        Nmat = zeros(2,nst);
        Bmat = zeros(3,nst);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(ll,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute rho, accel
                rho = rho0; %rho = f(x,y);
                accel = al*shl;
                
                for i = 1:nel
                    Bmat(:,2*i-1:2*i) = [Qxy(i,1) 0
                                         0 Qxy(i,2)
                                         Qxy(i,2) Qxy(i,1)];
                end
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                E(3) = E(3)/2; %for extra factor of 2 on shear strain
                Ekk = E(1) + E(2);
                lEkk = log(1+Ekk);
                sigma = aE/2*(lEkk+Ekk/(1+Ekk))*[1; 1; 0] + 3*bE*E;
                Ae = aE/2*(2+Ekk)/(1+Ekk)^2;
                D = [Ae+3*bE Ae 0 
                     Ae Ae+3*bE 0
                     0 0  3*bE/2];
                
                c1 = Wgt*Jdet*thick;
                c2 = Wgt*Jdet*thick*rho;
                %Evaluate Stiffness terms
                for o=1:nel
                    dox = Qxy(o,1)*c1;
                    doy = Qxy(o,2)*c1;
                    doN = shl(o)*c2;
                    BTsigma = [dox*sigma(1) + doy*sigma(3); dox*sigma(3) + doy*sigma(2)];

                    ElemF(ndf*o-1) = ElemF(ndf*o-1) - doN*accel(1) - coeff1*BTsigma(1);
                    ElemF(ndf*o) =   ElemF(ndf*o)   - doN*accel(2) - coeff1*BTsigma(2);
                    
                    ElemFn(ndf*o-1) = ElemFn(ndf*o-1) - Nalpha*BTsigma(1);
                    ElemFn(ndf*o) =   ElemFn(ndf*o)   - Nalpha*BTsigma(2);

                end %o
                ElemK = ElemK + c1*Bmat'*D*Bmat;

        end %je
        
    case 5 %Compute Mass
        
        ElemM = zeros(nst);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(ll,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                x = xl(1,:)*shl;
                y = xl(2,:)*shl;
%                 rho = rho0*(1-1/30*x); %rho = f(x,y);
                rho = rho0;
                
                c1 = Wgt*Jdet*thick;
                for o=1:nel
                    don = shl(o)*c1;

%                     for l=1:nel
                    for l=o:nel
                        dln = shl(l);
                        Mmat = don*rho*dln;
                        %(wx,ux)
                        ElemM(ndf*o-1,ndf*l-1) = ElemM(ndf*o-1,ndf*l-1) + ...
                                                 Mmat;
                        %(wy,uy)
                        ElemM(ndf*o,ndf*l)     = ElemM(ndf*o,ndf*l) + ...
                                                 Mmat;
                    end %l
                end %o

        end %je
        
        %Copy symmetric values
        
        for o=1:nel

            for l=o+1:nel
                
                Mmat = ElemM(ndf*o-1,ndf*l-1);
                
                %(wx,ux)
                ElemM(ndf*l-1,ndf*o-1) = Mmat;
                %(wy,uy)
                ElemM(ndf*l,ndf*o)     = Mmat;
                
            end %l
        end %o
% ElemM
    case 6 %Compute Residual
        
        ElemF = zeros(nst,1);
        ElemFn = zeros(nst,1);
        Nmat = zeros(2,nst);
        Bmat = zeros(3,nst);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(ll,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute rho, accel
                rho = rho0; %rho = f(x,y);
                accel = al*shl;
                
                for i = 1:nel
                    Bmat(:,2*i-1:2*i) = [Qxy(i,1) 0
                                         0 Qxy(i,2)
                                         Qxy(i,2) Qxy(i,1)];
                end
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                E(3) = E(3)/2; %for extra factor of 2 on shear strain
                Ekk = E(1) + E(2);
                lEkk = log(1+Ekk);
                sigma = aE/2*(lEkk+Ekk/(1+Ekk))*[1; 1; 0] + 3*bE*E;
                
                c1 = Wgt*Jdet*thick;
                c2 = Wgt*Jdet*thick*rho;
                %Evaluate Stiffness terms
                for o=1:nel
                    dox = Qxy(o,1)*c1;
                    doy = Qxy(o,2)*c1;
                    doN = shl(o)*c2;
                    BTsigma = [dox*sigma(1) + doy*sigma(3); dox*sigma(3) + doy*sigma(2)];

                    ElemF(ndf*o-1) = ElemF(ndf*o-1) - doN*accel(1) - coeff1*BTsigma(1);
                    ElemF(ndf*o) =   ElemF(ndf*o)   - doN*accel(2) - coeff1*BTsigma(2);
                    
                    ElemFn(ndf*o-1) = ElemFn(ndf*o-1) - Nalpha*BTsigma(1);
                    ElemFn(ndf*o) =   ElemFn(ndf*o)   - Nalpha*BTsigma(2);
                end %o

        end %je

    case 21 %Compute Stiffness
        
        ElemK = zeros(nst);
        Nmat = zeros(2,nst);
        Bmat = zeros(3,nst);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(ll,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute rho, accel
                rho = rho0; %rho = f(x,y);
                accel = al*shl;
                
                for i = 1:nel
                    Bmat(:,2*i-1:2*i) = [Qxy(i,1) 0
                                         0 Qxy(i,2)
                                         Qxy(i,2) Qxy(i,1)];
                end
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                E(3) = E(3)/2; %for extra factor of 2 on shear strain
                Ekk = E(1) + E(2);
                Ae = aE/2*(2+Ekk)/(1+Ekk)^2;
                D = [Ae+3*bE Ae 0 
                     Ae Ae+3*bE 0
                     0 0  3*bE/2];
                
                c1 = Wgt*Jdet*thick;
                
                ElemK = ElemK + c1*Bmat'*D*Bmat;

        end %je
        
    case 22
        
        ElemP = zeros(2*(3*ndf-3),nel);
        Bmat = zeros(3,nst);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(ll,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                %Compute strain and stress
                for i = 1:nel
                    Bmat(:,2*i-1:2*i) = [Qxy(i,1) 0
                                         0 Qxy(i,2)
                                         Qxy(i,2) Qxy(i,1)];
                end
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                E(3) = E(3)/2; %for extra factor of 2 on shear strain
                Ekk = E(1) + E(2);
                lEkk = log(1+Ekk);
                sigma = aE/2*(lEkk+Ekk/(1+Ekk))*[1; 1; 0] + 3*bE*E;
                
                ElemP(6*(ll-1)+1) = E(1);
                ElemP(6*(ll-1)+2) = E(2);
                ElemP(6*(ll-1)+3) = E(3);
                ElemP(6*(ll-1)+4) = sigma(1);
                ElemP(6*(ll-1)+5) = sigma(2);
                ElemP(6*(ll-1)+6) = sigma(3);

        end %je
        
    case 12 % energy
        
        ElemE = 0;
        Nmat = zeros(2,nst);
        Bmat = zeros(3,nst);
        
        thick = 1;
        
        % Load Guass Integration Points

        if nel == 3 || nel == 6
            lint = 13;
        elseif nel == 4
            lint = 4;
            % lint = 16;
        else
            lint = 9;
        end
        der = 0;
        bf = 0;
        ib = 0;

        %Integration Loop
        for ll = 1:lint

                %Evaluate first derivatives of basis functions at int. point
                if nel == 3 || nel == 6
                  [Wgt,litr,lits] =  intpntt(ll,lint,ib);
                  [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgt(xl,nel,shld,shls,nel,bf,der,be);
                else
                  [Wgt,litr,lits] =  intpntq(ll,lint,ib);
                  [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                  [Qxy, shgs, Jdet] = shgq(xl,nel,shld,shls,nel,bf,der,be);
                end
                
                for i = 1:nel
                    Bmat(:,2*i-1:2*i) = [Qxy(i,1) 0
                                         0 Qxy(i,2)
                                         Qxy(i,2) Qxy(i,1)];
                    Nmat(:,2*i-1:2*i) = [shl(i) 0
                                         0 shl(i)];
                end
                %Compute F, sigma, D
                E = Bmat*reshape(ul,nst,1);
                E(3) = E(3)/2; %for extra factor of 2 on shear strain
                Ekk = E(1) + E(2);
                lEkk = log(1+Ekk);
                Ue = aE/2*lEkk*Ekk + 3/2*bE*(E(1)^2 + E(2)^2 + 2*E(3)^2);
                
                c1 = Wgt*Jdet*thick;
                ElemE = ElemE + Ue*c1;
                
                rho = rho0;
                velo = vl*shl;
                De = 1/2*rho*(velo'*velo);
                ElemE = ElemE + De*c1;

        end %je

end