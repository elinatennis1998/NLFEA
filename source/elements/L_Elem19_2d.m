% Tim Truster
% 10/15/2013
%
% Vectorize Poisson equation in 2D
% Allows input of full material tensor
% For use in residual-free bubble calculations with acoustic tensor from
% finite deformations.

switch isw %Task Switch
    
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
%%
    case 3

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        Bmat = zeros(4,ndf*nel);

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
        
        thick = 1;
 
        Dmat = reshape(mateprop,4,4); % all 16 material moduli components
        % so that W_i,I*A_iIjJ*dU_j,J = Bmat'*Dmat*Bmat with the ordering
        % as below; use the file [A_iIjJ,C_IJKL] = CSFtoA(C,S,F,ndm)

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
                    
%  Bmat(:,3*ie-2:3*ie) = [QXY(ie,1) 0         0         
%                          0         QXY(ie,1) 0         
%                          0         0         QXY(ie,1) 
%                          QXY(ie,2) 0         0         
%                          0         QXY(ie,2) 0
%                          0         0         QXY(ie,2) 
%                          QXY(ie,3) 0         0         
%                          0         QXY(ie,3) 0
%                          0         0         QXY(ie,3) ];
 Bmat(:,2*ie-1:2*ie) = [shg(ie,1) 0         
                         0         shg(ie,1) 
                         shg(ie,2) 0         
                         0         shg(ie,2) ];
                 
            end
            
            sigma = Dmat*(Bmat*ulres(1:nel*ndf));
            
            ElemF = ElemF - c1*Bmat'*sigma;
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            
        end %je
ElemK;
%%
    case -1 % boundary tractions
        
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
        thick = 1;
%         % Load Gauss Points for quadrature
%         if enrich == 1
%             [rlist, rWgts, rnum] = GaussPoints(pr+2);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         else
%             [rlist, rWgts, rnum] = GaussPoints(pr+1);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         end

%         lamda = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
%         mu = PatchE/(2*(1+Patchv));
%         thick = 1;
        ideriv = 0;
        der = 0;
        bf = 0;

        Nmat = zeros(2,ndf*nel);
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            else %if nel == 4
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

            if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            Traction = zeros(1,2);
            
            if iprob == 3
                
                z = xl(2,1:nel)*shl;
                T1 = h*q/12/kT*(6*z/h - 1) + 1/(2*h*rhoT*cT)*qtau;
                T2 = q/(4*kT*h)*z^2;
                T = T1 + T2;
                NT = alph*PatchE/(rhoT*cT)*qtau;
                MT = alph*PatchE*q/(3*kT)*h^3;
                A = 2*h;
                I = 2*h^3/3;
                s_xx = (-alph*PatchE*T + NT/A + MT/I*z);
                s_xy = 0;
                s_yy = 0;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif iprob == 4
                Traction = zeros(2,1);
            elseif iprob == 5
                
                x = xl(1,1:nel)*shl;
                y = xl(2,1:nel)*shl;
                grav = 9.81d0;
                Cwid = D;
                Len = L/2;
                x2 = x-Len;
                s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len^2-15.d0*x2^2 +4.d0*Cwid^2 +10.d0*y^2)/Cwid^2;
                s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid^2-y^2)/Cwid^2;
                s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid^2-y^2)/Cwid^2;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            elseif iprob == 6
                x = 0;
                y = 0;

                for j = 1:nel
                   x = x + xl(1,j)*shl(j);
                   y = y + xl(2,j)*shl(j);
                end
                Cwid = D;
                Len = L;
%                 Load = 2560;
                s_xx = 3.d0*(x-Len)*y*Load/(2*Cwid^3);
                s_xy = 3.d0/4.d0*(Cwid^2-y^2)*Load/Cwid^3;
                s_yy = 0;
                Traction(1) = s_xx*tu3(1) + s_xy*tu3(2);
                Traction(2) = s_xy*tu3(1) + s_yy*tu3(2);
            else
                Traction = traction(1:2);
            end
            
            c1 = Wgt*tm3*drdr*thick;
            
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
                 
            end
            
            ElemF = ElemF + c1*Nmat'*Traction';

        end %ie
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);

end