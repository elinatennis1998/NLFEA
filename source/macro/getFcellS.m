function eF = getFcellS(xc,uc,d,iprob,xl,ul,ndf,ndfs,ndm,lint,nel,nen,nummat,MR,BR,MS,BS,t11,t12,t21,t22)
%
%     No modifications when copied to NLFEA ver2
%**********************************************************************
% 
% 
%       Implicit None
% 
%       logical debug
%       common /debugs/ debug
% 
% %     Input Variables
%       integer ndf,ndfs,ndm,nel,nen,lint,nummat
%       real*8  MR,BR,MS,BS
%       real*8  xl(ndm,nen),xc(ndm,nen),uc(ndfs,nen)
% 
% %     Output Variables
%       real*8 eF(nen*ndf)
% 
% %     Local Variables
%       integer i,j,l
% 
%       real*8 w,det,x,y,
%      >       ux_x,ux_y,uy_x,uy_y,djx,djy,djn,two,zero,
%      >       shl(3,10),shls(3,10),shg(3,10),shgs(3,10),
%      >       shel(3,10),shels(3,10),she(3,10),shes(3,10),
%      >       bigR,bigS,litr,lits,det2
%       real*8 bubble(3),sx(2,2)
%       integer ib
%       logical der,bf,dprt
% 
%       real*8 mu

%       two = 2.d0
%       zero = 0.d0
% 
% 	call pzero(eF,nen*ndf)
% 
% 
% % Determine formulation type of cell
% 
%       if(ndfs.eq.2) then  % Fine scale type cell
% 
% %   Initialize parameters
% 
%          ib = 0
% 	     der = .false.
% 	     bf = .false.
% 
% %.....---------------------------------------------------
% %     Integration of matrix and vector terms over element
% %.....---------------------------------------------------
%          do 250 l=1,lint
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call intpntt(l,lint,ib,w,litr,lits)
%               call shlt(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgt(xc,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call intpntq(l,lint,ib,w,litr,lits)
%               call shlq(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgq(xc,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             endif
% 
% 
% 
% %.....---------------------------------------------------
% %     Evaluate partition of unity, uno
% %.....---------------------------------------------------
% % Evaluate R,S of element in terms of r,s of cell
%             bigR = MR*litr + BR
%             bigS = MS*lits + BS
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'litr',litr
%          write(24,*) 'lits',lits
%          write(24,*) 'bigR',bigR
%          write(24,*) 'bigS',bigS
%          write(24,'(A)') ''
%       endif
%       endif
% 
% % Compute element shape functions Nbar(R,S)
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call shlt(bigR,bigS,nel,der,bf,shel,shels,bubble)
%               call shgt(xl,nel,shel,shels,nummat,nen,bf,der,det2,she,
%      &                shes,bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call shlq(bigR,bigS,nel,der,bf,shel,shels,bubble)
%               call shgq(xl,nel,shel,shels,nummat,nen,bf,der,det2,she,
%      &                shes,bubble,sx)
%             endif
% 
% 	if(debug) then
% %	if(l.eq.lint) then
% %	 dprt = .true.
% %	else
% 	 dprt = .false.
% %	endif
% 	if(dprt) then
% 	   write(24,*) 'shg'
% 	   do j = 1,nel
%             write(24,1002) (shg(i,j),i=1,3)
% 	   enddo
%  1002 format (1x,d10.4,1x,d10.4,1x,d10.4)
% 	   write(24,'(A)') ''
%          write(24,*) 'cell values'
%          do j = 1,nel
%             write(24,*) (uc(i,j),i=1,ndfs)
%          enddo
%          write(24,'(A)') ''
% 	endif
% 	endif
% 
% %.....---------------------------------------------------
% %     Compute residual of coarse scale equations
% %.....---------------------------------------------------
% 
%             x = zero
%             y = zero
%             ux_x = zero
%             ux_y = zero
%             uy_x = zero
%             uy_y = zero
% 
%             do j = 1,nel
%                x = x + shg(3,j)*xc(1,j)
%                y = y + shg(3,j)*xc(2,j)
%                ux_x = ux_x + shg(1,j)*uc(1,j)
%                ux_y = ux_y + shg(2,j)*uc(1,j)
%                uy_x = uy_x + shg(1,j)*uc(2,j)
%                uy_y = uy_y + shg(2,j)*uc(2,j)
%             enddo
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'x',x
%          write(24,*) 'y',y
%          write(24,*) 'ux_x',ux_x
%          write(24,*) 'ux_y',ux_y
%          write(24,*) 'uy_x',uy_x
%          write(24,*) 'uy_y',uy_y
%          write(24,'(A)') ''
%       endif
%       endif
% 
% %.....---------------------------------------------------
% %     Loop over nodes to compute contributions to rows of 
% %     the stiffness matrix
% %.....---------------------------------------------------
% 
%             do 200 j=1,nel
%                djx=she(1,j)*det*w
%                djy=she(2,j)*det*w
%                djn=she(3,j)*det*w
% 
%                eF(ndf*j-2) = eF(ndf*j-2) - mu*(two*ux_x*djx 
% 	&                       + uy_x*djy + ux_y*djy)
%                eF(ndf*j-1) = eF(ndf*j-1) - mu*(two*uy_y*djy
% 	&                       + uy_x*djx + ux_y*djx)
% 	         eF(ndf*j)   = eF(ndf*j) - djn*(ux_x + uy_y) 
% 
%  200        continue %j rows
%  250     continue %l int points
% 
%       else %Stabilized cell type
% 
%          %   Initialize parameters
% 
%          ib = 0
% 	     der = .true.
% 	     bf = .true.
%          
% 
% %.....---------------------------------------------------
% %     Integration of matrix and vector terms over element
% %.....---------------------------------------------------
% 
% 
%       endif
% 
%       return
% 
%       end

%%
      two = 2.d0;
      zero = 0.d0;

	  eF = zeros(nen*ndf,1);

%       if(iprob==5)
% 
%          grav = 9.81d0;
%          fx = zero;
%          fy = -rho*grav;
% 
%       else

         fx = zero;
         fy = zero;

%       end

% Determine formulation type of cell

      if(ndfs==2)  % Fine scale type cell

%   Initialize parameters

         ib = 0;
	     der = 1;
	     bf = 1;
        I1 = [1; 1; 0; 0];
        spvec0 = I1;
        cpmat1 = I1*I1';
        cpmat2 = diag([-2,-2,-1,0]);
        cpmat0 = cpmat1 + cpmat2;

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
         for l=1:lint

            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            else
              [w,litr,lits] =  intpntq(l,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end


%.....---------------------------------------------------
%     Evaluate partition of unity, uno
%.....---------------------------------------------------

% Evaluate R,S of element in terms of r,s of cell
            bigR = MR*litr + BR;
            bigS = MS*lits + BS;
            
% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shlp,shed,shels,bubble] = shlt(bigR,bigS,nel,nel,der,bf);
              [she] = shgt(xl,nel,shed,shels,nen,bf,der,bubble);
            else
              [shlp,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she] = shgq(xl,nel,shed,shels,nen,bf,der,bubble);
            end

            [F,JxX,fi,Qxy] = kine2d(she,ul,nel,1);
            Qxyc = shg*fi;
            press = ul(3,:)*shlp;
            [theta,theta1,theta2,theta3] = ThetaS(JxX,d);
            spvec = JxX*theta1*spvec0;
            cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
            sigmap = press*spvec;
            cmatp = press*cpmat;
            [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,d);
            sigma2 = sigmai + sigmap;
            Dmat = cmati + cmatp;
            Smat = [sigma2(1) 0  sigma2(3)/2 sigma2(3)/2
                    0 sigma2(2)  sigma2(3)/2 -sigma2(3)/2
                    sigma2(3)/2  sigma2(3)/2 (sigma2(2)+sigma2(1))/4 (sigma2(2)-sigma2(1))/4
                    sigma2(3)/2 -sigma2(3)/2 (sigma2(2)-sigma2(1))/4 (sigma2(2)+sigma2(1))/4];
            Dmat = Dmat + Smat;
            Dmat = [Dmat
                    spvec']; %#ok<AGROW>
                
            if nel == 3            
            % Form B matrix    
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1)];
            elseif nel == 4            
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        0         Qxy(2,1)  0        0         Qxy(3,1)  0        0         Qxy(4,1)  0        0         Qxy(5,1)  0        0         Qxy(6,1)  0        0         Qxy(7,1)  0        0         Qxy(8,1)  0        0         Qxy(9,1)  0        0
                    0         Qxy(1,2) 0         0         Qxy(2,2) 0         0         Qxy(3,2) 0         0         Qxy(4,2) 0         0         Qxy(5,2) 0         0         Qxy(6,2) 0         0         Qxy(7,2) 0         0         Qxy(8,2) 0         0         Qxy(9,2) 0
                    Qxy(1,2)  Qxy(1,1) 0         Qxy(2,2)  Qxy(2,1) 0         Qxy(3,2)  Qxy(3,1) 0         Qxy(4,2)  Qxy(4,1) 0         Qxy(5,2)  Qxy(5,1) 0         Qxy(6,2)  Qxy(6,1) 0         Qxy(7,2)  Qxy(7,1) 0         Qxy(8,2)  Qxy(8,1) 0         Qxy(9,2)  Qxy(9,1) 0
                    Qxy(1,2) -Qxy(1,1) 0         Qxy(2,2) -Qxy(2,1) 0         Qxy(3,2) -Qxy(3,1) 0         Qxy(4,2) -Qxy(4,1) 0         Qxy(5,2) -Qxy(5,1) 0         Qxy(6,2) -Qxy(6,1) 0         Qxy(7,2) -Qxy(7,1) 0         Qxy(8,2) -Qxy(8,1) 0         Qxy(9,2) -Qxy(9,1) 0
                    0        0         shlp(1,1) 0        0         shlp(2,1) 0        0         shlp(3,1) 0        0         shlp(4,1) 0        0         shlp(5,1) 0        0         shlp(6,1) 0        0         shlp(7,1) 0        0         shlp(8,1) 0        0         shlp(9,1)];
            end

%.....---------------------------------------------------
%     Compute residual of coarse scale equations
%.....---------------------------------------------------

            ux_x = uc(1,:)*Qxyc(:,1);
            ux_y = uc(1,:)*Qxyc(:,2);
            uy_x = uc(2,:)*Qxyc(:,1);
            uy_y = uc(2,:)*Qxyc(:,2);
            BF = [ux_x+0
                  0+uy_y
                  ux_y+uy_x
                  ux_y-uy_x];

%.....---------------------------------------------------
%     Loop over nodes to compute contributions to rows of 
%     the stiffness matrix
%.....---------------------------------------------------

            eF = eF - det*w*Bmat'*Dmat*BF;
            
         end %l int points

      else %Stabilized cell type

%   Initialize parameters

         ib = 0;
	     der = 1;
	     bf = 1;
        
        [t11,t12,t21,t22] = Tau(xl,mu,nel,nen,lint);

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
         for l=1:lint

            if(nel==3||nel==6)
              [w,litr,lits] = intpntt(l,lint,ib);
              [shl,shls,bubble] =  shlt(litr,lits,nel,der,bf);
              [shg,shgs,det] = shgt(xc,nel,shl,shls,nen,bf,der,bubble);
            else
              [w,litr,lits] = intpntq(l,lint,ib);
              [shl,shls,bubble] =  shlq(litr,lits,nel,der,bf);
              [shg,shgs,det] = shgq(xc,nel,shl,shls,nen,bf,der,bubble);
            end
            b = bubble(3);

%.....---------------------------------------------------
%     Evaluate partition of unity, uno
%.....---------------------------------------------------
% Evaluate R,S of element in terms of r,s of cell
            bigR = MR*litr + BR;
            bigS = MS*lits + BS;

% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shel,shels,bubble] =  shlt(bigR,bigS,nel,der,bf);
              [she,shes] = shgt(xl,nel,shel,shels,nen,bf,der,bubble);
            else
              [shel,shels,bubble] =  shlq(bigR,bigS,nel,der,bf);
              [she,shes] = shgq(xl,nel,shel,shels,nen,bf,der,bubble);
            end
%             b = bubble(3);

%.....---------------------------------------------------
%     Compute residual of coarse scale equations
%.....---------------------------------------------------

            x = zero;
            y = zero;
            p = zero;
            px = zero;
            py = zero;
            ux_x = zero;
            ux_y = zero;
            uy_x = zero;
            uy_y = zero;
            ux_xx = zero;
            ux_yy = zero;
            ux_xy = zero;
            uy_xx = zero;
            uy_yy = zero;
            uy_xy = zero;

            for j = 1:nel
               x = x + shg(3,j)*xc(1,j);
               y = y + shg(3,j)*xc(2,j);
               p = p + shg(3,j)*uc(3,j);
               px = px + shg(1,j)*uc(3,j);
               py = py + shg(2,j)*uc(3,j);
               ux_x = ux_x + shg(1,j)*uc(1,j);
               ux_y = ux_y + shg(2,j)*uc(1,j);
               uy_x = uy_x + shg(1,j)*uc(2,j);
               uy_y = uy_y + shg(2,j)*uc(2,j);
               ux_xx = ux_xx + shgs(1,j)*uc(1,j);
               ux_yy = ux_yy + shgs(2,j)*uc(1,j);
               ux_xy = ux_xy + shgs(3,j)*uc(1,j);
               uy_xx = uy_xx + shgs(1,j)*uc(2,j);
               uy_yy = uy_yy + shgs(2,j)*uc(2,j);
               uy_xy = uy_xy + shgs(3,j)*uc(2,j);
            end
            
            rx = fx + px + mu*(two*ux_xx + ux_yy + uy_xy);
            ry = fy + py + mu*(uy_xx + ux_xy + two*uy_yy);
            rp = p/lam - (ux_x + uy_y);
            uno = 1; unox = 0; unoy = 0;

%.....---------------------------------------------------
%     Loop over nodes to compute contributions to rows of 
%     the stiffness matrix
%.....---------------------------------------------------

            for j=1:nel
               djx=she(1,j)*det*w;
               djy=she(2,j)*det*w;
               djn=she(3,j)*det*w;

               djxx=shes(1,j)*det*w;
               djyy=shes(2,j)*det*w;
               djxy=shes(3,j)*det*w;

               eF(ndf*j-2) = eF(ndf*j-2) + ...
                          uno*(mu*(djyy*rx*b*t11+djyy*ry*b*t12 +...
                          2*djxx*(rx*b*t11 + ry*b*t12) ...
                          + djxy*rx*b*t21 + djxy*ry*b*t22)) ... 
                          + uno*djn*fx - ((uno*djx+unox*djn)*p ...
                          + two*(uno*djx+unox*djn)*mu*ux_x ...
                          + (uno*djy+unoy*djn)*mu*(ux_y+uy_x));
               eF(ndf*j-1) = eF(ndf*j-1) ...
                          + uno*(mu*(djxy*(rx*b*t11 + ry*b*t12) ...
                          + (djxx + 2*djyy)*(rx*b*t21 + ry*b*t22))) ...
                          + uno*djn*fy - ((uno*djy+unoy*djn)*p ...
                          + two*(uno*djy+unoy*djn)*mu*uy_y ...
                          + (uno*djx+unox*djn)*mu*(ux_y+uy_x));
               eF(ndf*j)   = eF(ndf*j) ...
                          + uno*(djx*rx*b*t11 + djx*ry*b*t12 ...
                          + djy*rx*b*t21 + djy*ry*b*t22 ...
                          + djn*rp);
%                djx=she(1,j)*det*w;
%                djy=she(2,j)*det*w;
%                djn=she(3,j)*det*w;
% 
%                eF(ndf*j-2) = eF(ndf*j-2) - mu*(two*ux_x*djx + uy_x*djy + ux_y*djy);
%                eF(ndf*j-1) = eF(ndf*j-1) - mu*(two*uy_y*djy + uy_x*djx + ux_y*djx);
% 	         eF(ndf*j)   = eF(ndf*j) - djn*(ux_x + uy_y) ;

            end %j rows
         end %l int points


      end
      