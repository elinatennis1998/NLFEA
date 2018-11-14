function [eK,eF] = getKFcellL3(xl,d,xe,ue,ndf,ndm,lint,nel,nen,nummat,strong,MR,BR,MS,BS,sslot,T11,T12,T21,T22)%)%
%
%...  Written by Tim Truster (Fall 2009)
%...  Program to compute stiffness matrix and force vector for cell 
%     of submesh. Stiffness matrix eK is formed using fine scale 
%     equation for ndf = 2, and using stabilized formulation for
%     ndf = 3. If strong is true, then the strong form of the residual
%     is computed for the coarse scale quantities; otherwise, the weak
%     form of the residual of the coarse scales is computed.
%
%**********************************************************************
% 
%       Implicit None
% 
%       logical debug
%       common /debugs/ debug
% 
% %     Input Variables
%       integer ndf,ndm,lint,nel,nen,nummat,sslot
%       real*8  MR,BR,MS,BS
%       real*8  d(*),xl(ndm,nen),xe(ndm,nen),ue(ndm+1,nen)
%       logical strong
% 
% %     Output Variables
%       real*8 eK(nen*ndf,nen*ndf),eF(nen*ndf)
% 
% %     Local Variables
%       integer i,j,l
%       integer iprob,IQTY
% 
%       real*8 w,det,x,y,px,py,ux_xx,ux_yy,ux_xy,
%      >       uy_xx,uy_yy,uy_xy,fx,fy,rx,ry,
%      >       p,ux_x,ux_y,uy_x,uy_y,
%      >       djx,djy,djn,djxx,djyy,djxy,
%      >       dix,diy,din,dixx,diyy,dixy,       
%      >       t11,t12,t21,t22,f11,f12,f21,f22,b,
%      >       two,zero,uno,
%      >       shl(3,10),shls(3,10),shg(3,10),shgs(3,10),
%      >       shel(3,10),shels(3,10),she(3,10),shes(3,10),
%      >       bigR,bigS,litr,lits,det2
%       real*8 bubble(3),sx(2,2)
%       integer ib
%       logical der,bf,dprt
% 
%       integer i_cond_field,Lsym,xinc,yinc,minc
%       real*8 vis,a1,a2,eps,rho,lam,young,nu,grav

%       two = 2.d0
%       zero = 0.d0
% 
% 	
%       call pzero(eK,nen*ndf*nen*ndf)
% 	call pzero(eF,nen*ndf)
% 
% 
% % Get material parameters
%       Call M_Prop(d,vis,a1,a2,eps,rho,Iprob,i_cond_field,Lsym,
%      *	        lam,young,nu,IQTY,xinc,yinc,minc)
% 
%       if(iprob.eq.2) then
% 
%          grav = 9.81d0
%          fx = zero
%          fy = -rho*grav
% 
%       else
% 
%          fx = zero
%          fy = zero
% 
%       endif
% 
% % Determine formulation type of cell
% 
%       if(ndf.eq.2) then  % Fine scale type cell
% 
%          if(strong) then % Strong residual form
% 
% %   Initialize parameters
% 
%          ib = 0
% 	     der = .true.
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
%               call shgt(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call intpntq(l,lint,ib,w,litr,lits)
%               call shlq(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgq(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             endif
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
% % Evaluate uno(R,S)
%             if(sslot.eq.0) then
% 	         uno = 1.d0
% 	      else
%                if(nel.eq.3.or.nel.eq.6) then
%                  call shlt(bigR,bigS,3,der,bf,shel,shels,bubble)
%                elseif(nel.eq.4.or.nel.eq.9) then
%                  call shlq(bigR,bigS,4,der,bf,shel,shels,bubble)
%                endif
%                uno = shel(3,sslot)
% 	      endif
% 
% % Compute element shape functions Nbar(R,S)
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call shlt(bigR,bigS,nel,der,bf,shel,shels,bubble)
%               call shgt(xe,nel,shel,shels,nummat,nen,bf,der,det2,she,
%      &                shes,bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call shlq(bigR,bigS,nel,der,bf,shel,shels,bubble)
%               call shgq(xe,nel,shel,shels,nummat,nen,bf,der,det2,she,
%      &                shes,bubble,sx)
%             endif
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'uno',uno
%          write(24,'(A)') ''
%       endif
%       endif

%%
      two = 2.d0;
      zero = 0.d0;

	
    eK = zeros(nen*ndf,nen*ndf);
	eF = zeros(nen*ndf,1);


% Get material parameters
cellE = d(1);
cellv = d(2);
rho = d(4);
PSPS = d(5);
iprob = d(6);
if PSPS == 0 %plane stress
    cellE = cellE*(1 + 2*cellv)/(1 + cellv)^2;
    cellv = cellv/(1 + cellv);
    lam = cellv*cellE/((1+cellv)*(1-2*cellv));
    vis = cellE/(2*(1+cellv));
else %plane strain
    lam = cellv*cellE/((1+cellv)*(1-2*cellv));
    vis = cellE/(2*(1+cellv));
end

      if(iprob==5)

         grav = 9.81d0;
         fx = zero;
         fy = -rho*grav;

      else

         fx = zero;
         fy = zero;

      end

% Determine formulation type of cell

      if(ndf==2)  % Fine scale type cell

         if(strong) % Strong residual form

%   Initialize parameters

         ib = 0;
	     der = 1;
	     bf = 0;

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
         for l=1:lint

            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,der,bf);
              [shg, shgs, det] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            else
              [w,litr,lits] =  intpntq(l,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,der,bf);
              [shg, shgs, det] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end

%.....---------------------------------------------------
%     Evaluate partition of unity, uno
%.....---------------------------------------------------
% Evaluate R,S of element in terms of r,s of cell
            bigR = MR*litr + BR;
            bigS = MS*lits + BS;

% Evaluate uno(R,S)
            if(sslot==0)
	         uno = 1.d0;
	        else
               if(nel==3||nel==6)
                 [shel] = shlt(bigR,bigS,3,der,bf);
               else
                 [shel] = shlq(bigR,bigS,4,der,bf);
               end
               uno = shel(3,sslot);
            end

% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shel,shed,shels,bubble] = shlt(bigR,bigS,nel,der,bf);
              [she,shes] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,der,bf);
              [she,shes] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end

% %.....---------------------------------------------------
% %     Compute residual of coarse scale equations
% %.....---------------------------------------------------
% 
%             x = zero
%             y = zero
%             px = zero
%             py = zero
%             ux_xx = zero
%             ux_yy = zero
%             ux_xy = zero
%             uy_xx = zero
%             uy_yy = zero
%             uy_xy = zero
% 
%             do j = 1,nel
%                x = x + she(3,j)*xe(1,j)
%                y = y + she(3,j)*xe(2,j)
%                px = px + she(1,j)*ue(3,j)
%                py = py + she(2,j)*ue(3,j)
%                ux_xx = ux_xx + shes(1,j)*ue(1,j)
%                ux_yy = ux_yy + shes(2,j)*ue(1,j)
%                ux_xy = ux_xy + shes(3,j)*ue(1,j)
%                uy_xx = uy_xx + shes(1,j)*ue(2,j)
%                uy_yy = uy_yy + shes(2,j)*ue(2,j)
%                uy_xy = uy_xy + shes(3,j)*ue(2,j)
%             enddo
%                             
%             rx = fx + px + vis*(two*ux_xx + ux_yy + uy_xy)
%             ry = fy + py + vis*(uy_xx + ux_xy + two*uy_yy)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'x',x
%          write(24,*) 'y',y
%          write(24,*) 'px',px
%          write(24,*) 'py',py
%          write(24,*) 'ux_xx',ux_xx
%          write(24,*) 'ux_yy',ux_yy
%          write(24,*) 'ux_xy',ux_xy
%          write(24,*) 'uy_xx',uy_xx
%          write(24,*) 'uy_yy',uy_yy
%          write(24,*) 'uy_xy',uy_xy
%          write(24,*) 'rx',rx
%          write(24,*) 'ry',ry
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
%                djx=shg(1,j)*det*w
%                djy=shg(2,j)*det*w
%                djn=shg(3,j)*det*w
% 
%                eF(ndf*j-1) = eF(ndf*j-1) + djn*rx*uno
%                eF(ndf*j)   = eF(ndf*j)   + djn*ry*uno
% 
% %.....---------------------------------------------------
% %     Loop over nodes to compute contributions to columns
% %     of the stiffness matrix
% %.....---------------------------------------------------   
% 
%             do 210 i=1,nel
%                dix=shg(1,i)
%                diy=shg(2,i)
% 
% %* (1,1)*************************************************************************************
%                eK(ndf*i-1,ndf*j-1) = eK(ndf*i-1,ndf*j-1) 
%      &         + vis*(2.d0*dix*djx + diy*djy)
% 
% %* (1,2)*************************************************************************************
%                eK(ndf*i-1,ndf*j)   = eK(ndf*i-1,ndf*j) 
%      &         + vis*diy*djx
% 
% %* (2,1)*************************************************************************************
%                eK(ndf*i,ndf*j-1)   = eK(ndf*i,ndf*j-1)
%      &         + vis*dix*djy
%      				    
% %* (2,2)*************************************************************************************
%                eK(ndf*i,ndf*j)     = eK(ndf*i,ndf*j)
%      &         + vis*(dix*djx + 2.d0*diy*djy)
% 
%  210        continue %i columns
%  200        continue %j rows
%  250     continue %l int points

%%
%.....---------------------------------------------------
%     Compute residual of coarse scale equations
%.....---------------------------------------------------

            x = zero;
            y = zero;
            px = zero;
            py = zero;
            ux_xx = zero;
            ux_yy = zero;
            ux_xy = zero;
            uy_xx = zero;
            uy_yy = zero;
            uy_xy = zero;

            for j = 1:nel
               x = x + she(3,j)*xe(1,j);
               y = y + she(3,j)*xe(2,j);
               px = px + she(1,j)*ue(3,j);
               py = py + she(2,j)*ue(3,j);
               ux_xx = ux_xx + shes(1,j)*ue(1,j);
               ux_yy = ux_yy + shes(2,j)*ue(1,j);
               ux_xy = ux_xy + shes(3,j)*ue(1,j);
               uy_xx = uy_xx + shes(1,j)*ue(2,j);
               uy_yy = uy_yy + shes(2,j)*ue(2,j);
               uy_xy = uy_xy + shes(3,j)*ue(2,j);
            end
                            
            rx = fx + px + vis*(two*ux_xx + ux_yy + uy_xy);
            ry = fy + py + vis*(uy_xx + ux_xy + two*uy_yy);

%.....---------------------------------------------------
%     Loop over nodes to compute contributions to rows of 
%     the stiffness matrix
%.....---------------------------------------------------

            for j=1:nel
               djx=shg(1,j)*det*w;
               djy=shg(2,j)*det*w;
               djn=shg(3,j)*det*w;

               eF(ndf*j-1) = eF(ndf*j-1) + djn*rx*uno;
               eF(ndf*j)   = eF(ndf*j)   + djn*ry*uno;

%.....---------------------------------------------------
%     Loop over nodes to compute contributions to columns
%     of the stiffness matrix
%.....---------------------------------------------------   

            for i=1:nel
               dix=shg(1,i);
               diy=shg(2,i);

%* (1,1)*************************************************************************************
               eK(ndf*i-1,ndf*j-1) = eK(ndf*i-1,ndf*j-1) + vis*(2.d0*dix*djx + diy*djy);

%* (1,2)*************************************************************************************
               eK(ndf*i-1,ndf*j)   = eK(ndf*i-1,ndf*j) + vis*diy*djx;

%* (2,1)*************************************************************************************
               eK(ndf*i,ndf*j-1)   = eK(ndf*i,ndf*j-1) + vis*dix*djy;
     				    
%* (2,2)*************************************************************************************
               eK(ndf*i,ndf*j)     = eK(ndf*i,ndf*j) + vis*(dix*djx + 2.d0*diy*djy);
               
            end %i columns
            end %j rows
         end %l int points

%          else % Weak residual form
% 
%  %   Initialize parameters
% 
%          ib = 0
% 	     der = .false.
% 	     bf = .false.
% 
% %.....---------------------------------------------------
% %     Integration of matrix and vector terms over element
% %.....---------------------------------------------------
%          do 350 l=1,lint
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call intpntt(l,lint,ib,w,litr,lits)
%               call shlt(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgt(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call intpntq(l,lint,ib,w,litr,lits)
%               call shlq(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgq(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             endif
% 
% %.....---------------------------------------------------
% %     Evaluate partition of unity, uno
% %.....---------------------------------------------------
% % Evaluate R,S of element in terms of r,s of cell
%             bigR = MR*litr + BR
%             bigS = MS*lits + BS
% 
%       if(debug) then
%       dprt = .false.
%       if(dprt) then
%          write(24,*) 'litr',litr
%          write(24,*) 'lits',lits
%          write(24,*) 'bigR',bigR
%          write(24,*) 'bigS',bigS
%          write(24,'(A)') ''
%       endif
%       endif
% 
% % Evaluate uno(R,S)
%             if(sslot.eq.0) then
% 	         uno = 1.d0
% 	      else
%                if(nel.eq.3.or.nel.eq.6) then
%                  call shlt(bigR,bigS,3,der,bf,shel,shels,bubble)
%                elseif(nel.eq.4.or.nel.eq.9) then
%                  call shlq(bigR,bigS,4,der,bf,shel,shels,bubble)
%                endif
%                uno = shel(3,sslot)
% 	      endif
%             
% % Compute element shape functions Nbar(R,S)
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call shlt(bigR,bigS,nel,der,bf,shel,shels,bubble)
%               call shgt(xe,nel,shel,shels,nummat,nen,bf,der,det2,she,
%      &                shes,bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call shlq(bigR,bigS,nel,der,bf,shel,shels,bubble)
%               call shgq(xe,nel,shel,shels,nummat,nen,bf,der,det2,she,
%      &                shes,bubble,sx)
%             endif
% 
%       if(debug) then
%       dprt = .false.
%       if(dprt) then
%          write(24,*) 'uno',uno
%          write(24,'(A)') ''
%       endif
%       endif
%             
% %.....---------------------------------------------------
% %     Compute residual of coarse scale equations
% %.....---------------------------------------------------
% 
%             x = zero
%             y = zero
%             p = zero
%             ux_x = zero
%             ux_y = zero
%             uy_x = zero
%             uy_y = zero
% 
%             do j = 1,nel
%                x = x + she(3,j)*xe(1,j)
%                y = y + she(3,j)*xe(2,j)
%                p = p + she(3,j)*ue(3,j)
%                ux_x = ux_x + she(1,j)*ue(1,j)
%                ux_y = ux_y + she(2,j)*ue(1,j)
%                uy_x = uy_x + she(1,j)*ue(2,j)
%                uy_y = uy_y + she(2,j)*ue(2,j)
%             enddo
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'x',x
%          write(24,*) 'y',y
%          write(24,*) 'p',p
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
%             do 300 j=1,nel
%                djx=shg(1,j)*det*w
%                djy=shg(2,j)*det*w
%                djn=shg(3,j)*det*w
% 
%                eF(ndf*j-1) = eF(ndf*j-1) + uno*(djn*fx - 
%      &           (djx*p + two*djx*vis*ux_x + djy*vis*(ux_y+uy_x)))
%                eF(ndf*j)   = eF(ndf*j)   + uno*(djn*fy - 
%      &           (djy*p + two*djy*vis*uy_y + djx*vis*(ux_y+uy_x)))
% 
% %.....---------------------------------------------------
% %     Loop over nodes to compute contributions to columns
% %     of the stiffness matrix
% %.....---------------------------------------------------   
% 
%             do 310 i=1,nel
%                dix=shg(1,i)
%                diy=shg(2,i)
% 
% %* (1,1)*************************************************************************************
%                eK(ndf*i-1,ndf*j-1) = eK(ndf*i-1,ndf*j-1) 
%      &         + vis*(2.d0*dix*djx + diy*djy)
% 
% %* (1,2)*************************************************************************************
%                eK(ndf*i-1,ndf*j)   = eK(ndf*i-1,ndf*j) 
%      &         + vis*diy*djx
% 
% %* (2,1)*************************************************************************************
%                eK(ndf*i,ndf*j-1)   = eK(ndf*i,ndf*j-1)
%      &         + vis*dix*djy
%      				    
% %* (2,2)*************************************************************************************
%                eK(ndf*i,ndf*j)     = eK(ndf*i,ndf*j)
%      &         + vis*(dix*djx + 2.d0*diy*djy)
% 
%  310        continue %i columns
%  300        continue %j rows
%  350     continue %l int points
% 
%          endif % Residual form

%%
         else % Weak residual form

 %   Initialize parameters

         ib = 0;
	     der = 1;
	     bf = 1;

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

% Evaluate uno(R,S)
            if(sslot==0)
	         uno = 1.d0;
             unox = 0;
             unoy = 0;
	        else
               if(nel==3||nel==6)
                 [shel,shed] = shlt(bigR,bigS,3,3,der,bf);
                 shed = shgt(xe,3,shed,shls,nen,bf,der,be);
               else
                 [shel,shed] = shlq(bigR,bigS,4,4,der,bf);
                 shed = shgq(xe,4,shed,shls,nen,bf,der,be);
               end
               uno = shel(sslot);
               unox = shed(sslot,1);
               unoy = shed(sslot,2);
            end

% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shel,shed,shels,bubble] = shlt(bigR,bigS,nel,nel,der,bf);
              [she,shes,ted,bubble] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she,shes,ted,bubble] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end
            
%.....---------------------------------------------------
%     Compute residual of coarse scale equations
%.....---------------------------------------------------

%             x = zero;
%             y = zero;
%             p = zero;
%             ux_x = zero;
%             ux_y = zero;
%             uy_x = zero;
%             uy_y = zero;
%             px = zero;
%             py = zero;
%             ux_xx = zero;
%             ux_yy = zero;
%             ux_xy = zero;
%             uy_xx = zero;
%             uy_yy = zero;
%             uy_xy = zero;

            x = xe(1,:)*shel;
            y = xe(2,:)*shel;
            p = ue(3,:)*shel;
            ux_x = ue(1,:)*she(:,1);
            ux_y = ue(1,:)*she(:,2);
            uy_x = ue(2,:)*she(:,1);
            uy_y = ue(2,:)*she(:,2);
            px = ue(3,:)*she(:,1);
            py = ue(3,:)*she(:,2);
            ux_xx = ue(1,:)*shes(:,1);
            ux_yy = ue(1,:)*shes(:,2);
            ux_xy = ue(1,:)*shes(:,3);
            uy_xx = ue(2,:)*shes(:,1);
            uy_yy = ue(2,:)*shes(:,2);
            uy_xy = ue(2,:)*shes(:,3);
%             for j = 1:nel
%                x = x + she(j)*xe(1,j);
%                y = y + she(j)*xe(2,j);
%                p = p + she(j)*ue(3,j);
%                ux_x = ux_x + she(1,j)*ue(1,j);
%                ux_y = ux_y + she(2,j)*ue(1,j);
%                uy_x = uy_x + she(1,j)*ue(2,j);
%                uy_y = uy_y + she(2,j)*ue(2,j);
%                px = px + she(1,j)*ue(3,j);
%                py = py + she(2,j)*ue(3,j);
%                ux_xx = ux_xx + shes(1,j)*ue(1,j);
%                ux_yy = ux_yy + shes(2,j)*ue(1,j);
%                ux_xy = ux_xy + shes(3,j)*ue(1,j);
%                uy_xx = uy_xx + shes(1,j)*ue(2,j);
%                uy_yy = uy_yy + shes(2,j)*ue(2,j);
%                uy_xy = uy_xy + shes(3,j)*ue(2,j);
%             end
                            
            rx = fx + px + vis*(two*ux_xx + ux_yy + uy_xy);
            ry = fy + py + vis*(uy_xx + ux_xy + two*uy_yy);
            
            ux_x = (T11*rx+T12*ry)*bubble(1) + ux_x;
            ux_y = (T11*rx+T12*ry)*bubble(2) + ux_y;
            uy_x = (T21*rx+T22*ry)*bubble(1) + uy_x;
            uy_y = (T21*rx+T22*ry)*bubble(2) + uy_y;
            
%.....---------------------------------------------------
%     Loop over nodes to compute contributions to rows of 
%     the stiffness matrix
%.....---------------------------------------------------

            for j=1:nel
               djx=shg(j,1)*det*w;
               djy=shg(j,2)*det*w;
               djn=shl(j)*det*w;

%                eF(ndf*j-1) = eF(ndf*j-1) + uno*(djn*fx - (djx*p + two*djx*vis*ux_x + djy*vis*(ux_y+uy_x)));
%                eF(ndf*j)   = eF(ndf*j)   + uno*(djn*fy - (djy*p + two*djy*vis*uy_y + djx*vis*(ux_y+uy_x)));
               eF(ndf*j-1) = eF(ndf*j-1) + uno*djn*fx - ((uno*djx+unox*djn)*p + ...
                             two*(uno*djx+unox*djn)*vis*ux_x + ...
                             (uno*djy+unoy*djn)*vis*(ux_y+uy_x));
               eF(ndf*j)   = eF(ndf*j)   + uno*djn*fy - ((uno*djy+unoy*djn)*p + ...
                             two*(uno*djy+unoy*djn)*vis*uy_y + ...
                             (uno*djx+unox*djn)*vis*(ux_y+uy_x));

%.....---------------------------------------------------
%     Loop over nodes to compute contributions to columns
%     of the stiffness matrix
%.....---------------------------------------------------   

            for i=1:nel
               dix=shg(i,1);
               diy=shg(i,2);

%* (1,1)*************************************************************************************
               eK(ndf*i-1,ndf*j-1) = eK(ndf*i-1,ndf*j-1) + vis*(2.d0*dix*djx + diy*djy);

%* (1,2)*************************************************************************************
               eK(ndf*i-1,ndf*j)   = eK(ndf*i-1,ndf*j) + vis*diy*djx;

%* (2,1)*************************************************************************************
               eK(ndf*i,ndf*j-1)   = eK(ndf*i,ndf*j-1) + vis*dix*djy;
     				    
%* (2,2)*************************************************************************************
               eK(ndf*i,ndf*j)     = eK(ndf*i,ndf*j) + vis*(dix*djx + 2.d0*diy*djy);

            end %i columns
            end %j rows
         end %l int points

         end % Residual form

%       else %Stabilized cell type
% 
%          %   Initialize parameters
% 
%          ib = 0
% 	     der = .true.
% 	     bf = .true.
%          
%          %.....Stabilization Tensor Tau (Consistent Formulation) 5.25.08
% 
% 	     Call Tau(xl,vis,nel,ndm,t11,t12,t21,t22,f11,f12,f21,f22)
% 
% %.....---------------------------------------------------
% %     Integration of matrix and vector terms over element
% %.....---------------------------------------------------
% 	
%          do 450 l=1,lint
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call intpntt(l,lint,ib,w,litr,lits)
%               call shlt(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgt(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call intpntq(l,lint,ib,w,litr,lits)
%               call shlq(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgq(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             endif
% 	        b = bubble(3)
% 
% %.....---------------------------------------------------
% %     Loop over nodes to compute contributions to force
% %     vector and rows of the stiffness matrix
% %.....---------------------------------------------------
% 
%             do 400 j=1,nel
%                djx=shg(1,j)*det*w
%                djy=shg(2,j)*det*w
%                djn=shg(3,j)*det*w
% 
%                djxx=shgs(1,j)*det*w
%                djyy=shgs(2,j)*det*w
%                djxy=shgs(3,j)*det*w
%                
% %.....---------------------------------------------------
% %     Loop over nodes to compute contributions to columns
% %     of the stiffness matrix
% %.....---------------------------------------------------   
% 
%             do 410 i=1,nel
%                dix=shg(1,i)
%                diy=shg(2,i)
%                din=shg(3,i)
% 
%                dixx=shgs(1,i)
%                diyy=shgs(2,i)
%                dixy=shgs(3,i)
% 
% %* (1,1)*************************************************************************************
%             eK(ndf*i-2,ndf*j-2)=eK(ndf*i-2,ndf*j-2) 
% %Term6x:Diffusion-Diffusion
%      &    - vis*vis*(2.d0*dixx*(2.d0*djxx*b*t11+djyy*b*t11+djxy*b*t12) 
%      &    + diyy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) 
%      &    + dixy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) 
% %Term20x:Standard Diffusion
%      &    + vis*(2.d0*dix*djx + diy*djy)
% 
% %* (1,2)*************************************************************************************
%             eK(ndf*i-2,ndf*j-1)=eK(ndf*i-2,ndf*j-1)  
% %Term6x:Diffusion-Diffusion
%      &    - vis*vis*(diyy*(djxy*b*t11 + djxx*b*t12 + 2.d0*djyy*b*t12) 
%      &    + 2.d0*dixx*(djxy*b*t11 + (djxx+2.d0*djyy)*b*t12) 
%      &    + dixy*(djxy*b*t21 + djxx*b*t22 + 2.d0*djyy*b*t22)) 
% %Term20x:Standard Diffusion
%      &    + vis*diy*djx
% 
% %* (1,3)*************************************************************************************
%             eK(ndf*i-2,ndf*j  )=eK(ndf*i-2,ndf*j  ) 
% %Term7x:Diffusion-Pressure
%      &    - vis*(diyy*djx*b*t11 + diyy*djy*b*t12 + 2.d0*dixx*(djx*b*t11 
%      &    + djy*b*t12) + dixy*djx*b*t21 +  dixy*djy*b*t22) 	
% %Term21x:Standard Pressure
%      &    + dix*djn
%       
% %* (2,1)*************************************************************************************
%             eK(ndf*i-1,ndf*j-2)=eK(ndf*i-1,ndf*j-2)  
% %Term6y:Diffusion-Diffusion
%      &    - vis*vis*(dixy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) 
%      &    + (dixx + 2.d0*diyy)*(2.d0*djxx*b*t21+djyy*b*t21+djxy*b*t22))     
% %Term20y:Standard Diffusion
%      &    + vis*dix*djy
%      				    
% %* (2,2)*************************************************************************************
%             eK(ndf*i-1,ndf*j-1)=eK(ndf*i-1,ndf*j-1) 
% %Term6y:Diffusion-Diffusion
%      &    - vis*vis*(dixy*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) 
%      &    + (dixx + 2.d0*diyy)*(djxy*b*t21 + (djxx + 2.d0*djyy)*b*t22)) 
% %Term20y:Standard Diffusion
%      &    + vis*(dix*djx + 2.d0*diy*djy)
% 
% %* (2,3)*************************************************************************************
%             eK(ndf*i-1,ndf*j  )=eK(ndf*i-1,ndf*j  ) 
% %Term7y:Diffusion-Pressure
%      &    - vis*(dixy*(djx*b*t11 + djy*b*t12) 
%      &    + (dixx + 2.d0*diyy)*(djx*b*t21 + djy*b*t22)) 
% %Term21y:Standard Pressure
%      &    + diy*djn
%      			 
% %* (3,1)*************************************************************************************
%             eK(ndf*i  ,ndf*j-2)=eK(ndf*i  ,ndf*j-2) 
% %Term5:Continuity-Diffusion
%      &    - vis*(dix*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) 
%      &    + diy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) 
% %Term22x:Standard Continuity
%      &    + din*djx
%        
% %* (3,2)*************************************************************************************
%             eK(ndf*i  ,ndf*j-1)=eK(ndf*i  ,ndf*j-1) 
% %Term5:Continuity-Diffusion
%      &    - vis*(dix*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) 
%      &    + diy*(djxy*b*t21 + (djxx + 2.d0*djyy)*b*t22)) 
% %Term22y:Standard Continuity
%      &    + din*djy
% 
% %* (3,3)*************************************************************************************
%             eK(ndf*i  ,ndf*j  )=eK(ndf*i  ,ndf*j  ) 
% %Term1:Continuity-Pressure
%      &    -dix*djx*b*t11 - dix*djy*b*t12 - diy*djx*b*t21 - diy*djy*b*t22
% %Term2:Nearly Incompressible Elasticity-Pressure
%      &    -din*djn*(1.d0/lam)
% 
%  410        continue %i columns
%  400        continue %j rows
%  450     continue %l int points
% 
%       endif
% 
%       return
% 
%       end

      else %Stabilized cell type

         if(strong) % Strong residual form

         %   Initialize parameters

         ib = 0;
	     der = 1;
	     bf = 1;
         
         %.....Stabilization Tensor Tau (Consistent Formulation) 5.25.08

	     [t11,t12,t21,t22] = Tau3_2d(xl,vis,nel,nen,lint);

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
	
         for l=1:lint

            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,lint,ib);
              [shl,shld,shls,bubble] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, bubble] = shgt(xl,nel,shld,shls,nen,bf,der,bubble);
            else
              [w,litr,lits] =  intpntq(l,lint,ib);
              [shl,shld,shls,bubble] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, bubble] = shgq(xl,nel,shld,shls,nen,bf,der,bubble);
            end
	        b = bubble(3);
            
%.....---------------------------------------------------
%     Evaluate partition of unity, uno
%.....---------------------------------------------------
% Evaluate R,S of element in terms of r,s of cell
            bigR = MR*litr + BR;
            bigS = MS*lits + BS;

% Evaluate uno(R,S)
            if(sslot==0)
	         uno = 1.d0;
	        else
               if(nel==3||nel==6)
                 [shel] = shlt(bigR,bigS,3,3,der,bf);
               else
                 [shel] = shlq(bigR,bigS,4,4,der,bf);
               end
               uno = shel(sslot);
            end

% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shel,shed,shels,bubble] = shlt(bigR,bigS,nel,nel,der,bf);
              [she,shes] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she,shes] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end

% %.....---------------------------------------------------
% %     Compute residual of coarse scale equations
% %.....---------------------------------------------------
% 
%             x = zero
%             y = zero
%             px = zero
%             py = zero
%             ux_xx = zero
%             ux_yy = zero
%             ux_xy = zero
%             uy_xx = zero
%             uy_yy = zero
%             uy_xy = zero
% 
%             do j = 1,nel
%                x = x + she(3,j)*xe(1,j)
%                y = y + she(3,j)*xe(2,j)
%                px = px + she(1,j)*ue(3,j)
%                py = py + she(2,j)*ue(3,j)
%                ux_xx = ux_xx + shes(1,j)*ue(1,j)
%                ux_yy = ux_yy + shes(2,j)*ue(1,j)
%                ux_xy = ux_xy + shes(3,j)*ue(1,j)
%                uy_xx = uy_xx + shes(1,j)*ue(2,j)
%                uy_yy = uy_yy + shes(2,j)*ue(2,j)
%                uy_xy = uy_xy + shes(3,j)*ue(2,j)
%             enddo
%                             
%             rx = fx + px + vis*(two*ux_xx + ux_yy + uy_xy)
%             ry = fy + py + vis*(uy_xx + ux_xy + two*uy_yy)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'x',x
%          write(24,*) 'y',y
%          write(24,*) 'px',px
%          write(24,*) 'py',py
%          write(24,*) 'ux_xx',ux_xx
%          write(24,*) 'ux_yy',ux_yy
%          write(24,*) 'ux_xy',ux_xy
%          write(24,*) 'uy_xx',uy_xx
%          write(24,*) 'uy_yy',uy_yy
%          write(24,*) 'uy_xy',uy_xy
%          write(24,*) 'rx',rx
%          write(24,*) 'ry',ry
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
%                djx=shg(1,j)*det*w
%                djy=shg(2,j)*det*w
%                djn=shg(3,j)*det*w
% 
%                eF(ndf*j-1) = eF(ndf*j-1) + djn*rx*uno
%                eF(ndf*j)   = eF(ndf*j)   + djn*ry*uno
% 
% %.....---------------------------------------------------
% %     Loop over nodes to compute contributions to columns
% %     of the stiffness matrix
% %.....---------------------------------------------------   
% 
%             do 210 i=1,nel
%                dix=shg(1,i)
%                diy=shg(2,i)
% 
% %* (1,1)*************************************************************************************
%                eK(ndf*i-1,ndf*j-1) = eK(ndf*i-1,ndf*j-1) 
%      &         + vis*(2.d0*dix*djx + diy*djy)
% 
% %* (1,2)*************************************************************************************
%                eK(ndf*i-1,ndf*j)   = eK(ndf*i-1,ndf*j) 
%      &         + vis*diy*djx
% 
% %* (2,1)*************************************************************************************
%                eK(ndf*i,ndf*j-1)   = eK(ndf*i,ndf*j-1)
%      &         + vis*dix*djy
%      				    
% %* (2,2)*************************************************************************************
%                eK(ndf*i,ndf*j)     = eK(ndf*i,ndf*j)
%      &         + vis*(dix*djx + 2.d0*diy*djy)
% 
%  210        continue %i columns
%  200        continue %j rows
%  250     continue %l int points

%%
%.....---------------------------------------------------
%     Compute residual of coarse scale equations
%.....---------------------------------------------------

            x = zero;
            y = zero;
            p = zero;
            ux_x = zero;
            uy_y = zero;
            px = zero;
            py = zero;
            ux_xx = zero;
            ux_yy = zero;
            ux_xy = zero;
            uy_xx = zero;
            uy_yy = zero;
            uy_xy = zero;

            for j = 1:nel
               x = x + she(3,j)*xe(1,j);
               y = y + she(3,j)*xe(2,j);
               p = p + she(3,j)*ue(3,j);
               ux_x = ux_x + she(1,j)*ue(1,j);
               uy_y = uy_y + she(2,j)*ue(2,j);
               px = px + she(1,j)*ue(3,j);
               py = py + she(2,j)*ue(3,j);
               ux_xx = ux_xx + shes(1,j)*ue(1,j);
               ux_yy = ux_yy + shes(2,j)*ue(1,j);
               ux_xy = ux_xy + shes(3,j)*ue(1,j);
               uy_xx = uy_xx + shes(1,j)*ue(2,j);
               uy_yy = uy_yy + shes(2,j)*ue(2,j);
               uy_xy = uy_xy + shes(3,j)*ue(2,j);
            end
                            
            rx = fx + px + vis*(two*ux_xx + ux_yy + uy_xy);
            ry = fy + py + vis*(uy_xx + ux_xy + two*uy_yy);
            rp = p/lam - (ux_x + uy_y);

%.....---------------------------------------------------
%     Loop over nodes to compute contributions to force
%     vector and rows of the stiffness matrix
%.....---------------------------------------------------

            for j=1:nel
               djx=shg(1,j)*det*w;
               djy=shg(2,j)*det*w;
               djn=shg(3,j)*det*w;

               djxx=shgs(1,j)*det*w;
               djyy=shgs(2,j)*det*w;
               djxy=shgs(3,j)*det*w;

               eF(ndf*j-2) = eF(ndf*j-2) + ...
                          uno*(vis*(djyy*rx*b*t11+djyy*ry*b*t12 +...
                          2*djxx*(rx*b*t11 + ry*b*t12) ...
                          + djxy*rx*b*t21 + djxy*ry*b*t22) ... 
                          + djn*rx);
               eF(ndf*j-1) = eF(ndf*j-1) ...
                          + uno*(vis*(djxy*(rx*b*t11 + ry*b*t12) ...
                          + (djxx + 2*djyy)*(rx*b*t21 + ry*b*t22)) ...
                          + djn*ry);
               eF(ndf*j)   = eF(ndf*j) ...
                          + uno*(djx*rx*b*t11 + djx*ry*b*t12 ...
                          + djy*rx*b*t21 + djy*ry*b*t22 ...
                          + djn*rp);
               
%.....---------------------------------------------------
%     Loop over nodes to compute contributions to columns
%     of the stiffness matrix
%.....---------------------------------------------------   

            for i=1:nel
               dix=shg(1,i);
               diy=shg(2,i);
               din=shg(3,i);

               dixx=shgs(1,i);
               diyy=shgs(2,i);
               dixy=shgs(3,i);

%* (1,1)*************************************************************************************
            eK(ndf*i-2,ndf*j-2)=eK(ndf*i-2,ndf*j-2) ...
          - vis*vis*(2.d0*dixx*(2.d0*djxx*b*t11+djyy*b*t11+djxy*b*t12) ...
          + diyy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) ...
          + dixy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) ...
          + vis*(2.d0*dix*djx + diy*djy);

%* (1,2)*************************************************************************************
            eK(ndf*i-2,ndf*j-1)=eK(ndf*i-2,ndf*j-1) ...
          - vis*vis*(diyy*(djxy*b*t11 + djxx*b*t12 + 2.d0*djyy*b*t12) ...
          + 2.d0*dixx*(djxy*b*t11 + (djxx+2.d0*djyy)*b*t12) ...
          + dixy*(djxy*b*t21 + djxx*b*t22 + 2.d0*djyy*b*t22)) ...
          + vis*diy*djx;

%* (1,3)*************************************************************************************
            eK(ndf*i-2,ndf*j  )=eK(ndf*i-2,ndf*j  ) ...
          - vis*(diyy*djx*b*t11 + diyy*djy*b*t12 + 2.d0*dixx*(djx*b*t11 ...
          + djy*b*t12) + dixy*djx*b*t21 +  dixy*djy*b*t22) ...
          + dix*djn;
      
%* (2,1)*************************************************************************************
            eK(ndf*i-1,ndf*j-2)=eK(ndf*i-1,ndf*j-2) ...
          - vis*vis*(dixy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) ...
          + (dixx + 2.d0*diyy)*(2.d0*djxx*b*t21+djyy*b*t21+djxy*b*t22)) ...
          + vis*dix*djy;
     				    
%* (2,2)*************************************************************************************
            eK(ndf*i-1,ndf*j-1)=eK(ndf*i-1,ndf*j-1) ...
          - vis*vis*(dixy*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) ...
          + (dixx + 2.d0*diyy)*(djxy*b*t21 + (djxx + 2.d0*djyy)*b*t22)) ...
          + vis*(dix*djx + 2.d0*diy*djy);

%* (2,3)*************************************************************************************
            eK(ndf*i-1,ndf*j  )=eK(ndf*i-1,ndf*j  ) ...
          - vis*(dixy*(djx*b*t11 + djy*b*t12) ...
          + (dixx + 2.d0*diyy)*(djx*b*t21 + djy*b*t22)) ...
          + diy*djn;
     			 
%* (3,1)*************************************************************************************
            eK(ndf*i  ,ndf*j-2)=eK(ndf*i  ,ndf*j-2) ...
          - vis*(dix*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) ...
          + diy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) ...
          + din*djx;
       
%* (3,2)*************************************************************************************
            eK(ndf*i  ,ndf*j-1)=eK(ndf*i  ,ndf*j-1) ...
          - vis*(dix*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) ...
          + diy*(djxy*b*t21 + (djxx + 2.d0*djyy)*b*t22)) ...
          + din*djy;

%* (3,3)*************************************************************************************
            eK(ndf*i  ,ndf*j  )=eK(ndf*i  ,ndf*j  ) ...
          -dix*djx*b*t11 - dix*djy*b*t12 - diy*djx*b*t21 - diy*djy*b*t22 ...
          -din*djn*(1.d0/lam);

            end
            end
         end
         
         else % Weak residual form
             
             

         %   Initialize parameters

         ib = 0;
	     der = 1;
	     bf = 1;
         
         %.....Stabilization Tensor Tau (Consistent Formulation) 5.25.08

	     [t11,t12,t21,t22] = Tau3_2d(xl,vis,nel,nen,lint);

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
	
         for l=1:lint

            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,lint,ib);
              [shl,shld,shls,bubble] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, bubble] = shgt(xl,nel,shld,shls,nen,bf,der,bubble);
            else
              [w,litr,lits] =  intpntq(l,lint,ib);
              [shl,shld,shls,bubble] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, bubble] = shgq(xl,nel,shld,shls,nen,bf,der,bubble);
            end
	        b = bubble(3);
            
%.....---------------------------------------------------
%     Evaluate partition of unity, uno
%.....---------------------------------------------------
% Evaluate R,S of element in terms of r,s of cell
            bigR = MR*litr + BR;
            bigS = MS*lits + BS;

% Evaluate uno(R,S)
            if(sslot==0)
	         uno = 1.d0;
             unox = 0;
             unoy = 0;
	        else
               if(nel==3||nel==6)
                 [shel,shed,shels,be] = shlt(bigR,bigS,3,3,der,bf);
                 shed = shgt(xe,3,shed,shls,nen,bf,der,be);
               else
                 [shel,shed,shels,be] = shlq(bigR,bigS,4,4,der,bf);
                 shed = shgq(xe,4,shed,shls,nen,bf,der,be);
               end
               uno = shel(sslot);
               unox = shed(1,sslot);
               unoy = shed(2,sslot);
            end

% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shel,shed,shels,bubble] = shlt(bigR,bigS,nel,nel,der,bf);
              [she,shes,ted,bubble] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she,shes,ted,bubble] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end

% %.....---------------------------------------------------
% %     Compute residual of coarse scale equations
% %.....---------------------------------------------------
% 
%             x = zero
%             y = zero
%             px = zero
%             py = zero
%             ux_xx = zero
%             ux_yy = zero
%             ux_xy = zero
%             uy_xx = zero
%             uy_yy = zero
%             uy_xy = zero
% 
%             do j = 1,nel
%                x = x + she(3,j)*xe(1,j)
%                y = y + she(3,j)*xe(2,j)
%                px = px + she(1,j)*ue(3,j)
%                py = py + she(2,j)*ue(3,j)
%                ux_xx = ux_xx + shes(1,j)*ue(1,j)
%                ux_yy = ux_yy + shes(2,j)*ue(1,j)
%                ux_xy = ux_xy + shes(3,j)*ue(1,j)
%                uy_xx = uy_xx + shes(1,j)*ue(2,j)
%                uy_yy = uy_yy + shes(2,j)*ue(2,j)
%                uy_xy = uy_xy + shes(3,j)*ue(2,j)
%             enddo
%                             
%             rx = fx + px + vis*(two*ux_xx + ux_yy + uy_xy)
%             ry = fy + py + vis*(uy_xx + ux_xy + two*uy_yy)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'x',x
%          write(24,*) 'y',y
%          write(24,*) 'px',px
%          write(24,*) 'py',py
%          write(24,*) 'ux_xx',ux_xx
%          write(24,*) 'ux_yy',ux_yy
%          write(24,*) 'ux_xy',ux_xy
%          write(24,*) 'uy_xx',uy_xx
%          write(24,*) 'uy_yy',uy_yy
%          write(24,*) 'uy_xy',uy_xy
%          write(24,*) 'rx',rx
%          write(24,*) 'ry',ry
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
%                djx=shg(1,j)*det*w
%                djy=shg(2,j)*det*w
%                djn=shg(3,j)*det*w
% 
%                eF(ndf*j-1) = eF(ndf*j-1) + djn*rx*uno
%                eF(ndf*j)   = eF(ndf*j)   + djn*ry*uno
% 
% %.....---------------------------------------------------
% %     Loop over nodes to compute contributions to columns
% %     of the stiffness matrix
% %.....---------------------------------------------------   
% 
%             do 210 i=1,nel
%                dix=shg(1,i)
%                diy=shg(2,i)
% 
% %* (1,1)*************************************************************************************
%                eK(ndf*i-1,ndf*j-1) = eK(ndf*i-1,ndf*j-1) 
%      &         + vis*(2.d0*dix*djx + diy*djy)
% 
% %* (1,2)*************************************************************************************
%                eK(ndf*i-1,ndf*j)   = eK(ndf*i-1,ndf*j) 
%      &         + vis*diy*djx
% 
% %* (2,1)*************************************************************************************
%                eK(ndf*i,ndf*j-1)   = eK(ndf*i,ndf*j-1)
%      &         + vis*dix*djy
%      				    
% %* (2,2)*************************************************************************************
%                eK(ndf*i,ndf*j)     = eK(ndf*i,ndf*j)
%      &         + vis*(dix*djx + 2.d0*diy*djy)
% 
%  210        continue %i columns
%  200        continue %j rows
%  250     continue %l int points

%%
%.....---------------------------------------------------
%     Compute residual of coarse scale equations
%.....---------------------------------------------------

            x = zero;
            y = zero;
            p = zero;
            ux_x = zero;
            ux_y = zero;
            uy_x = zero;
            uy_y = zero;
            px = zero;
            py = zero;
            ux_xx = zero;
            ux_yy = zero;
            ux_xy = zero;
            uy_xx = zero;
            uy_yy = zero;
            uy_xy = zero;

            for j = 1:nel
               x = x + she(3,j)*xe(1,j);
               y = y + she(3,j)*xe(2,j);
               p = p + she(3,j)*ue(3,j);
               ux_x = ux_x + she(1,j)*ue(1,j);
               ux_y = ux_y + she(2,j)*ue(1,j);
               uy_x = uy_x + she(1,j)*ue(2,j);
               uy_y = uy_y + she(2,j)*ue(2,j);
               px = px + she(1,j)*ue(3,j);
               py = py + she(2,j)*ue(3,j);
               ux_xx = ux_xx + shes(1,j)*ue(1,j);
               ux_yy = ux_yy + shes(2,j)*ue(1,j);
               ux_xy = ux_xy + shes(3,j)*ue(1,j);
               uy_xx = uy_xx + shes(1,j)*ue(2,j);
               uy_yy = uy_yy + shes(2,j)*ue(2,j);
               uy_xy = uy_xy + shes(3,j)*ue(2,j);
            end
                            
            rx = fx + px + vis*(two*ux_xx + ux_yy + uy_xy);
            ry = fy + py + vis*(uy_xx + ux_xy + two*uy_yy);
            
            ux_x = (T11*rx+T12*ry)*bubble(1) + ux_x;
            ux_y = (T11*rx+T12*ry)*bubble(2) + ux_y;
            uy_x = (T21*rx+T22*ry)*bubble(1) + uy_x;
            uy_y = (T21*rx+T22*ry)*bubble(2) + uy_y;
            ux_xx = (T11*rx+T12*ry)*bubble(4);
            ux_xy = (T11*rx+T12*ry)*bubble(6);
            ux_yy = (T11*rx+T12*ry)*bubble(5);
            uy_xx = (T21*rx+T22*ry)*bubble(4);
            uy_xy = (T21*rx+T22*ry)*bubble(6);
            uy_yy = (T21*rx+T22*ry)*bubble(5);
            rx = rx + vis*(two*ux_xx + ux_yy + uy_xy);
            ry = ry + vis*(uy_xx + ux_xy + two*uy_yy);

            rp = p/lam - (ux_x + uy_y);

%.....---------------------------------------------------
%     Loop over nodes to compute contributions to force
%     vector and rows of the stiffness matrix
%.....---------------------------------------------------

            for j=1:nel
               djx=shg(1,j)*det*w;
               djy=shg(2,j)*det*w;
               djn=shg(3,j)*det*w;

               djxx=shgs(1,j)*det*w;
               djyy=shgs(2,j)*det*w;
               djxy=shgs(3,j)*det*w;

               eF(ndf*j-2) = eF(ndf*j-2) + ...
                          uno*(vis*(djyy*rx*b*t11+djyy*ry*b*t12 +...
                          2*djxx*(rx*b*t11 + ry*b*t12) ...
                          + djxy*rx*b*t21 + djxy*ry*b*t22)) ... 
                          + uno*djn*fx - ((uno*djx+unox*djn)*p ...
                          + two*(uno*djx+unox*djn)*vis*ux_x ...
                          + (uno*djy+unoy*djn)*vis*(ux_y+uy_x));
               eF(ndf*j-1) = eF(ndf*j-1) ...
                          + uno*(vis*(djxy*(rx*b*t11 + ry*b*t12) ...
                          + (djxx + 2*djyy)*(rx*b*t21 + ry*b*t22))) ...
                          + uno*djn*fy - ((uno*djy+unoy*djn)*p ...
                          + two*(uno*djy+unoy*djn)*vis*uy_y ...
                          + (uno*djx+unox*djn)*vis*(ux_y+uy_x));
               eF(ndf*j)   = eF(ndf*j) ...
                          + uno*(djx*rx*b*t11 + djx*ry*b*t12 ...
                          + djy*rx*b*t21 + djy*ry*b*t22 ...
                          + djn*rp);
               
%.....---------------------------------------------------
%     Loop over nodes to compute contributions to columns
%     of the stiffness matrix
%.....---------------------------------------------------   

            for i=1:nel
               dix=shg(1,i);
               diy=shg(2,i);
               din=shg(3,i);

               dixx=shgs(1,i);
               diyy=shgs(2,i);
               dixy=shgs(3,i);

%* (1,1)*************************************************************************************
            eK(ndf*i-2,ndf*j-2)=eK(ndf*i-2,ndf*j-2) ...
          - vis*vis*(2.d0*dixx*(2.d0*djxx*b*t11+djyy*b*t11+djxy*b*t12) ...
          + diyy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) ...
          + dixy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) ...
          + vis*(2.d0*dix*djx + diy*djy);

%* (1,2)*************************************************************************************
            eK(ndf*i-2,ndf*j-1)=eK(ndf*i-2,ndf*j-1) ...
          - vis*vis*(diyy*(djxy*b*t11 + djxx*b*t12 + 2.d0*djyy*b*t12) ...
          + 2.d0*dixx*(djxy*b*t11 + (djxx+2.d0*djyy)*b*t12) ...
          + dixy*(djxy*b*t21 + djxx*b*t22 + 2.d0*djyy*b*t22)) ...
          + vis*diy*djx;

%* (1,3)*************************************************************************************
            eK(ndf*i-2,ndf*j  )=eK(ndf*i-2,ndf*j  ) ...
          - vis*(diyy*djx*b*t11 + diyy*djy*b*t12 + 2.d0*dixx*(djx*b*t11 ...
          + djy*b*t12) + dixy*djx*b*t21 +  dixy*djy*b*t22) ...
          + dix*djn;
      
%* (2,1)*************************************************************************************
            eK(ndf*i-1,ndf*j-2)=eK(ndf*i-1,ndf*j-2) ...
          - vis*vis*(dixy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) ...
          + (dixx + 2.d0*diyy)*(2.d0*djxx*b*t21+djyy*b*t21+djxy*b*t22)) ...
          + vis*dix*djy;
     				    
%* (2,2)*************************************************************************************
            eK(ndf*i-1,ndf*j-1)=eK(ndf*i-1,ndf*j-1) ...
          - vis*vis*(dixy*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) ...
          + (dixx + 2.d0*diyy)*(djxy*b*t21 + (djxx + 2.d0*djyy)*b*t22)) ...
          + vis*(dix*djx + 2.d0*diy*djy);

%* (2,3)*************************************************************************************
            eK(ndf*i-1,ndf*j  )=eK(ndf*i-1,ndf*j  ) ...
          - vis*(dixy*(djx*b*t11 + djy*b*t12) ...
          + (dixx + 2.d0*diyy)*(djx*b*t21 + djy*b*t22)) ...
          + diy*djn;
     			 
%* (3,1)*************************************************************************************
            eK(ndf*i  ,ndf*j-2)=eK(ndf*i  ,ndf*j-2) ...
          - vis*(dix*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) ...
          + diy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) ...
          + din*djx;
       
%* (3,2)*************************************************************************************
            eK(ndf*i  ,ndf*j-1)=eK(ndf*i  ,ndf*j-1) ...
          - vis*(dix*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) ...
          + diy*(djxy*b*t21 + (djxx + 2.d0*djyy)*b*t22)) ...
          + din*djy;

%* (3,3)*************************************************************************************
            eK(ndf*i  ,ndf*j  )=eK(ndf*i  ,ndf*j  ) ...
          -dix*djx*b*t11 - dix*djy*b*t12 - diy*djx*b*t21 - diy*djy*b*t22 ...
          -din*djn*(1.d0/lam);

            end
            end
         end

         end % Residual form

      end
