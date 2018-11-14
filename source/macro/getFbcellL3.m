function eF = getFbcellL3(xl,d,xe,ue,ndf,ndm,nel,nen,nummat,strong,MR,BR,MS,BS,sslot,tractyn)
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
%       integer ndf,ndm,nel,nen,nummat,sslot
%       real*8  MR,BR,MS,BS
%       real*8  d(*),xl(ndm,nen),xe(ndm,nen),ue(ndm+1,nen)
%       logical strong,tractyn
% 
% %     Output Variables
%       real*8 eF(nen*ndf)
% 
% %     Local Variables
%       integer i,j,l
%       integer iprob,IQTY
% 
%       real*8 w,det,x,y,hx,hy,
%      >       p,ux_x,ux_y,uy_x,uy_y,rx,ry,x2,
%      >       djx,djy,djn,t1(3),t2(3),t3(3),nv(3),
%      >       two,zero,uno,Cwid,Len,load,s_xx,s_xy,s_yy,
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
% 	hx = 0.d0
% 	hy = 0.d0
% 
% 	call pzero(eF,nen*ndf)
% 
% 
%       if(debug) then
% 	dprt = .false.
% 	if(dprt) then
%          write(24,*) 'sslot',sslot
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
% % Get material parameters
%       Call M_Prop(d,vis,a1,a2,eps,rho,Iprob,i_cond_field,Lsym,
%      *	        lam,young,nu,IQTY,xinc,yinc,minc)
% 
% % Determine formulation type of cell
% 
%       if(ndf.eq.2) then  % Fine scale type cell
% 
%          if(strong) then % Strong residual form
% 
% %   Initialize parameters
% 
%          ib = 1
% 	     der = .false.
% 	     bf = .false.
% 
% %.....---------------------------------------------------
% %     Integration of matrix and vector terms over element
% %.....---------------------------------------------------
%          do 250 l=1,4
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call intpntt(l,4,ib,w,litr,lits)
%               call shlt(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgt(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call intpntq(l,4,ib,w,litr,lits)
%               call shlq(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgq(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             endif
% 
% 	call pzero(t1,3)
% 	call pzero(t2,3)
% 	call pzero(t3,3)
% 
%       t1(1) = sx(1,1)
% 	t1(2) = sx(2,1)
%       
%       t2(3) = 1.d0
% 
%       t3(1) = t1(2)*t2(3) - t1(3)*t2(2)
%       t3(2) = t1(3)*t2(1) - t1(1)*t2(3)
%       t3(3) = t1(1)*t2(2) - t1(2)*t2(1)
%       call VecNormalize(t3,det,nv)
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
%          write(24,*) 'w',w
% 	   write(24,*) 'nel',nel
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 	if(debug) then
% %	if(l.eq.lint) then
% %	 dprt = .true.
% %	else
% 	 dprt = .false.
% %	endif
% 	if(dprt) then
% 	   write(24,*) 'shl'
% 	   do j = 1,nel
%             write(24,1002) (shl(i,j),i=1,3)
% 	   enddo
%        write(24,*) 'shg'
% 	   do j = 1,nel
%             write(24,1002) (shg(i,j),i=1,3)
% 	   enddo
%  1002 format (1x,d10.4,1x,d10.4,1x,d10.4)
% 	   write(24,'(A)') ''
% 	endif
% 	endif
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
% 	if(debug) then
% %	if(l.eq.lint) then
% %	 dprt = .true.
% %	else
% 	 dprt = .false.
% %	endif
% 	if(dprt) then
% 	   write(24,*) 'shel'
% 	   do j = 1,nel
%             write(24,1002) (shel(i,j),i=1,3)
% 	   enddo
%        write(24,*) 'she'
% 	   do j = 1,nel
%             write(24,1002) (she(i,j),i=1,3)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'uno',uno
% 	   write(24,*) 'det',det
% 	   write(24,*) 't1(1)',t1(1)
% 	   write(24,*) 't1(2)',t1(2)
% 	   write(24,*) 't3(1)',t3(1)
% 	   write(24,*) 't3(2)',t3(2)
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
% % Compute traction
%       if(tractyn) then
% 
%       if(iprob.eq.1) then
% 	   Cwid = 1.d0
%          Len = 10.d0
%          load = 2560.d0
%          s_xx = 3.d0*(x-Len)*y*load/(2*Cwid**2)
%          s_xy = 3.d0/4.d0*(Cwid**2-y**2)*load/Cwid**2
%          s_yy = -1.d0/4.d0*load*y*(Cwid**2-y**2)/Cwid**2
%       elseif(iprob.eq.2) then
%          grav = 9.81d0
%          Cwid = 1.d0
%          Len = 5.d0
%          x2 = x-Len
%          s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len**2
%      &              -15.d0*x2**2 +4.d0*Cwid**2 +10.d0*y**2)/Cwid**2
%          s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid**2-y**2)/Cwid**2
%          s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid**2
%      &              -y**2)/Cwid**2;
%       endif
% 
% 	   hx = s_xx*nv(1) + s_xy*nv(2)
% 	   hy = s_xy*nv(1) + s_yy*nv(2)
% 
% 	endif
% 
%             rx = hx-(p*nv(1) + vis*(uy_x*nv(2) 
%      &           + ux_y*nv(2) + two*ux_x*nv(1)))
%             ry = hy-(p*nv(2) + vis*(uy_x*nv(1) 
%      &           + ux_y*nv(1) + two*uy_y*nv(2)))
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
%          write(24,*) 's_xx',s_xx
%          write(24,*) 's_xy',s_xy
%          write(24,*) 's_yy',s_yy
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
%  200        continue %j rows
%  250     continue %l int points

%%
      two = 2.d0;
      zero = 0.d0;
	hx = 0.d0;
	hy = 0.d0;

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
%     lam = cellv*cellE/((1+cellv)*(1-2*cellv));
    vis = cellE/(2*(1+cellv));
else %plane strain
%     lam = cellv*cellE/((1+cellv)*(1-2*cellv));
    vis = cellE/(2*(1+cellv));
end

% Determine formulation type of cell

      if(ndf==2)  % Fine scale type cell

         if(strong) % Strong residual form

%   Initialize parameters

         ib = 1;
	     der = 0;
	     bf = 0;

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
         for l=1:4

            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,4,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            else
              [w,litr,lits] =  intpntq(l,4,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end

	t1 = zeros(1,3);
	t2 = t1;
	t3 = t1;

      t1(1) = sx(1,1);
	t1(2) = sx(2,1);
      
      t2(3) = 1.d0;

      t3(1) = t1(2)*t2(3) - t1(3)*t2(2);
      t3(2) = t1(3)*t2(1) - t1(1)*t2(3);
      t3(3) = t1(1)*t2(2) - t1(2)*t2(1);
      [det,nv] = VecNormalize(t3);

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
               uno = shel(3,sslot);
            end

% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shel,shed,shels,bubble] = shlt(bigR,bigS,nel,nel,der,bf);
              [she] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end

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

            for j = 1:nel
               x = x + she(3,j)*xe(1,j);
               y = y + she(3,j)*xe(2,j);
               p = p + she(3,j)*ue(3,j);
               ux_x = ux_x + she(1,j)*ue(1,j);
               ux_y = ux_y + she(2,j)*ue(1,j);
               uy_x = uy_x + she(1,j)*ue(2,j);
               uy_y = uy_y + she(2,j)*ue(2,j);
            end

% Compute traction
      if(tractyn)

      if(iprob==4)
	   Cwid = 1.d0;
         Len = 10.d0;
         load = 2560.d0;
         s_xx = 3.d0*(x-Len)*y*load/(2*Cwid^2);
         s_xy = 3.d0/4.d0*(Cwid^2-y^2)*load/Cwid^2;
         s_yy = -1.d0/4.d0*load*y*(Cwid^2-y^2)/Cwid^2;
      elseif(iprob==5)
         grav = 9.81d0;
         Cwid = 1.d0;
         Len = 5.d0;
         x2 = x-Len;
         s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len^2 ...
                    -15.d0*x2^2 +4.d0*Cwid^2 +10.d0*y^2)/Cwid^2;
         s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid^2-y^2)/Cwid^2;
         s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid^2 ...
                    -y^2)/Cwid^2;
      elseif(iprob==2)
          s_xx = 72*y;
          s_xy = 0;
          s_yy = 0;
      elseif(iprob==7)
          lam = 0.54448373678246398;
          Q = 0.54307557883673652;
	      r = sqrt(x^2+y^2);
	     if(y>=0.d0)
           theta = acos(x/r);
         else
	       theta = -acos(x/r);
         end
	     r = lam*r^(lam-1.d0);
         s_xx = r*((2.d0 - Q*(lam+1.d0))*cos((lam-1.d0)*theta) ...
              - (lam-1.d0)*cos((lam-3.d0)*theta));
         s_yy = r*((2.d0 + Q*(lam+1.d0))*cos((lam-1.d0)*theta) ...
              + (lam-1.d0)*cos((lam-3.d0)*theta));
         s_xy = r*((lam-1.d0)*sin((lam-3.d0)*theta) ...
              + Q*(lam+1.d0)*sin((lam-1.d0)*theta));
      end

	   hx = s_xx*nv(1) + s_xy*nv(2);
	   hy = s_xy*nv(1) + s_yy*nv(2);

      end

            rx = hx-(p*nv(1) + vis*(uy_x*nv(2) + ux_y*nv(2) + two*ux_x*nv(1)));
            ry = hy-(p*nv(2) + vis*(uy_x*nv(1) + ux_y*nv(1) + two*uy_y*nv(2)));

%.....---------------------------------------------------
%     Loop over nodes to compute contributions to rows of 
%     the stiffness matrix
%.....---------------------------------------------------

            for j=1:nel
               djn=shg(3,j)*det*w;

               eF(ndf*j-1) = eF(ndf*j-1) + djn*rx*uno;
               eF(ndf*j)   = eF(ndf*j)   + djn*ry*uno;

            end %j rows
         end %l int points

%          else % Weak residual form
% 	   if(tractyn) then %integrate only if on traction boundary
% 
%  %   Initialize parameters
% 
%          ib = 1
% 	     der = .false.
% 	     bf = .false.
% 
% %.....---------------------------------------------------
% %     Integration of matrix and vector terms over element
% %.....---------------------------------------------------
%          do 350 l=1,4
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call intpntt(l,4,ib,w,litr,lits)
%               call shlt(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgt(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call intpntq(l,4,ib,w,litr,lits)
%               call shlq(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgq(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             endif
% 	
% 	call pzero(t1,3)
% 	call pzero(t2,3)
% 	call pzero(t3,3)
% 
%       t1(1) = sx(1,1)
% 	t1(2) = sx(2,1)
%       
%       t2(3) = 1.d0
% 
%       t3(1) = t1(2)*t2(3) - t1(3)*t2(2)
%       t3(2) = t1(3)*t2(1) - t1(1)*t2(3)
%       t3(3) = t1(1)*t2(2) - t1(2)*t2(1)
%       call VecNormalize(t3,det,nv)
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
%          write(24,*) 'w',w
% 	   write(24,*) 'nel',nel
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 	if(debug) then
% %	if(l.eq.lint) then
% %	 dprt = .true.
% %	else
% 	 dprt = .false.
% %	endif
% 	if(dprt) then
% 	   write(24,*) 'shl'
% 	   do j = 1,nel
%             write(24,1002) (shl(i,j),i=1,3)
% 	   enddo
%        write(24,*) 'shg'
% 	   do j = 1,nel
%             write(24,1002) (shg(i,j),i=1,3)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
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
% 	if(debug) then
% %	if(l.eq.lint) then
% %	 dprt = .true.
% %	else
% 	 dprt = .false.
% %	endif
% 	if(dprt) then
% 	   write(24,*) 'shel'
% 	   do j = 1,nel
%             write(24,1002) (shel(i,j),i=1,3)
% 	   enddo
%        write(24,*) 'she'
% 	   do j = 1,nel
%             write(24,1002) (she(i,j),i=1,3)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 	
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'uno',uno
% 	   write(24,*) 'det',det
% 	   write(24,*) 't1(1)',t1(1)
% 	   write(24,*) 't1(2)',t1(2)
% 	   write(24,*) 't3(1)',t3(1)
% 	   write(24,*) 't3(2)',t3(2)
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
% 
%             do j = 1,nel
%                x = x + she(3,j)*xe(1,j)
%                y = y + she(3,j)*xe(2,j)
%             enddo
% 
%       if(iprob.eq.1) then
% 	   Cwid = 1.d0
%          Len = 10.d0
%          load = 2560.d0
%          s_xx = 3.d0*(x-Len)*y*load/(2*Cwid**2)
%          s_xy = 3.d0/4.d0*(Cwid**2-y**2)*load/Cwid**2
%          s_yy = -1.d0/4.d0*load*y*(Cwid**2-y**2)/Cwid**2
%       elseif(iprob.eq.2) then
%          grav = 9.81d0
%          Cwid = 1.d0
%          Len = 5.d0
%          x2 = x-Len
%          s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len**2
%      &              -15.d0*x2**2 +4.d0*Cwid**2 +10.d0*y**2)/Cwid**2
%          s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid**2-y**2)/Cwid**2
%          s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid**2
%      &              -y**2)/Cwid**2;
%       endif
% 
% 	   rx = s_xx*nv(1) + s_xy*nv(2)
% 	   ry = s_xy*nv(1) + s_yy*nv(2)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'x',x
%          write(24,*) 'y',y
%          write(24,*) 's_xx',s_xx
%          write(24,*) 's_xy',s_xy
%          write(24,*) 's_yy',s_yy
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
%             do 300 j=1,nel
%                djx=shg(1,j)*det*w
%                djy=shg(2,j)*det*w
%                djn=shg(3,j)*det*w
% 
%                eF(ndf*j-1) = eF(ndf*j-1) + djn*rx*uno
%                eF(ndf*j)   = eF(ndf*j)   + djn*ry*uno
% 
%  300        continue %j rows
%  350     continue %l int points
% 
%          endif
%          endif % Residual form
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
% 

%%
         else % Weak residual form
         if(tractyn) %integrate only if on traction boundary

 %   Initialize parameters

         ib = 1;
	     der = 0;
	     bf = 0;

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
         for l=1:4

            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,4,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            else
              [w,litr,lits] =  intpntq(l,4,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
	
    t1 = zeros(1,3);
	t2 = t1;
	t3 = t1;

      t1(1) = sx(1,1);
	t1(2) = sx(2,1);
      
      t2(3) = 1.d0;

      t3(1) = t1(2)*t2(3) - t1(3)*t2(2);
      t3(2) = t1(3)*t2(1) - t1(1)*t2(3);
      t3(3) = t1(1)*t2(2) - t1(2)*t2(1);
      [det,nv] = VecNormalize(t3);

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
              [she] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end
            
%.....---------------------------------------------------
%     Compute residual of coarse scale equations
%.....---------------------------------------------------

%             x = zero;
%             y = zero;

%             for j = 1:nel
%                x = x + she(3,j)*xe(1,j);
%                y = y + she(3,j)*xe(2,j);
%             end
            x = xe(1,:)*shel;
            y = xe(2,:)*shel;

      if(iprob==4)
	   Cwid = 1.d0;
         Len = 10.d0;
         load = 2560.d0;
         s_xx = 3.d0*(x-Len)*y*load/(2*Cwid^2);
         s_xy = 3.d0/4.d0*(Cwid^2-y^2)*load/Cwid^2;
         s_yy = -1.d0/4.d0*load*y*(Cwid^2-y^2)/Cwid^2;
      elseif(iprob==5)
         grav = 9.81d0;
         Cwid = 1.d0;
         Len = 5.d0;
         x2 = x-Len;
         s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len^2 ...
                    -15.d0*x2^2 +4.d0*Cwid^2 +10.d0*y^2)/Cwid^2;
         s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid^2-y^2)/Cwid^2;
         s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid^2 ...
                    -y^2)/Cwid^2;
      elseif(iprob==2)
          s_xx = 72*y;
          s_xy = 0;
          s_yy = 0;
      elseif(iprob==7)
          lam = 0.54448373678246398;
          Q = 0.54307557883673652;
	      r = sqrt(x^2+y^2);
	     if(y>=0.d0)
           theta = acos(x/r);
         else
	       theta = -acos(x/r);
         end
	     r = lam*r^(lam-1.d0);
         s_xx = r*((2.d0 - Q*(lam+1.d0))*cos((lam-1.d0)*theta) ...
              - (lam-1.d0)*cos((lam-3.d0)*theta));
         s_yy = r*((2.d0 + Q*(lam+1.d0))*cos((lam-1.d0)*theta) ...
              + (lam-1.d0)*cos((lam-3.d0)*theta));
         s_xy = r*((lam-1.d0)*sin((lam-3.d0)*theta) ...
              + Q*(lam+1.d0)*sin((lam-1.d0)*theta));
      elseif(iprob==8)
          s_xx = 0;
          s_xy = 6.25;
          s_yy = 0;
      end

	   rx = s_xx*nv(1) + s_xy*nv(2);
	   ry = s_xy*nv(1) + s_yy*nv(2);
            
%.....---------------------------------------------------
%     Loop over nodes to compute contributions to rows of 
%     the stiffness matrix
%.....---------------------------------------------------

            for j=1:nel
               djn=shl(j)*det*w;

               eF(ndf*j-1) = eF(ndf*j-1) + djn*rx*uno;
               eF(ndf*j)   = eF(ndf*j)   + djn*ry*uno;

            end %j rows
         end %l int points

         end
         end % Residual form

      else %Stabilized cell type

         if(strong) % Strong residual form

%   Initialize parameters

         ib = 1;
	     der = 0;
	     bf = 0;

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
         for l=1:4

            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,4,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            else
              [w,litr,lits] =  intpntq(l,4,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end

	t1 = zeros(1,3);
	t2 = t1;
	t3 = t1;

      t1(1) = sx(1,1);
	t1(2) = sx(2,1);
      
      t2(3) = 1.d0;

      t3(1) = t1(2)*t2(3) - t1(3)*t2(2);
      t3(2) = t1(3)*t2(1) - t1(1)*t2(3);
      t3(3) = t1(1)*t2(2) - t1(2)*t2(1);
      [det,nv] = VecNormalize(t3);

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
               uno = shel(3,sslot);
            end

% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shel,shed,shels,bubble] = shlt(bigR,bigS,nel,nel,der,bf);
              [she] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end

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

            for j = 1:nel
               x = x + she(3,j)*xe(1,j);
               y = y + she(3,j)*xe(2,j);
               p = p + she(3,j)*ue(3,j);
               ux_x = ux_x + she(1,j)*ue(1,j);
               ux_y = ux_y + she(2,j)*ue(1,j);
               uy_x = uy_x + she(1,j)*ue(2,j);
               uy_y = uy_y + she(2,j)*ue(2,j);
            end

% Compute traction
      if(tractyn)

      if(iprob==4)
	   Cwid = 1.d0;
         Len = 10.d0;
         load = 2560.d0;
         s_xx = 3.d0*(x-Len)*y*load/(2*Cwid^2);
         s_xy = 3.d0/4.d0*(Cwid^2-y^2)*load/Cwid^2;
         s_yy = -1.d0/4.d0*load*y*(Cwid^2-y^2)/Cwid^2;
      elseif(iprob==5)
         grav = 9.81d0;
         Cwid = 1.d0;
         Len = 5.d0;
         x2 = x-Len;
         s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len^2 ...
                    -15.d0*x2^2 +4.d0*Cwid^2 +10.d0*y^2)/Cwid^2;
         s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid^2-y^2)/Cwid^2;
         s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid^2 ...
                    -y^2)/Cwid^2;
      elseif(iprob==2)
          s_xx = 72*y;
          s_xy = 0;
          s_yy = 0;
      end

	   hx = s_xx*nv(1) + s_xy*nv(2);
	   hy = s_xy*nv(1) + s_yy*nv(2);

      end

            rx = hx-(p*nv(1) + vis*(uy_x*nv(2) + ux_y*nv(2) + two*ux_x*nv(1)));
            ry = hy-(p*nv(2) + vis*(uy_x*nv(1) + ux_y*nv(1) + two*uy_y*nv(2)));

%.....---------------------------------------------------
%     Loop over nodes to compute contributions to rows of 
%     the stiffness matrix
%.....---------------------------------------------------

            for j=1:nel
               djn=shg(3,j)*det*w;

               eF(ndf*j-2) = eF(ndf*j-2) + djn*rx*uno;
               eF(ndf*j-1) = eF(ndf*j-1) + djn*ry*uno;

            end %j rows
         end %l int points

%          else % Weak residual form
% 	   if(tractyn) then %integrate only if on traction boundary
% 
%  %   Initialize parameters
% 
%          ib = 1
% 	     der = .false.
% 	     bf = .false.
% 
% %.....---------------------------------------------------
% %     Integration of matrix and vector terms over element
% %.....---------------------------------------------------
%          do 350 l=1,4
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call intpntt(l,4,ib,w,litr,lits)
%               call shlt(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgt(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call intpntq(l,4,ib,w,litr,lits)
%               call shlq(litr,lits,nel,der,bf,shl,shls,bubble)
%               call shgq(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             endif
% 	
% 	call pzero(t1,3)
% 	call pzero(t2,3)
% 	call pzero(t3,3)
% 
%       t1(1) = sx(1,1)
% 	t1(2) = sx(2,1)
%       
%       t2(3) = 1.d0
% 
%       t3(1) = t1(2)*t2(3) - t1(3)*t2(2)
%       t3(2) = t1(3)*t2(1) - t1(1)*t2(3)
%       t3(3) = t1(1)*t2(2) - t1(2)*t2(1)
%       call VecNormalize(t3,det,nv)
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
%          write(24,*) 'w',w
% 	   write(24,*) 'nel',nel
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 	if(debug) then
% %	if(l.eq.lint) then
% %	 dprt = .true.
% %	else
% 	 dprt = .false.
% %	endif
% 	if(dprt) then
% 	   write(24,*) 'shl'
% 	   do j = 1,nel
%             write(24,1002) (shl(i,j),i=1,3)
% 	   enddo
%        write(24,*) 'shg'
% 	   do j = 1,nel
%             write(24,1002) (shg(i,j),i=1,3)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
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
% 	if(debug) then
% %	if(l.eq.lint) then
% %	 dprt = .true.
% %	else
% 	 dprt = .false.
% %	endif
% 	if(dprt) then
% 	   write(24,*) 'shel'
% 	   do j = 1,nel
%             write(24,1002) (shel(i,j),i=1,3)
% 	   enddo
%        write(24,*) 'she'
% 	   do j = 1,nel
%             write(24,1002) (she(i,j),i=1,3)
% 	   enddo
% 	   write(24,'(A)') ''
% 	endif
% 	endif
% 	
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'uno',uno
% 	   write(24,*) 'det',det
% 	   write(24,*) 't1(1)',t1(1)
% 	   write(24,*) 't1(2)',t1(2)
% 	   write(24,*) 't3(1)',t3(1)
% 	   write(24,*) 't3(2)',t3(2)
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
% 
%             do j = 1,nel
%                x = x + she(3,j)*xe(1,j)
%                y = y + she(3,j)*xe(2,j)
%             enddo
% 
%       if(iprob.eq.1) then
% 	   Cwid = 1.d0
%          Len = 10.d0
%          load = 2560.d0
%          s_xx = 3.d0*(x-Len)*y*load/(2*Cwid**2)
%          s_xy = 3.d0/4.d0*(Cwid**2-y**2)*load/Cwid**2
%          s_yy = -1.d0/4.d0*load*y*(Cwid**2-y**2)/Cwid**2
%       elseif(iprob.eq.2) then
%          grav = 9.81d0
%          Cwid = 1.d0
%          Len = 5.d0
%          x2 = x-Len
%          s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len**2
%      &              -15.d0*x2**2 +4.d0*Cwid**2 +10.d0*y**2)/Cwid**2
%          s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid**2-y**2)/Cwid**2
%          s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid**2
%      &              -y**2)/Cwid**2;
%       endif
% 
% 	   rx = s_xx*nv(1) + s_xy*nv(2)
% 	   ry = s_xy*nv(1) + s_yy*nv(2)
% 
%       if(debug) then
%       dprt = .true.
%       if(dprt) then
%          write(24,*) 'x',x
%          write(24,*) 'y',y
%          write(24,*) 's_xx',s_xx
%          write(24,*) 's_xy',s_xy
%          write(24,*) 's_yy',s_yy
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
%             do 300 j=1,nel
%                djx=shg(1,j)*det*w
%                djy=shg(2,j)*det*w
%                djn=shg(3,j)*det*w
% 
%                eF(ndf*j-1) = eF(ndf*j-1) + djn*rx*uno
%                eF(ndf*j)   = eF(ndf*j)   + djn*ry*uno
% 
%  300        continue %j rows
%  350     continue %l int points
% 
%          endif
%          endif % Residual form
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
% 

%%
         else % Weak residual form
         if(tractyn) %integrate only if on traction boundary

 %   Initialize parameters

         ib = 1;
	     der = 0;
	     bf = 0;

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
         for l=1:4

            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,4,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            else
              [w,litr,lits] =  intpntq(l,4,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
	
    t1 = zeros(1,3);
	t2 = t1;
	t3 = t1;

      t1(1) = sx(1,1);
	t1(2) = sx(2,1);
      
      t2(3) = 1.d0;

      t3(1) = t1(2)*t2(3) - t1(3)*t2(2);
      t3(2) = t1(3)*t2(1) - t1(1)*t2(3);
      t3(3) = t1(1)*t2(2) - t1(2)*t2(1);
      [det,nv] = VecNormalize(t3);

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
               uno = shel(3,sslot);
            end

% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shel,shed,shels,bubble] = shlt(bigR,bigS,nel,nel,der,bf);
              [she] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end
            
%.....---------------------------------------------------
%     Compute residual of coarse scale equations
%.....---------------------------------------------------

            x = zero;
            y = zero;

            for j = 1:nel
               x = x + she(3,j)*xe(1,j);
               y = y + she(3,j)*xe(2,j);
            end

      if(iprob==4)
	   Cwid = 1.d0;
         Len = 10.d0;
         load = 2560.d0;
         s_xx = 3.d0*(x-Len)*y*load/(2*Cwid^2);
         s_xy = 3.d0/4.d0*(Cwid^2-y^2)*load/Cwid^2;
         s_yy = -1.d0/4.d0*load*y*(Cwid^2-y^2)/Cwid^2;
      elseif(iprob==5)
         grav = 9.81d0;
         Cwid = 1.d0;
         Len = 5.d0;
         x2 = x-Len;
         s_xx = rho*grav*y - 1.d0/10.d0*rho*grav*y*(15.d0*Len^2 ...
                    -15.d0*x2^2 +4.d0*Cwid^2 +10.d0*y^2)/Cwid^2;
         s_xy = 3.d0/2.d0*rho*grav*x2*(Cwid^2-y^2)/Cwid^2;
         s_yy = rho*grav*y - 1.d0/2.d0*rho*grav*y*(3.d0*Cwid^2 ...
                    -y^2)/Cwid^2;
      elseif(iprob==2)
          s_xx = 72*y;
          s_xy = 0;
          s_yy = 0;
      elseif(iprob==7)
          lam = 0.54448373678246398;
          Q = 0.54307557883673652;
	      r = sqrt(x^2+y^2);
	     if(y>=0.d0)
           theta = acos(x/r);
         else
	       theta = -acos(x/r);
         end
	     r = lam*r^(lam-1.d0);
         s_xx = r*((2.d0 - Q*(lam+1.d0))*cos((lam-1.d0)*theta) ...
              - (lam-1.d0)*cos((lam-3.d0)*theta));
         s_yy = r*((2.d0 + Q*(lam+1.d0))*cos((lam-1.d0)*theta) ...
              + (lam-1.d0)*cos((lam-3.d0)*theta));
         s_xy = r*((lam-1.d0)*sin((lam-3.d0)*theta) ...
              + Q*(lam+1.d0)*sin((lam-1.d0)*theta));
      end

	   rx = s_xx*nv(1) + s_xy*nv(2);
	   ry = s_xy*nv(1) + s_yy*nv(2);
            
%.....---------------------------------------------------
%     Loop over nodes to compute contributions to rows of 
%     the stiffness matrix
%.....---------------------------------------------------

            for j=1:nel
               djn=shg(3,j)*det*w;

               eF(ndf*j-2) = eF(ndf*j-2) + djn*rx*uno;
               eF(ndf*j-1) = eF(ndf*j-1) + djn*ry*uno;

            end %j rows
         end %l int points

         end
         end % Residual form


      end
