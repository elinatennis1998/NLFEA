function eK = getKcellNL5(xl,d,xe,ue,MR,BR,MS,BS,ndf,ndm,lint,nel,nen,nummat,symmns)
%
%...  Written by Tim Truster (Fall 2009)
%     No modifications when copied to NLFEA ver2
%...  Program to compute stiffness matrix for cell of submesh
%     Stiffness matrix eK is formed using fine scale equation for
%     ndf = 2, and using stabilized formulation for ndf = 3.
%
%**********************************************************************
% 
%       Implicit None
% 
%       logical debug
%       common /debugs/ debug
% 
% %     Input Variables
%       integer ndf,ndm,lint,nel,nen,nummat
%       real*8 d(*),xl(ndm,nen)
% 
% %     Output Variables
%       real*8 eK(nen*ndf,nen*ndf)
% 
% %     Local Variables
%       integer i,j,l
%       integer iprob,IQTY
% 
%       real*8 w,det,
%      >       djx,djy,djn,djxx,djyy,djxy,
%      >       dix,diy,din,dixx,diyy,dixy,       
%      >       t11,t12,t21,t22,f11,f12,f21,f22,b,r,s,
%      >       shl(3,10),shls(3,10),shg(3,10),shgs(3,10)
%       real*8 bubble(3),sx(2,2)
%       integer ib
%       logical der,bf,dprt
% 
%       integer i_cond_field,Lsym,xinc,yinc,minc
%       real*8 mu,a1,a2,eps,rho,lam,young,nu

eK = zeros(nen*ndf,nen*ndf);

% cellE = d(2);
% cellv = d(3);
% cellrho = d(4);
% PSPS = d(5);
% % iprob = d(6);
% if PSPS == 0 %plane stress
%     cellE = cellE*(1 + 2*cellv)/(1 + cellv)^2;
%     cellv = cellv/(1 + cellv);
%     lam = cellv*cellE/((1+cellv)*(1-2*cellv));
%     mu = cellE/(2*(1+cellv));
% else %plane strain
%     lam = cellv*cellE/((1+cellv)*(1-2*cellv));
%     mu = cellE/(2*(1+cellv));
% end
% % Get material parameters
%       Call M_Prop(d,mu,a1,a2,eps,rho,Iprob,i_cond_field,Lsym,
%      *	        lam,young,nu,IQTY,xinc,yinc,minc)
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
%          write(24,*) 'mu',mu
%          write(24,*) 'lam',lam
%          write(24,*) 'young',young
%          write(24,*) 'nu',nu
% 	   write(24,'(A)') ''
% 	endif
% 	endif

% % Determine formulation type of cell
% 
%       if(ndf.eq.2) then  % Fine scale type cell
% 
% %   Initialize parameters
% 
%          ib = 0
% 	   der = .false.
% 	   bf = .false.
% 
%          call pzero(eK,nen*ndf*nen*ndf)
% 
% %.....---------------------------------------------------
% %     Integration of matrix terms over cell
% %.....---------------------------------------------------
%          do 250 l=1,lint
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call intpntt(l,lint,ib,w,r,s)
%               call shlt(r,s,nel,der,bf,shl,shls,bubble)
%               call shgt(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call intpntq(l,lint,ib,w,r,s)
%               call shlq(r,s,nel,der,bf,shl,shls,bubble)
%               call shgq(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             endif
% 
% 	if(debug) then
% 	dprt = .false.
% 	if(dprt) then
%          write(24,*) 'det',det
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
% %     Loop over weighting functions of cell to compute  
% %     contributions to rows of the stiffness matrix
% %.....---------------------------------------------------
% 
%             do 200 j=1,nel
%                djx=shg(1,j)*det*w
%                djy=shg(2,j)*det*w
% 
% %.....---------------------------------------------------
% %     Loop over fine scale functions of cell to compute  
% %     contributions to columns of the stiffness matrix
% %.....---------------------------------------------------   
% 
%             do 210 i=1,nel
%                dix=shg(1,i)
%                diy=shg(2,i)
% 
% %* (1,1)*************************************************************************************
%                eK(ndf*i-1,ndf*j-1) = eK(ndf*i-1,ndf*j-1) 
%      &         + mu*(2.d0*dix*djx + diy*djy)
% 
% %* (1,2)*************************************************************************************
%                eK(ndf*i-1,ndf*j)   = eK(ndf*i-1,ndf*j) 
%      &         + mu*diy*djx
% 
% %* (2,1)*************************************************************************************
%                eK(ndf*i,ndf*j-1)   = eK(ndf*i,ndf*j-1)
%      &         + mu*dix*djy
%      				    
% %* (2,2)*************************************************************************************
%                eK(ndf*i,ndf*j)     = eK(ndf*i,ndf*j)
%      &         + mu*(dix*djx + 2.d0*diy*djy)
% 
% %          -----------------------------------------------------
% %           End loop over fine scale functions
% %          -----------------------------------------------------
%  210        continue
% %        -----------------------------------------------------
% %         End loop over weighting functions
% %        -----------------------------------------------------
%  200        continue
% 
% %	if(debug) then
% %	dprt = .true.
% %	if(dprt) then
% %        write(24,*) 'eK11',eK(1,2)
% %	   write(24,'(A)') ''
% %	endif
% %	endif
% %      -----------------------------------------------------
% %       End loop over integration points
% %      -----------------------------------------------------
%  250     continue

%%
% Determine formulation type of cell

      if(ndf==2)  % Fine scale type cell

%   Initialize parameters

         ib = 0;
	   der = 0;
	   bf = 0;
        
        I1 = [1; 1; 0; 0];
        spvec0 = I1;
        cpmat1 = I1*I1';
        cpmat2 = diag([-2,-2,-1,0]);
        cpmat0 = cpmat1 + cpmat2;

%.....---------------------------------------------------
%     Integration of matrix terms over cell
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
              [shel,shed,shels,bubble] = shlt(bigR,bigS,nel,nel,der,bf);
              [she] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end

            [F,JxX,fi] = kine2d(she,ue,nel,0);
            Qxy = shg*fi;
            if symmns == 1 %NS
            spvec = JxX*spvec0;
            cpmat = JxX*cpmat0;
            press = ue(3,:)*shel;
            else %S
            [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
            spvec = JxX*theta1*spvec0;
            cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
            press = ue(3,:)*shel;
            end
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
            
            if nel == 3            
            % Form B matrix    
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1)];
            elseif nel == 4            
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1)];
            elseif nel == 6
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0        Qxy(5,1)  0        Qxy(6,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 0         Qxy(5,2) 0         Qxy(6,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1) Qxy(5,2)  Qxy(5,1) Qxy(6,2)  Qxy(6,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1) Qxy(5,2) -Qxy(5,1) Qxy(6,2) -Qxy(6,1)];
            elseif nel == 9
            % Form B matrix
            Bmat = [Qxy(1,1)  0        Qxy(2,1)  0        Qxy(3,1)  0        Qxy(4,1)  0        Qxy(5,1)  0        Qxy(6,1)  0        Qxy(7,1)  0        Qxy(8,1)  0        Qxy(9,1)  0       
                    0         Qxy(1,2) 0         Qxy(2,2) 0         Qxy(3,2) 0         Qxy(4,2) 0         Qxy(5,2) 0         Qxy(6,2) 0         Qxy(7,2) 0         Qxy(8,2) 0         Qxy(9,2)
                    Qxy(1,2)  Qxy(1,1) Qxy(2,2)  Qxy(2,1) Qxy(3,2)  Qxy(3,1) Qxy(4,2)  Qxy(4,1) Qxy(5,2)  Qxy(5,1) Qxy(6,2)  Qxy(6,1) Qxy(7,2)  Qxy(7,1) Qxy(8,2)  Qxy(8,1) Qxy(9,2)  Qxy(9,1)
                    Qxy(1,2) -Qxy(1,1) Qxy(2,2) -Qxy(2,1) Qxy(3,2) -Qxy(3,1) Qxy(4,2) -Qxy(4,1) Qxy(5,2) -Qxy(5,1) Qxy(6,2) -Qxy(6,1) Qxy(7,2) -Qxy(7,1) Qxy(8,2) -Qxy(8,1) Qxy(9,2) -Qxy(9,1)];
            end
            
            eK = eK + det*w*Bmat'*Dmat*Bmat;
            
% %.....---------------------------------------------------
% %     Loop over weighting functions of cell to compute  
% %     contributions to rows of the stiffness matrix
% %.....---------------------------------------------------
% 
%             for j=1:nel
%                djx=shg(1,j)*det*w;
%                djy=shg(2,j)*det*w;
% 
% %.....---------------------------------------------------
% %     Loop over fine scale functions of cell to compute  
% %     contributions to columns of the stiffness matrix
% %.....---------------------------------------------------   
% 
%             for i=1:nel
%                dix=shg(1,i);
%                diy=shg(2,i);
% 
% %* (1,1)*************************************************************************************
%                eK(ndf*i-1,ndf*j-1) = eK(ndf*i-1,ndf*j-1) + mu*(2.d0*dix*djx + diy*djy);
% 
% %* (1,2)*************************************************************************************
%                eK(ndf*i-1,ndf*j)   = eK(ndf*i-1,ndf*j) + mu*diy*djx;
% 
% %* (2,1)*************************************************************************************
%                eK(ndf*i,ndf*j-1)   = eK(ndf*i,ndf*j-1) + mu*dix*djy;
%      				    
% %* (2,2)*************************************************************************************
%                eK(ndf*i,ndf*j)     = eK(ndf*i,ndf*j) + mu*(dix*djx + 2.d0*diy*djy);
% 
% %          -----------------------------------------------------
% %           End loop over fine scale functions
% %          -----------------------------------------------------
%             end
% %        -----------------------------------------------------
% %         End loop over weighting functions
% %        -----------------------------------------------------
%             end

%	if(debug) then
%	dprt = .true.
%	if(dprt) then
%        write(24,*) 'eK11',eK(1,2)
%	   write(24,'(A)') ''
%	endif
%	endif
%      -----------------------------------------------------
%       End loop over integration points
%      -----------------------------------------------------
         end


      else %Stabilized cell type

%          %   Initialize parameters
% 
%          ib = 0
% 	     der = .true.
% 	     bf = .true.
% 
%          call pzero(eK,nen*ndf*nen*ndf)
%          
%          %.....Stabilization Tensor Tau (Consistent Formulation) 5.25.08
% 
% 	     Call Tau(xl,mu,nel,ndm,t11,t12,t21,t22,f11,f12,f21,f22)
% 
% %.....---------------------------------------------------
% %     Integration of matrix and vector terms over element
% %.....---------------------------------------------------
% 	
%          do 350 l=1,lint
% 
%             if(nel.eq.3.or.nel.eq.6) then
%               call intpntt(l,lint,ib,w,r,s)
%               call shlt(r,s,nel,der,bf,shl,shls,bubble)
%               call shgt(xl,nel,shl,shls,nummat,nen,bf,der,det,shg,shgs,
%      &                bubble,sx)
%             elseif(nel.eq.4.or.nel.eq.9) then
%               call intpntq(l,lint,ib,w,r,s)
%               call shlq(r,s,nel,der,bf,shl,shls,bubble)
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
%             do 300 j=1,nel
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
%             do 310 i=1,nel
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
%      &    - mu*mu*(2.d0*dixx*(2.d0*djxx*b*t11+djyy*b*t11+djxy*b*t12) 
%      &    + diyy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) 
%      &    + dixy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) 
% %Term20x:Standard Diffusion
%      &    + mu*(2.d0*dix*djx + diy*djy)
% 
% %* (1,2)*************************************************************************************
%             eK(ndf*i-2,ndf*j-1)=eK(ndf*i-2,ndf*j-1)  
% %Term6x:Diffusion-Diffusion
%      &    - mu*mu*(diyy*(djxy*b*t11 + djxx*b*t12 + 2.d0*djyy*b*t12) 
%      &    + 2.d0*dixx*(djxy*b*t11 + (djxx+2.d0*djyy)*b*t12) 
%      &    + dixy*(djxy*b*t21 + djxx*b*t22 + 2.d0*djyy*b*t22)) 
% %Term20x:Standard Diffusion
%      &    + mu*diy*djx
% 
% %* (1,3)*************************************************************************************
%             eK(ndf*i-2,ndf*j  )=eK(ndf*i-2,ndf*j  ) 
% %Term7x:Diffusion-Pressure
%      &    - mu*(diyy*djx*b*t11 + diyy*djy*b*t12 + 2.d0*dixx*(djx*b*t11 
%      &    + djy*b*t12) + dixy*djx*b*t21 +  dixy*djy*b*t22) 	
% %Term21x:Standard Pressure
%      &    + dix*djn
%       
% %* (2,1)*************************************************************************************
%             eK(ndf*i-1,ndf*j-2)=eK(ndf*i-1,ndf*j-2)  
% %Term6y:Diffusion-Diffusion
%      &    - mu*mu*(dixy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) 
%      &    + (dixx + 2.d0*diyy)*(2.d0*djxx*b*t21+djyy*b*t21+djxy*b*t22))     
% %Term20y:Standard Diffusion
%      &    + mu*dix*djy
%      				    
% %* (2,2)*************************************************************************************
%             eK(ndf*i-1,ndf*j-1)=eK(ndf*i-1,ndf*j-1) 
% %Term6y:Diffusion-Diffusion
%      &    - mu*mu*(dixy*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) 
%      &    + (dixx + 2.d0*diyy)*(djxy*b*t21 + (djxx + 2.d0*djyy)*b*t22)) 
% %Term20y:Standard Diffusion
%      &    + mu*(dix*djx + 2.d0*diy*djy)
% 
% %* (2,3)*************************************************************************************
%             eK(ndf*i-1,ndf*j  )=eK(ndf*i-1,ndf*j  ) 
% %Term7y:Diffusion-Pressure
%      &    - mu*(dixy*(djx*b*t11 + djy*b*t12) 
%      &    + (dixx + 2.d0*diyy)*(djx*b*t21 + djy*b*t22)) 
% %Term21y:Standard Pressure
%      &    + diy*djn
%      			 
% %* (3,1)*************************************************************************************
%             eK(ndf*i  ,ndf*j-2)=eK(ndf*i  ,ndf*j-2) 
% %Term5:Continuity-Diffusion
%      &    - mu*(dix*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) 
%      &    + diy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) 
% %Term22x:Standard Continuity
%      &    + din*djx
%        
% %* (3,2)*************************************************************************************
%             eK(ndf*i  ,ndf*j-1)=eK(ndf*i  ,ndf*j-1) 
% %Term5:Continuity-Diffusion
%      &    - mu*(dix*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) 
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
%  310        continue %i columns
%  300        continue %j rows
%  350     continue %l int points
% 
%       endif
% 
%       return
% 
%       end
% 

%%
         %   Initialize parameters

         ib = 0;
	     der = 1;
	     bf = 1;
         
         %.....Stabilization Tensor Tau (Consistent Formulation) 5.25.08

	     [t11,t12,t21,t22] = Tau(xl,mu,nel,nen,lint);

%.....---------------------------------------------------
%     Integration of matrix and vector terms over element
%.....---------------------------------------------------
	
         for l=1:lint

            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,lint,ib);
              [shl,shls,bubble] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, bubble] = shgt(xl,nel,shl,shls,nen,bf,der,bubble);
            else
              [w,litr,lits] =  intpntq(l,lint,ib);
              [shl,shls,bubble] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, bubble] = shgq(xl,nel,shl,shls,nen,bf,der,bubble);
            end
	        b = bubble(3);

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
          - mu*mu*(2.d0*dixx*(2.d0*djxx*b*t11+djyy*b*t11+djxy*b*t12) ...
          + diyy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) ...
          + dixy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) ...
          + mu*(2.d0*dix*djx + diy*djy);

%* (1,2)*************************************************************************************
            eK(ndf*i-2,ndf*j-1)=eK(ndf*i-2,ndf*j-1) ...
          - mu*mu*(diyy*(djxy*b*t11 + djxx*b*t12 + 2.d0*djyy*b*t12) ...
          + 2.d0*dixx*(djxy*b*t11 + (djxx+2.d0*djyy)*b*t12) ...
          + dixy*(djxy*b*t21 + djxx*b*t22 + 2.d0*djyy*b*t22)) ...
          + mu*diy*djx;

%* (1,3)*************************************************************************************
            eK(ndf*i-2,ndf*j  )=eK(ndf*i-2,ndf*j  ) ...
          - mu*(diyy*djx*b*t11 + diyy*djy*b*t12 + 2.d0*dixx*(djx*b*t11 ...
          + djy*b*t12) + dixy*djx*b*t21 +  dixy*djy*b*t22) ...
          + dix*djn;
      
%* (2,1)*************************************************************************************
            eK(ndf*i-1,ndf*j-2)=eK(ndf*i-1,ndf*j-2) ...
          - mu*mu*(dixy*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) ...
          + (dixx + 2.d0*diyy)*(2.d0*djxx*b*t21+djyy*b*t21+djxy*b*t22)) ...
          + mu*dix*djy;
     				    
%* (2,2)*************************************************************************************
            eK(ndf*i-1,ndf*j-1)=eK(ndf*i-1,ndf*j-1) ...
          - mu*mu*(dixy*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) ...
          + (dixx + 2.d0*diyy)*(djxy*b*t21 + (djxx + 2.d0*djyy)*b*t22)) ...
          + mu*(dix*djx + 2.d0*diy*djy);

%* (2,3)*************************************************************************************
            eK(ndf*i-1,ndf*j  )=eK(ndf*i-1,ndf*j  ) ...
          - mu*(dixy*(djx*b*t11 + djy*b*t12) ...
          + (dixx + 2.d0*diyy)*(djx*b*t21 + djy*b*t22)) ...
          + diy*djn;
     			 
%* (3,1)*************************************************************************************
            eK(ndf*i  ,ndf*j-2)=eK(ndf*i  ,ndf*j-2) ...
          - mu*(dix*(2.d0*djxx*b*t11 + djyy*b*t11 + djxy*b*t12) ...
          + diy*(2.d0*djxx*b*t21 + djyy*b*t21 + djxy*b*t22)) ...
          + din*djx;
       
%* (3,2)*************************************************************************************
            eK(ndf*i  ,ndf*j-1)=eK(ndf*i  ,ndf*j-1) ...
          - mu*(dix*(djxy*b*t11 + (djxx + 2.d0*djyy)*b*t12) ...
          + diy*(djxy*b*t21 + (djxx + 2.d0*djyy)*b*t22)) ...
          + din*djy;

%* (3,3)*************************************************************************************
            eK(ndf*i  ,ndf*j  )=eK(ndf*i  ,ndf*j  ) ...
          -dix*djx*b*t11 - dix*djy*b*t12 - diy*djx*b*t21 - diy*djy*b*t22 ...
          -din*djn*(1.d0/lam);

            end
            end
         end

      end

