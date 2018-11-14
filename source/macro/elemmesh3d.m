function [xc,ixc] = elemmesh3d(m,xe,ixe,ndm,nel,nen,mat,cel,celn,ncel)
%
%...  Written by Tim Truster (Fall 2009)
%...  Program to compute cell nodes and connectivity; program calls
%     FEAP routines for generation surfaces (blocks) of elements;
%     for triangles, the element is extended to an equivalent
%     quadrilateral with the intention of ignoring half of the
%     generated data
%
%**********************************************************************

%       implicit none
% 
%       logical debug
%       common /debugs/ debug
% 
% %     Input Variables
%       integer m,ixe(*),ndm,nel,nen,mat,cel,celn,ncel
%       real*8 xe(ndm,*)
% 
% %     Output Variables
%       real*8 xc(ndm,celn)
% 	integer ixc(nen+1,celn)
% 
% %     Local Variables
%       integer	nr,ns,ixl(9),no,ntyp,nm,ilr(2)
% 	integer i,j
% 	real*8 xl(3,9),shp(3,9)
% 	real*8 dr,ds
% 	character ctype*15
% 	logical dprt

  ixl = zeros(27,1);
  ilr = zeros(2,1);
  xl = zeros(3,27);
  shp = zeros(3,27);
  xc = zeros(0);
  ixc = zeros(0);

% 	ctype = 'cart'
% 	do j = 1,celn
% 	do i = 1,nen+1
% 	ixc(i,j) = 0
% 	enddo
% 	enddo
% 
% %Set up local element array in proper format for FEAP subroutines
% 	if(nel.eq.3) then
% 
% %  Create equivalent 4-node quadrilateral
% 	   ntyp = 2 %for triangles with diagonals ul to lr
% 	   nm = 4
% 	   do i = 1,2
% 	      do j = 1,2
% 	         xl(j,i) = xe(j,i)
% 	      enddo
% 	      xl(3,i) = 0.d0
% 	      ixl(i) = i
% 	   enddo
% %  Store third vertex as fourth corner
% 	   do j = 1,2
% 	      xl(j,4) = xe(j,3)
% 	   enddo
% 	   xl(3,4) = 0.d0
% 	   ixl(4) = 4
% %  Compute equivalent parallelogram corner
% 	   do j = 1,2
% 	      xl(j,3) = xe(j,2) + xe(j,3) - xe(j,1)
% 	   enddo
% 	   xl(3,3) = 0.d0
% 	   ixl(3) = 3
% %  Zero remaining node flags
% 	   do i = 5,9
% 	      ixl(i) = 0
% 	   enddo
% 
% %  Compute increments for cell generation
% 	   nr = m+1
% 	   ns = m+1
% 	   dr = 2.d0/m
% 	   ds = 2.d0/m
% 
% 	elseif(nel.eq.4) then
% 
% %  Transfer element values into local array
% 	   ntyp = 0
% 	   nm = 4
% 	   do i = 1,4
% 	      do j = 1,2
% 	         xl(j,i) = xe(j,i)
% 	      enddo
% 	      xl(3,i) = 0.d0
% 	      ixl(i) = i
% 	   enddo
% %  Zero remaining node flags
% 	   do i = 5,9
% 	      ixl(i) = 0
% 	   enddo
% 
% %  Compute increments for cell generation
% 	   nr = m+1
% 	   ns = m+1
% 	   dr = 2.d0/m
% 	   ds = 2.d0/m

%%	
%     ctype = 'cart';
% 	do j = 1,celn
% 	do i = 1,nen+1
% 	ixc(i,j) = 0
% 	enddo
% 	enddo

%Set up local element array in proper format for FEAP subroutines
	if(nel==3)

%  Create equivalent 4-node quadrilateral
	   ntyp = 2; %for triangles with diagonals ul to lr
	   nm = 4;
       for i = 1:2
          for j = 1:2
	         xl(j,i) = xe(j,i);
          end
	      xl(3,i) = 0.d0;
	      ixl(i) = i;
       end
%  Store third vertex as fourth corner
       for j = 1:2
	      xl(j,4) = xe(j,3);
       end
	   xl(3,4) = 0.d0;
	   ixl(4) = 4;
%  Compute equivalent parallelogram corner
       for j = 1:2
	      xl(j,3) = xe(j,2) + xe(j,3) - xe(j,1);
       end
	   xl(3,3) = 0.d0;
	   ixl(3) = 3;
%  Zero remaining node flags
       for i = 5:9
	      ixl(i) = 0;
       end

%  Compute increments for cell generation
	   nr = m+1;
	   ns = m+1;
	   dr = 2.d0/m;
	   ds = 2.d0/m;

% 	elseif(nel.eq.4) then
% 
% %  Transfer element values into local array
% 	   ntyp = 0
% 	   nm = 4
% 	   do i = 1,4
% 	      do j = 1,2
% 	         xl(j,i) = xe(j,i)
% 	      enddo
% 	      xl(3,i) = 0.d0
% 	      ixl(i) = i
% 	   enddo
% %  Zero remaining node flags
% 	   do i = 5,9
% 	      ixl(i) = 0
% 	   enddo
% 
% %  Compute increments for cell generation
% 	   nr = m+1
% 	   ns = m+1
% 	   dr = 2.d0/m
% 	   ds = 2.d0/m
%        
% 	elseif(nel.eq.6) then
% 
% %  Transfer element values into local array
% 	   ntyp	= 7
% 	   nm = 9
% 	   do i = 1,2
% 	      do j = 1,2
% 	         xl(j,i) = xe(j,i)
% 	      enddo
% 	      xl(3,i) = 0.d0
% 	      ixl(i) = i
% 	   enddo
% %  Store third vertex as fourth corner
% 	   do j = 1,2
% 	      xl(j,4) = xe(j,3)
% 	   enddo
% 	   xl(3,4) = 0.d0
% 	   ixl(4) = 4
% %  Compute equivalent parallelogram corner
% 	   do j = 1,2
% 	      xl(j,3) = xe(j,2) + xe(j,3) - xe(j,1)
% 	   enddo
% 	   xl(3,3) = 0.d0
% 	   ixl(3) = 3
% 	   do j = 1,2
% 	      xl(j,5) = xe(j,4)
% 	   enddo
% 	   xl(3,5) = 0.d0
% 	   ixl(5) = 5
% 	   do j = 1,2
% 	      xl(j,6) = (xl(j,2) + xl(j,3))/2.d0
% 	   enddo
% 	   xl(3,6) = 0.d0
% 	   ixl(6) = 6
% 	   do j = 1,2
% 	      xl(j,7) = (xl(j,3) + xl(j,4))/2.d0
% 	   enddo
% 	   xl(3,7) = 0.d0
% 	   ixl(7) = 7
% 	   do j = 1,2
% 	      xl(j,8) = xe(j,6)
% 	   enddo
% 	   xl(3,8) = 0.d0
% 	   ixl(8) = 8
% 	   do j = 1,2
% 	      xl(j,9) = xe(j,5)
% 	   enddo
% 	   xl(3,9) = 0.d0
% 	   ixl(9) = 9
% 
% %  Compute increments for cell generation
% 	   nr = 2*m+1
% 	   ns = 2*m+1
% 	   dr = 1.d0/m
% 	   ds = 1.d0/m

%%
	elseif(nel==8)

%  Transfer element values into local array
	   ntyp = 10;
       for i = 1:8
          for j = 1:3
	         xl(j,i) = xe(j,i);
          end
	      ixl(i) = i;
       end
%  Zero remaining node flags
       for i = 9:27
	      ixl(i) = 0;
       end

%  Compute increments for cell generation
	   nr = m+1;
	   ns = m+1;
       nt = m+1;
	   dr = 2.d0/m;
	   ds = 2.d0/m;
	   dt = 2.d0/m;
       
	elseif(nel==10)

%  Transfer element values into local array
	   ntyp	= 7;
	   nm = 9;
       for i = 1:2
          for j = 1:2
	         xl(j,i) = xe(j,i);
          end
	      xl(3,i) = 0.d0;
	      ixl(i) = i;
       end
%  Store third vertex as fourth corner
       for j = 1:2
	      xl(j,4) = xe(j,3);
       end
	   xl(3,4) = 0.d0;
	   ixl(4) = 4;
%  Compute equivalent parallelogram corner
       for j = 1:2
	      xl(j,3) = xe(j,2) + xe(j,3) - xe(j,1);
       end
	   xl(3,3) = 0.d0;
	   ixl(3) = 3;
       for j = 1:2
	      xl(j,5) = xe(j,4);
       end
	   xl(3,5) = 0.d0;
	   ixl(5) = 5;
       for j = 1:2
	      xl(j,6) = (xl(j,2) + xl(j,3))/2.d0;
       end
	   xl(3,6) = 0.d0;
	   ixl(6) = 6;
       for j = 1:2
	      xl(j,7) = (xl(j,3) + xl(j,4))/2.d0;
       end
	   xl(3,7) = 0.d0;
	   ixl(7) = 7;
       for j = 1:2
	      xl(j,8) = xe(j,6);
       end
	   xl(3,8) = 0.d0;
	   ixl(8) = 8;
       for j = 1:2
	      xl(j,9) = xe(j,5);
       end
	   xl(3,9) = 0.d0;
	   ixl(9) = 9;

%  Compute increments for cell generation
	   nr = 2*m+1;
	   ns = 2*m+1;
	   dr = 1.d0/m;
	   ds = 1.d0/m;

% 	elseif(nel.eq.9) then
% 
% %  Transfer element values into local array
% 	   ntyp	= 9
% 	   nm = 9
% 	   do i = 1,9
% 	      do j = 1,2
% 	         xl(j,i) = xe(j,i)
% 	      enddo
% 	      xl(3,i) = 0.d0
% 	      ixl(i) = i
% 	   enddo
% 
% %  Compute increments for cell generation
% 	   nr = 2*m+1
% 	   ns = 2*m+1
% 	   dr = 1.d0/m
% 	   ds = 1.d0/m
% 
% 	endif
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	WRITE(24,*) 'nr',nr
% 	WRITE(24,*) 'ns',ns
% 	WRITE(24,*) 'dr',dr
% 	WRITE(24,*) 'ds',ds
% 	WRITE(24,*) 'ntyp',ntyp
% 	WRITE(24,*) 'ctype',ctype
% 	WRITE(24,*) 'nm',nm
% 	WRITE(24,*) 'xl'
% 	do i=1,9	
% 	   WRITE(24,*) (xl(j,i),j=1,3)    
% 	enddo
% 	write(24,*) 'ixl'
% 	do i=1,9
% 	   WRITE(24,*) ixl(i)
% 	enddo
% 	endif
% 	endif
% 
% %Compute cell nodes
% 	call sblkn(nr,ns,xl,ixl,shp,xc,dr,ds,1,no,ndm,0,ntyp,
%      &                 nm,ctype,.false.,.false.)
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then 
% 	WRITE(24,*) 'no',no
% 	WRITE(24,*) 'xc'
% 	do j=1,ncel	
% 		
% 		   WRITE(24,*) (xc(i,j),i=1,ndm)    
% 
% 	enddo
% 	endif
% 	endif
% 
% %Compute cell connectivity
%       if(nel.eq.6) then
% 
% 	call sblke6(nr,ns,xc,ixc,1,1,no,ndm,nen+1,0,ntyp,nm,mat,
%      &                 0,ilr,ctype)
% 
% 	else
% 
% 	call sblke(nr,ns,xc,ixc,1,1,no,ndm,nen+1,0,ntyp,nm,mat,
%      &                 0,ilr,ctype)
% 
% 	endif
% 
% 	if(debug) then
% 	dprt = .true.
% 	if(dprt) then
% 	WRITE(24,*) 'ixc'
% 	do j=1,ncel	
% 		
% 		   WRITE(24,*) (ixc(i,j),i=1,nen+1)    
% 
% 	enddo
% 	write(24,'(A)') ''
% 	endif
% 	endif

%%
	elseif(nel==27)

%  Transfer element values into local array
	   ntyp	= 12;
% 	   for i = 1:27
%           for j = 1:3
% 	         xl(j,i) = xe(j,i);
%           end
% 	      ixl(i) = i;
% 	   end
       for i = 1:8
          for j = 1:3
	         xl(j,i) = xe(j,i);
          end
	      ixl(i) = i;
       end
%  Zero remaining node flags
       for i = 9:27
	      ixl(i) = 0;
       end

%  Compute increments for cell generation
	   nr = 2*m+1;
	   ns = 2*m+1;
	   nt = 2*m+1;
	   dr = 1.d0/m;
	   ds = 1.d0/m;
	   dt = 1.d0/m;

	end
    
    nr1 = nr - 1;
    ns1 = ns - 1;
    nt1 = nt - 1;
    if(ntyp==10)     % Linear    hexahedron
      nf = nr1*ns1*nt1 - 1;
    elseif(ntyp==11) % Linear    tetrahedron
      nf = 6*nr1*ns1*nt1 - 1;
    elseif(ntyp==12) % Quadratic hexahedron
      nf = (nr1*ns1*nt1)/8 - 1;
    elseif(ntyp==13) % Quadratic tetrahedron
      nf = 6*(nr1*ns1*nt1)/8 - 1;
    end
    numel11 = nf;
    ng1 = nr*ns*nt -1;
    numnp11 = ng1;
    
%Compute cell nodes
     xc = vblkn(nr,ns,nt,xl,xc,ixl,dr,ds,dt,1,ndm,'cart',numnp11,[0 0 0]);

%Compute cell connectivity
    ixc = vblke(nr,ns,nt,ixc,1,1,nen+1,mat,ntyp,numel11);
