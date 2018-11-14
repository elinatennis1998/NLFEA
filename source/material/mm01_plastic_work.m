function [cgn1] = mm01_plastic_work(cgn,cgn1,yield,deps,nu,ym,shear_mod)
% c     ****************************************************************
% c     *                                                              *
% c     *                    subroutine mm01_plastic_work              *
% c     *                                                              *
% c     *                       written by : rhd                       *
% c     *                                                              *
% c     *                   last modified: 10/1/00                     *
% c     *                                                              *
% c     *                                                              *
% c     ****************************************************************
% c
one = 1.d0; two = 2.d0; three = 3.d0; %zero = 0.d0; 
half = one/two; root2 = sqrt(two);

% c
% c                       for plastic points, compute the updated
% c                       plastic work density and accumulated
% c                       (uniaxial) plastic strain
% c
% c
%       do i = 1, span
       cgn1(8)        = cgn(8);
       cgn1(9)        = cgn(9);
       if ( yield ) %then  
%          deps_plas_bar    = zero;
         dsig(1:6) = cgn1(1:6) - cgn(1:6);
%          dsig(2) = cgn1(i,2) - cgn(i,2)
%          dsig(3) = cgn1(i,3) - cgn(i,3)
%          dsig(4) = cgn1(i,4) - cgn(i,4)
%          dsig(5) = cgn1(i,5) - cgn(i,5)
%          dsig(6) = cgn1(i,6) - cgn(i,6)
         deps_plas(1) = deps(1) - ...
                        (dsig(1) - nu*(dsig(2)+dsig(3)))/ym;
         deps_plas(2) = deps(2) - ...
                        (dsig(2) - nu*(dsig(1)+dsig(3)))/ym;
         deps_plas(3) = deps(3) - ...
                        (dsig(3) - nu*(dsig(1)+dsig(2)))/ym;
         deps_plas(4) = deps(4) - dsig(4) / shear_mod;
         deps_plas(5) = deps(5) - dsig(5) / shear_mod;
         deps_plas(6) = deps(6) - dsig(6) / shear_mod;
         cgn1(8) = cgn(8) +  half * ( ...
             deps_plas(1) * (cgn1(1) + cgn(1)) ...
           + deps_plas(2) * (cgn1(2) + cgn(2)) ...
           + deps_plas(3) * (cgn1(3) + cgn(3)) ...
           + deps_plas(4) * (cgn1(4) + cgn(4)) ...
           + deps_plas(5) * (cgn1(5) + cgn(5)) ...
           + deps_plas(6) * (cgn1(6) + cgn(6)) );
         factor1 = ( deps_plas(1) - deps_plas(2) )^2  + ...
                   ( deps_plas(2) - deps_plas(3) )^2  + ...
                   ( deps_plas(1) - deps_plas(3) )^2 ;
         factor2 = deps_plas(4)^2 +  deps_plas(5)^2 + ...
                   deps_plas(6)^2;
         deps_plas_bar =  (root2/three) * sqrt( factor1 +  ...
                          (three/two)*factor2 ) ;  
         cgn1(9) = cgn(9) + deps_plas_bar;
        end %if
%       end do 