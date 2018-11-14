function cep = cnst1(rtsg,nu,e,kn1,hprime,beta,ldt,dj,w,yield)
% c     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
% c     *                                                              *
% c     *                      subroutine cnst1                        *
% c     *                                                              *
% c     *                       written by : bh, rhd                   *
% c     *                                                              *
% c     *                   last modified : 06/29/91, 10/8/97          *
% c     *                                                              *
% c     *      this subroutine computes the consistent                 *
% c     *      tangent operator matrix for simple rate independent     *
% c     *      mises plasticity with constant (mixed) hardening        *
% c     *      all elements in the block are processed in full         *
% c     *      vector mode. this is our fastest material model.        *
% c     *                                                              *
% c     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% c
% c
cep = zeros(6,6);
one = 1.d0; two = 2.d0; three = 3.d0; root2 = sqrt(two); zero = 0.d0; 
%       do i = 1, span
% c           
% c                       branch on whether or not the material
% c                       has yielded.
% c
         if( ~yield ) %then
% c
% c                       linear isotropic elastic matrix. there
% c                       is no yield.
% c
            cep(1,4)= zero;
            cep(1,5)= zero;
            cep(1,6)= zero;
            cep(2,4)= zero;
            cep(2,5)= zero;
            cep(2,6)= zero;
            cep(3,4)= zero;
            cep(3,5)= zero;
            cep(3,6)= zero;
            cep(4,1)= zero;
            cep(4,2)= zero;
            cep(4,3)= zero;
            cep(4,5)= zero;
            cep(4,6)= zero;
            cep(5,1)= zero;
            cep(5,2)= zero;
            cep(5,3)= zero;
            cep(5,4)= zero;
            cep(5,6)= zero;
            cep(6,1)= zero;
            cep(6,2)= zero;
            cep(6,3)= zero;
            cep(6,4)= zero;
            cep(6,5)= zero;
            c1= (e/((one+nu)*(one-two*nu)))*dj*w;
            c2= (one-nu)*c1   ;
            c3= ((one-two*nu)/two)*c1;
            c4= nu*c1;
            cep(1,1)= c2;
            cep(2,2)= c2;
            cep(3,3)= c2;
            cep(4,4)= c3;
            cep(5,5)= c3;
            cep(6,6)= c3;
            cep(1,2)= c4;
            cep(1,3)= c4;
            cep(2,1)= c4;
            cep(3,1)= c4;
            cep(2,3)= c4;
            cep(3,2)= c4;
         else
% c
% c                       the material point is in the plastic
% c                       range. compute the consistent constituitive
% c                       matrix.
% c
            g = e/(two*(one+nu));
            l = (e*nu)/((one+nu)*(one-two*nu));
            k = (three*l+two*g)/three;                  
            mrtsq = rtsg(1)^2+rtsg(2)^2+rtsg(3)^2+two* ...
                      (rtsg(4)^2+rtsg(5)^2+rtsg(6)^2);
            bb = (root2*kn1+(two/three)*(one-beta)*hprime ...
                    *ldt)/sqrt(mrtsq);
            gamma= one/(one+hprime/(three*g));
            gambar=  gamma-one+bb;
            gbar= g*bb;
            albar= k- two*gbar/three;
            thbar= two*g*gambar;
            cep(1,1)=  (albar+two*gbar-thbar*(rtsg(1)^2)/ ...
                          mrtsq)*dj*w;
            cep(2,2)=  (albar+two*gbar-thbar*(rtsg(2)^2)/ ...
                          mrtsq)*dj*w;
            cep(3,3)=  (albar+two*gbar-thbar*(rtsg(3)^2)/ ...
                          mrtsq)*dj*w;
            cep(4,4)=  (gbar-thbar*(rtsg(4)^2)/mrtsq)* ...
                          dj*w;
            cep(5,5)=  (gbar-thbar*(rtsg(5)^2)/mrtsq)* ...
                          dj*w;
            cep(6,6)=  (gbar-thbar*(rtsg(6)^2)/mrtsq)* ...
                          dj*w;
            cep(2,1)=  (albar-thbar*rtsg(1)*rtsg(2)/ ...
                          mrtsq)*dj*w;
            cep(3,1)=  (albar-thbar*rtsg(1)*rtsg(3)/ ...
                          mrtsq)*dj*w;
            cep(4,1)= -(thbar*rtsg(1)*rtsg(4)/mrtsq)*dj*w;
            cep(5,1)= -(thbar*rtsg(1)*rtsg(5)/mrtsq)*dj*w;
            cep(6,1)= -(thbar*rtsg(1)*rtsg(6)/mrtsq)*dj*w;
            cep(3,2)=  (albar-thbar*rtsg(3)*rtsg(2)/ ...
                          mrtsq)*dj*w;
            cep(4,2)= -(thbar*rtsg(2)*rtsg(4)/mrtsq)*dj*w;
            cep(5,2)= -(thbar*rtsg(2)*rtsg(5)/mrtsq)*dj*w;
            cep(6,2)= -(thbar*rtsg(2)*rtsg(6)/mrtsq)*dj*w;
            cep(4,3)= -(thbar*rtsg(3)*rtsg(4)/mrtsq)*dj*w;
            cep(5,3)= -(thbar*rtsg(3)*rtsg(5)/mrtsq)*dj*w;
            cep(6,3)= -(thbar*rtsg(3)*rtsg(6)/mrtsq)*dj*w;
            cep(5,4)= -(thbar*rtsg(4)*rtsg(5)/mrtsq)*dj*w;
            cep(6,4)= -(thbar*rtsg(4)*rtsg(6)/mrtsq)*dj*w;
            cep(6,5)= -(thbar*rtsg(5)*rtsg(6)/mrtsq)*dj*w;
% c                  
            cep(1,2)=cep(2,1);
            cep(1,3)=cep(3,1);
            cep(1,4)=cep(4,1);
            cep(1,5)=cep(5,1);
            cep(1,6)=cep(6,1);
            cep(2,3)=cep(3,2);
            cep(2,4)=cep(4,2);
            cep(2,5)=cep(5,2);
            cep(2,6)=cep(6,2);
            cep(3,4)=cep(4,3);
            cep(3,5)=cep(5,3);
            cep(3,6)=cep(6,3);
            cep(4,5)=cep(5,4);
            cep(4,6)=cep(6,4);
            cep(5,6)=cep(6,5);
         end %if
%       end do
%       return
      end
