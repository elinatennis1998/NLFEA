function cep = ctran1(cep,qn1,cs,qbar,dj,w)
% c     ****************************************************************
% c     *                                                              *
% c     *                      subroutine ctran1                       *
% c     *                                                              *
% c     *                       written by : bh                        *
% c     *                                                              *
% c     *                   last modified : 06/18/12 rhd               *
% c     *                                                              *
% c     *     transform [Dt] from a form relating the unrotated stress *
% c     *     rate and the unrotated rate of deformation tensor to     *
% c     *     one relating the cauchy stress rate and the rate of      *
% c     *     (spatial) deformation strain for a block of elements.    *
% c     *                                                              *
% c     ****************************************************************
% c
% c

one = 1.d0; two = 2.d0; half = one/two;
% c
% c
% c             [cep] (mxvl x 6 x 6) relates increments
% c             of unrotated cauchy stress to increments
% c             of the unrotated deformation. transform [cep]
% c             so it relates increments of cauchy stress to
% c             increments of the deformation, both on the
% c             spatial coordinates.
% c
% c             [cep*] = [qn1] * [cep] * trans([qn1])
% c
% c             [qn1] is a rotation matrix constructed from the
% c             [R] obtained by polar decomposition of the deformation
% c             gradient, [F] =[R][U].
% c
% c             for UMATs the UMAT computed [cep] may already refer to
% c             the Cauchy stress. No rotation to be done.
% c
% c             For crystal plasticity model,
% c
%       do_transform = .true.
%       if( is_crys_pls ) do_transform = .true.
%       if( is_umat ) then
%         if( umat_stress_type .eq. 1 ) do_transform = .false.
%       end if
%       if( .not. do_transform ) go to 1000
% c
% c
% c             perform multiplication of [tc] = [qn1] * [cep]
% c
%       do j = 1, nstr
%          do i = 1, span
% c
%             tc(i,j,1)= (qn1(i,j,1)*cep(i,1,1)+
%      &                  qn1(i,j,2)*cep(i,2,1)+
%      &                  qn1(i,j,3)*cep(i,3,1)+
%      &                  qn1(i,j,4)*cep(i,4,1)+
%      &                  qn1(i,j,5)*cep(i,5,1)+
%      &                  qn1(i,j,6)*cep(i,6,1))
% c
%             tc(i,j,2)= (qn1(i,j,1)*cep(i,1,2)+
%      &                  qn1(i,j,2)*cep(i,2,2)+
%      &                  qn1(i,j,3)*cep(i,3,2)+
%      &                  qn1(i,j,4)*cep(i,4,2)+
%      &                  qn1(i,j,5)*cep(i,5,2)+
%      &                  qn1(i,j,6)*cep(i,6,2))
% c
%             tc(i,j,3)= (qn1(i,j,1)*cep(i,1,3)+
%      &                  qn1(i,j,2)*cep(i,2,3)+
%      &                  qn1(i,j,3)*cep(i,3,3)+
%      &                  qn1(i,j,4)*cep(i,4,3)+
%      &                  qn1(i,j,5)*cep(i,5,3)+
%      &                  qn1(i,j,6)*cep(i,6,3))
% c
%             tc(i,j,4)= (qn1(i,j,1)*cep(i,1,4)+
%      &                  qn1(i,j,2)*cep(i,2,4)+
%      &                  qn1(i,j,3)*cep(i,3,4)+
%      &                  qn1(i,j,4)*cep(i,4,4)+
%      &                  qn1(i,j,5)*cep(i,5,4)+
%      &                  qn1(i,j,6)*cep(i,6,4))
% c
%             tc(i,j,5)= (qn1(i,j,1)*cep(i,1,5)+
%      &                  qn1(i,j,2)*cep(i,2,5)+
%      &                  qn1(i,j,3)*cep(i,3,5)+
%      &                  qn1(i,j,4)*cep(i,4,5)+
%      &                  qn1(i,j,5)*cep(i,5,5)+
%      &                  qn1(i,j,6)*cep(i,6,5))
% c
%             tc(i,j,6)= (qn1(i,j,1)*cep(i,1,6)+
%      &                  qn1(i,j,2)*cep(i,2,6)+
%      &                  qn1(i,j,3)*cep(i,3,6)+
%      &                  qn1(i,j,4)*cep(i,4,6)+
%      &                  qn1(i,j,5)*cep(i,5,6)+
%      &                  qn1(i,j,6)*cep(i,6,6))
% c
%          end do
%       end do
% c
% c
% c                       perform multiplication of
% c                       [cep*] =  [tc] * transpose([qn1])
% c
%       do j = 1, nstr
%          do i = 1, span
% c
%             cep(i,j,1)= tc(i,j,1)*qn1(i,1,1)+
%      &                  tc(i,j,2)*qn1(i,1,2)+
%      &                  tc(i,j,3)*qn1(i,1,3)+
%      &                  tc(i,j,4)*qn1(i,1,4)+
%      &                  tc(i,j,5)*qn1(i,1,5)+
%      &                  tc(i,j,6)*qn1(i,1,6)
% c
%             cep(i,j,2)= tc(i,j,1)*qn1(i,2,1)+
%      &                  tc(i,j,2)*qn1(i,2,2)+
%      &                  tc(i,j,3)*qn1(i,2,3)+
%      &                  tc(i,j,4)*qn1(i,2,4)+
%      &                  tc(i,j,5)*qn1(i,2,5)+
%      &                  tc(i,j,6)*qn1(i,2,6)
% c
%             cep(i,j,3)= tc(i,j,1)*qn1(i,3,1)+
%      &                  tc(i,j,2)*qn1(i,3,2)+
%      &                  tc(i,j,3)*qn1(i,3,3)+
%      &                  tc(i,j,4)*qn1(i,3,4)+
%      &                  tc(i,j,5)*qn1(i,3,5)+
%      &                  tc(i,j,6)*qn1(i,3,6)
% c
%             cep(i,j,4)= tc(i,j,1)*qn1(i,4,1)+
%      &                  tc(i,j,2)*qn1(i,4,2)+
%      &                  tc(i,j,3)*qn1(i,4,3)+
%      &                  tc(i,j,4)*qn1(i,4,4)+
%      &                  tc(i,j,5)*qn1(i,4,5)+
%      &                  tc(i,j,6)*qn1(i,4,6)
% c
%             cep(i,j,5)= tc(i,j,1)*qn1(i,5,1)+
%      &                  tc(i,j,2)*qn1(i,5,2)+
%      &                  tc(i,j,3)*qn1(i,5,3)+
%      &                  tc(i,j,4)*qn1(i,5,4)+
%      &                  tc(i,j,5)*qn1(i,5,5)+
%      &                  tc(i,j,6)*qn1(i,5,6)
% c
%             cep(i,j,6)= tc(i,j,1)*qn1(i,6,1)+
%      &                  tc(i,j,2)*qn1(i,6,2)+
%      &                  tc(i,j,3)*qn1(i,6,3)+
%      &                  tc(i,j,4)*qn1(i,6,4)+
%      &                  tc(i,j,5)*qn1(i,6,5)+
%      &                  tc(i,j,6)*qn1(i,6,6)
% c
%          end do
%       end do
cep = qn1*cep*qn1';
% c
% c            subtract the [Q-bar] matrix from the transformed
% c            [cep]. this is the "initial stress" at the material
% c            point level. this remains an option indicated by qbar.
% c            note: we must multiply in the gauss weight factor and
% c            gauss point det[J] for the subtracted terms. the [cep]
% c            passed in had these factors included by the cnst...
% c            routines. the [Q-bar] 6x6 comes from the tensor
% c            expression -2 (de.De):s, where, s is the stress tensor,
% c            de is the rate of deformation tensor and De is the virtual
% c            rate of deformation tensor. this expression in matrix form
% c            is: - trans([B]) * [Q-bar] * [B]. this modification of [cep]
% c            is essential for convergence of nearly homogeneous
% c            deformation problems.
% c
%  1000 continue
      if ( qbar ) %then
%         do i = 1, span
         wf    = dj  * w;
         halfw = half * wf;
         cep( 1,1) = cep( 1,1) - two * cs( 1) * wf;
         cep( 2,2) = cep( 2,2) - two * cs( 2) * wf;
         cep( 3,3) = cep( 3,3) - two * cs( 3) * wf;
         cep( 4,1) = cep( 4,1) - cs( 4) * wf;
         cep( 6,1) = cep( 6,1) - cs( 6) * wf;
         cep( 4,2) = cep( 4,2) - cs( 4) * wf;
         cep( 5,2) = cep( 5,2) - cs( 5) * wf;
         cep( 5,3) = cep( 5,3) - cs( 5) * wf;
         cep( 6,3) = cep( 6,3) - cs( 6) * wf;
         cep( 4,4) = cep( 4,4) - halfw * ( cs( 1)+cs( 2) );
         cep( 5,5) = cep( 5,5) - halfw * ( cs( 2)+cs( 3) );
         cep( 6,6) = cep( 6,6) - halfw * ( cs( 1)+cs( 3) );
         cep( 5,4) = cep( 5,4) - halfw * cs( 6);
         cep( 6,4) = cep( 6,4) - halfw * cs( 5);
         cep( 6,5) = cep( 6,5) - halfw * cs( 4);
         cep( 1,4) = cep( 4,1);
         cep( 1,6) = cep( 6,1);
         cep( 2,4) = cep( 4,2);
         cep( 2,5) = cep( 5,2);
         cep( 3,5) = cep( 5,3);
         cep( 3,6) = cep( 6,3);
         cep( 4,5) = cep( 5,4);
         cep( 4,6) = cep( 6,4);
         cep( 5,6) = cep( 6,5);
%         end do
      end %if
% c
%       return
      end