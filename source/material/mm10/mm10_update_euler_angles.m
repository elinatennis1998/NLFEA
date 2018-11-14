% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_update_euler_angles          *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 3/28/12                     *
% c     *                                                              *
% c     *    update euler angles to the new rotation                   *
% c     *                                                              *
% c     ****************************************************************
% c
function [props,np1,n] = mm10_update_euler_angles(props,np1,n)
%             use mm10_defs
%             implicit none
%             type (crystal_props) :: props
%             type (crystal_state) :: n, np1
% c
%             double precision, dimension(3,3) :: full_rot
%             double precision :: psiK, phiK, thetaK, psi, phi, theta,
%      &            pi, pps, pms, tol
%             double precision, external :: atan2_zp
%             parameter(tol=1.0E-16)
% c
            pi = 2.0*acos(0.0);
% c
% c                 Note: This subroutine needs a major fix if I'm
% c                 ever going to support anything other than degrees+
% c                 Kocks convention
            full_rot = props.g*(np1.Rp*transpose(np1.R));
% c
            psiK = atan2_zp(full_rot(3,2),full_rot(3,1));
            phiK = atan2_zp(full_rot(2,3),full_rot(1,3));
            thetaK = acos(max(min(full_rot(3,3),1),-1)); %ensure real output
%             
            if (props.angle_convention == 1)
                  psi = psiK;
                  phi = phiK;
                  theta = thetaK;
            else
%                   write (*,*) "Angle convention not implemented."
%                   call die_gracefully
            end
% c
            if (props.angle_type == 1)
                  np1.euler_angles(1) = 180.0/pi*psi;
                  np1.euler_angles(2) = 180.0/pi*theta;
                  np1.euler_angles(3) = 180.0/pi*phi;
            elseif (props.angle_type == 2)
                  np1.euler_angles(1) = psi;
                  np1.euler_angles(2) = theta;
                  np1.euler_angles(3) = phi;
            else
%                   write (*,*) "Unrecognized angle convention."
%                   call die_gracefully
            end
% c
%             return
end