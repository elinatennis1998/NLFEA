% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_rotation_matrix                                                       *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 3/22/12 mcm                                      *
% c *                                                                          *
% c *         Given euler angles, an angle convention, and an angle type       *
% c *         send back the correct rotation matrix.                           *
% c *                                                                          *
% c ****************************************************************************
% c
function R = ...
    mm10_rotation_matrix(angles, aconv, atype)
%             implicit none
%             double precision, dimension(3), intent(in) :: angles
%             character, intent(in) :: aconv*5
%             character, intent(in) :: atype*7
%             integer, intent(in) :: out
% c
%             double precision, dimension(3,3), intent(out) :: R
% c
%             double precision :: a, b, c, psi, theta, phi, pi
% c
            pi = 2*acos(0.0);
% c
            a = angles(1);
            b = angles(2);
            c = angles(3);
% 
            if (strcmp(atype, 'degrees') == 1)
                  a = a*pi/180;
                  b = b*pi/180;
                  c = c*pi/180;
            elseif (strcmp(atype, 'radians') == 1)
            else
%                   write (out,9000) % commented out since it's a write function
            end
%             
            if (strcmp(aconv, 'kocks') == 1)
                  psi = a;
                  theta = b;
                  phi = c;
            elseif (strcmp(aconv, 'bunge') == 1)
                  psi = a - pi/2;
                  theta = b;
                  phi = pi/2 - c;
            elseif (strcmp(aconv, 'roe') == 1)
                  psi = a;
                  theta = b;
                  phi = 3*pi/2-c;
            else
%                   write (out,9001) % commented out since it's a write function
            end
% 
% c           Pretty sure the listing in Kocks, Tome, and Wenk is wrong
% 
            R(1,1) = -sin(psi)*sin(phi)-cos(psi)*cos(phi)*cos(theta);
            R(1,2) = cos(psi)*sin(phi)-sin(psi)*cos(phi)*cos(theta);
            R(1,3) = cos(phi)*sin(theta);
            R(2,1) = sin(psi)*cos(phi)-cos(psi)*sin(phi)*cos(theta);
            R(2,2) = -cos(psi)*cos(phi)-sin(psi)*sin(phi)*cos(theta);
            R(2,3) = sin(phi)*sin(theta);
            R(3,1) = cos(psi)*sin(theta);
            R(3,2) = sin(psi)*sin(theta);
            R(3,3) = cos(theta);
% 
%             return
%  9000 format(/'Danger: Unknown angle type passed to rotation_matrix'/)
%  9001 format(/'Danger: Unknown angle convention passed to',
%      &        ' rotation_matrix'/)
end