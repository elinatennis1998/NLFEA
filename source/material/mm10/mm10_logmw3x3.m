% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_logmw3x3                                                         *
% c *                                                                          *
% c *         written by : tjt                                                 *
% c *         last modified : 9/10/14 tjt                                      *
% c *                                                                          *
% c *         Calculates log(W) where W is a 3x3 rotation matrix.              *
% c *         Returns full 3x3; uses Wikipedia formula                         *
% c *                                                                          *
% c ****************************************************************************
% c
% 30 times faster than logm(W)
function A = mm10_logmw3x3(W)
%             implicit none
%             double precision, dimension(3,3), intent(in) :: W
%             double precision, dimension(3,3), intent(out) :: A
%             double precision :: alpha
%             integer :: i
% c
% c           Compute alpha := norm(vector(W))
            theta = acos((trace(W)-1)/2);
% c
% c           Algorithm will fail with alpha = 0 (=> W=0)
% c           Correct result is expm(0) = I, so program that in
            if (abs(theta) < 1.0E-16)
                  A = zeros(3,3);
            elseif abs(abs(theta) - pi) < 1e-14
                  error('not sure what to do for logm(pi)')
            else
                  A = theta/(2*sin(theta))*(W-W');
            end
% 
%             return
end