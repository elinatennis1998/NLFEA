% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_expmw3x3                                                              *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 3/22/12 mcm                                      *
% c *                                                                          *
% c *         Calculates exp(W) where W is a 3x3 skew matrix.                  *
% c *         Returns full 3x3 because the result                              *
% c *         is only orthogonal (not skew or symmetric)                       *
% c *                                                                          *
% c ****************************************************************************
% c
% 09/10/2014: fixed a bug in the Rodrigues' formula
% Also, this function is 10 times faster than expm for rotation matrices.
function A = mm10_expmw3x3(W)
%             implicit none
%             double precision, dimension(3,3), intent(in) :: W
%             double precision, dimension(3,3), intent(out) :: A
%             double precision :: alpha
%             integer :: i
% c
% c           Compute alpha := norm(vector(W))
            alpha = sqrt(W(2,3)^2+W(1,3)^2+W(1,2)^2);
% c
% c           Algorithm will fail with alpha = 0 (=> W=0)
% c           Correct result is expm(0) = I, so program that in
            if (alpha < 1.0E-16)
                  A = zeros(3,3);
            else
                  A = W;
                  A = (1-cos(alpha))/(alpha^2)*W*W + sin(alpha)/alpha*A;
            end
% 
% c           Add the identity
                  A = A + eye(3);
% 
%             return
end