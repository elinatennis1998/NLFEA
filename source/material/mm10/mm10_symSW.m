% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_symSW                                                            *    *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 12/11/12 mcm                                     *
% c *                                                                          *
% c *   Take the symmetric part of S*W for some stress and skew tensor         *
% c *                                                                          *
% c ****************************************************************************
% c
function SW = mm10_symSW(S,W)
%             implicit none
%             double precision, dimension(6), intent(in) :: S
%             double precision, dimension(3), intent(in) :: W
%             double precision, dimension(6), intent(out) :: SW
% c
            SW = zeros(6,size(W,2));
            SW(1,:) = S(4)*W(3,:) - S(6)*W(2,:);
            SW(2,:) = S(4)*W(3,:) - S(5)*W(1,:);
            SW(3,:) = S(6)*W(2,:) + S(5)*W(1,:);
            SW(4,:) = 0.5*(W(3,:)*(S(1)-S(2)) + W(1,:)*S(6) - W(2,:)*S(5));
            SW(5,:) = 0.5*(W(1,:)*(S(2)-S(3)) + W(2,:)*S(4) + W(3,:)*S(6));
            SW(6,:) = 0.5*(W(2,:)*(S(1)-S(3)) + W(1,:)*S(4) - W(3,:)*S(5));
% 
%             return
% 
end