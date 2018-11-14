% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_IW                                                               *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 12/11/12 mcm                                     *
% c *                                                                          *
% c *   Get Iik*Wlj - Wik*Ijl in our voigt notation                            *
% c *                                                                          *
% c ****************************************************************************
% c
function IW = mm10_IW(W)
%            implicit none
%             double precision, dimension(3), intent(in) :: W
%             double precision, dimension(6,6), intent(out) :: IW
% c
            IW = 0.0;
            IW(1,4) =  2.0*W(3);
            IW(1,6) = -2.0*W(2);
            IW(2,4) =  2.0*W(3);
            IW(2,5) = -2.0*W(1);
            IW(3,5) =  2.0*W(1);
            IW(3,6) =  2.0*W(2);
            IW(4,1) =  W(3);
            IW(4,2) = -W(3);
            IW(4,5) = -W(2);
            IW(4,6) =  W(1);
            IW(5,2) =  W(1);
            IW(5,3) = -W(1);
            IW(5,4) =  W(2);
            IW(5,6) =  W(3);
            IW(6,1) =  W(2);
            IW(6,3) = -W(2);
            IW(6,4) =  W(1);
            IW(6,5) = -W(3);
% 
%             return
end