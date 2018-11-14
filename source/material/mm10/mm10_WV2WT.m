% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_WV2WT                                                                 *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 3/22/12 mcm                                      *
% c *                                                                          *
% c *         Skew vector to skew tensor                                       *
% c *                                                                          *
% c ****************************************************************************
% c
function WT = mm10_WV2WT(WV)
%             implicit none
%             double precision, dimension(3,3), intent(out) :: WT
%             double precision, dimension(3), intent(in) :: WV
%
WT = zeros(3,3);
           WT(1,1) = 0.0;
            WT(1,2) = WV(3);
            WT(1,3) = WV(2);
            WT(2,1) = -WV(3);
            WT(2,2) = 0;
            WT(2,3) = WV(1);
            WT(3,1) = -WV(2);
            WT(3,2) = -WV(1);
            WT(3,3) = 0;
%            
end