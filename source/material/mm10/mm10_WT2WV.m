% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_WT2WV                                                                 *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 3/22/12 mcm                                      *
% c *                                                                          *
% c *         Skew   tensor to skew   vector                                   *
% c *                                                                          *
% c ****************************************************************************
% c
function WV = mm10_WT2WV(WT)
%             implicit none
%             double precision, dimension(3,3), intent(in) :: WT
%             double precision, dimension(3), intent(out) :: WV
%
            WV(1) = WT(2,3);
            WV(2) = WT(1,3);
            WV(3) = WT(1,2);
%            return
end