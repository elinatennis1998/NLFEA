% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_ET2EV                                                                 *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 3/22/12 mcm                                      *
% c *                                                                          *
% c *         Strain tensor to strain vector                                   *
% c *                                                                          *
% c ****************************************************************************
% c
function EV = mm10_ET2EV(ET)
%             implicit none
%             double precision, dimension(3,3), intent(in) :: ET
%             double precision, dimension(6), intent(out) :: EV
%
            EV(1) = ET(1,1);
            EV(2) = ET(2,2);
            EV(3) = ET(3,3);
            EV(4) = 2*ET(1,2);
            EV(5) = 2*ET(2,3);
            EV(6) = 2*ET(1,3);
%
end
