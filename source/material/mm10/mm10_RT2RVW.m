% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_RT2RVW                                                                *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 3/22/12 mcm                                      *
% c *                                                                          *
% c *         Takes a 3x3 rotation tensor and returns it in a 3x3 form         *
% c *         suitable for rotating my 3x1 skew vectors                        *
% c *                                                                          *
% c ****************************************************************************
% c
function RV = mm10_RT2RVW(RT)
%             implicit none
%             double precision, dimension(3,3), intent(in) :: RT
%             double precision, dimension(3,3), intent(out) :: RV
% 
            RV(1,1)=RT(2,2)*RT(3,3)-RT(2,3)*RT(3,2);
            RV(1,2)=RT(2,1)*RT(3,3)-RT(2,3)*RT(3,1);
            RV(1,3)=RT(2,1)*RT(3,2)-RT(2,2)*RT(3,1);
            RV(2,1)=RT(1,2)*RT(3,3)-RT(1,3)*RT(3,2);
            RV(2,2)=RT(1,1)*RT(3,3)-RT(1,3)*RT(3,1);
            RV(2,3)=RT(1,1)*RT(3,2)-RT(1,2)*RT(3,1);
            RV(3,1)=RT(1,2)*RT(2,3)-RT(1,3)*RT(2,2);
            RV(3,2)=RT(1,1)*RT(2,3)-RT(1,3)*RT(2,1);
            RV(3,3)=RT(1,1)*RT(2,2)-RT(1,2)*RT(2,1);
%             
%             return
end