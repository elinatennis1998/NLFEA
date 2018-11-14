% c
% c ****************************************************************************
% c *                                                                          *
% c *    mm10_RT2RVE                                                                *
% c *                                                                          *
% c *         written by : mcm                                                 *
% c *         last modified : 3/22/12 mcm                                      *
% c *                                                                          *
% c *         Takes a 3x3 rotation tensor and returns it in a 6x6 form         *
% c *         suitable for rotating Voigt-type strain vectors                  *
% c *                                                                          *
% c ****************************************************************************
% c
function RV = mm10_RT2RVE(RT)
%             implicit none
%             double precision, dimension(3,3), intent(in) :: RT
%             double precision, dimension(6,6), intent(out) :: RV
%        
            RV(1,1)=RT(1,1)^2;
            RV(1,2)=RT(1,2)^2;
            RV(1,3)=RT(1,3)^2;
            RV(1,4)=2*RT(1,1)*RT(1,2);
            RV(1,5)=2*RT(1,3)*RT(1,2);
            RV(1,6)=2*RT(1,1)*RT(1,3);
            RV(2,1)=RT(2,1)^2;
            RV(2,2)=RT(2,2)^2;
            RV(2,3)=RT(2,3)^2;
            RV(2,4)=2*RT(2,1)*RT(2,2);
            RV(2,5)=2*RT(2,3)*RT(2,2);
            RV(2,6)=2*RT(2,1)*RT(2,3);
            RV(3,1)=RT(3,1)^2;
            RV(3,2)=RT(3,2)^2;
            RV(3,3)=RT(3,3)^2;
            RV(3,4)=2*RT(3,1)*RT(3,2);
            RV(3,5)=2*RT(3,3)*RT(3,2);
            RV(3,6)=2*RT(3,1)*RT(3,3);
            RV(4,1)=RT(1,1)*RT(2,1);
            RV(4,2)=RT(1,2)*RT(2,2);
            RV(4,3)=RT(1,3)*RT(2,3);
            RV(4,4)=RT(1,1)*RT(2,2)+RT(2,1)*RT(1,2);
            RV(4,5)=RT(1,2)*RT(2,3)+RT(1,3)*RT(2,2);
            RV(4,6)=RT(1,1)*RT(2,3)+RT(1,3)*RT(2,1);
            RV(5,1)=RT(2,1)*RT(3,1);
            RV(5,2)=RT(3,2)*RT(2,2);
            RV(5,3)=RT(2,3)*RT(3,3);
            RV(5,4)=RT(2,1)*RT(3,2)+RT(2,2)*RT(3,1);
            RV(5,5)=RT(2,2)*RT(3,3)+RT(3,2)*RT(2,3);
            RV(5,6)=RT(2,1)*RT(3,3)+RT(2,3)*RT(3,1);
            RV(6,1)=RT(1,1)*RT(3,1);
            RV(6,2)=RT(1,2)*RT(3,2);
            RV(6,3)=RT(1,3)*RT(3,3);
            RV(6,4)=RT(1,1)*RT(3,2)+RT(1,2)*RT(3,1);
            RV(6,5)=RT(1,2)*RT(3,3)+RT(1,3)*RT(3,2);
            RV(6,6)=RT(1,1)*RT(3,3)+RT(3,1)*RT(1,3);
% 
%             return
end