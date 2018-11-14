function QMATX = ATMDFB(   AMATX            ,STRES      )
      MADIM=4;
      QMATX = zeros(MADIM,MADIM);
      
R0 = 0;  
RP5 = 1/2;
% C***********************************************************************
% C COMPUTE THE ADDITIONAL TANGENT MODULUS "q" REQUIRED BY F-BAR
% C ELEMENTS:
% C                       1                 2
% C               q  :=  --- a:(I (x) I) - --- [sigma] (x) I
% C                       3                 3
% C FOR AXISYMMETRIC CASE, AND
% C                       1                 1
% C               q  :=  --- a:(I (x) I) - --- [sigma] (x) I
% C                       2                 2
% C FOR PLANE STRAIN.
% C
% C REFERENCE: Expressions (15.11) and (15.22)
% C***********************************************************************
% C
%       IF(NTYPE.EQ.2)THEN
% C Plane strain
        A=RP5;
        B=-RP5;
        QMATX(1,1)=A*(AMATX(1,1)+AMATX(1,4))+B*STRES(1);
        QMATX(2,1)=A*(AMATX(2,1)+AMATX(2,4))+B*STRES(3);
        QMATX(3,1)=A*(AMATX(3,1)+AMATX(3,4))+B*STRES(3);
        QMATX(4,1)=A*(AMATX(4,1)+AMATX(4,4))+B*STRES(2);
        QMATX(1,2)=R0;
        QMATX(2,2)=R0;
        QMATX(3,2)=R0;
        QMATX(4,2)=R0;
        QMATX(1,3)=R0;
        QMATX(2,3)=R0;
        QMATX(3,3)=R0;
        QMATX(4,3)=R0;
        QMATX(1,4)=QMATX(1,1);
        QMATX(2,4)=QMATX(2,1);
        QMATX(3,4)=QMATX(3,1);
        QMATX(4,4)=QMATX(4,1);
