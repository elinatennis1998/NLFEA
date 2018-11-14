function lint = IntPoint3(nel)
% !
% !**********************************************************************
% 	Implicit None
% 	Integer Nel,Lint
          if    (nel==4)
              lint=4;14;   
          elseif(nel==10)
              lint=13;
          elseif(nel==8)
              lint=64;
          elseif(nel==27)
              lint=64;
          end