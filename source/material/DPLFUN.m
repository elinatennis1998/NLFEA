function dplfun = DPLFUN(X, NPOINT, XFX)
% C
%       INTEGER NPOINT, I
%       DOUBLE PRECISION X, XFX(2,NPOINT), R0
%       DATA R0 / 0.0D0 /
% C***********************************************************************
% C DERIVATIVE OF THE PIECEWISE LINEAR FUNCTION 'PLFUN' DEFINED BY A SET
% C OF NPOINT PAIRS {X,F(X)} STORED IN THE MATRIX XFX (DIM. 2*NPOINT).
% C***********************************************************************
R0 = 0.d0;
      for I=1:NPOINT 
        if (X>=XFX(1,I)) 
%           GOTO 100
        else
          if (I==1) 
% C           -- x < x1   --> f(x)=f(x1) --> df(x)/dx=0 --- 
            dplfun=R0;
%             GOTO 999
            break
          else
% C           -- x(i-1) <= x < x(i) ---
            dplfun=(XFX(2,I)-XFX(2,I-1))/ ...
                   (XFX(1,I)-XFX(1,I-1));
%             GOTO 999
            break
          end %IF
        end %IF
      end
% C     ---- x >= x(npoint) --> f(x) = f(x(npoint)) --> df/dx=0 ---
      if I==NPOINT
      dplfun=R0;
      end
%  999  CONTINUE
%       RETURN
%       END
