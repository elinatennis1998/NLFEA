function plfun = PLFUN(X, NPOINT, XFX)
% C
%       INTEGER NPOINT, I
%       DOUBLE PRECISION X, XFX(2,*)
% C***********************************************************************
% C PIECEWISE LINEAR FUNCTION DEFINED BY A SET OF NPOINT PAIRS
% C {X,F(X)} STORED IN THE MATRIX XFX (DIM. 2*NPOINT).
% C***********************************************************************
      for I=1:NPOINT 
        if (X>=XFX(1,I))
%           GOTO 100
        else  
          if (I==1)
% C           -- x < x1 --> f(x)=f(x1) --- 
            plfun=XFX(2,1);
%             GOTO 999
            break
          else
% C           -- x(i-1) <= x < x(i) ---
            plfun=XFX(2,I-1)+(X-XFX(1,I-1))* ...
                       (XFX(2,I)-XFX(2,I-1))/ ...
                       (XFX(1,I)-XFX(1,I-1));
%             GOTO 999
            break
          end %IF
        end %IF
      end
% C     ----  x >= x(npoint) --> f(x) = f(x(npoint))  ---
      if I==NPOINT
      plfun=XFX(2,NPOINT);
      end
%  999  CONTINUE
%       RETURN
%       END
