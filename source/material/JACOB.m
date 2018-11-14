function [D,V] = JACOB(A,N)
%       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
%       PARAMETER (MJITER=50,NMAX=100)
%       DIMENSION
%      1    A(N,N)     ,D(N)      ,V(N,N)
%       DIMENSION
%      1    B(NMAX)    ,Z(NMAX)
%       DATA R0   ,RP2  ,RP5  ,R1   ,R100   /
%      1     0.0D0,0.2D0,0.5D0,1.0D0,100.0D0/
%       DATA TOLER  /
%      1     1.0D-12/
% C***********************************************************************
% C JACOBI ITERATIVE PROCEDURE FOR SPECTRAL DECOMPOSITION OF A
% C N-DIMENSIONAL SYMMETRIC MATRIX
% C
% C REFERENCE: WH Press, SA Teukolsky, WT Vetting & BP Flannery. Numerical
% C            recipes in FORTRAN: The art of scientific computing. 2nd
% C            Edn., Cambridge University Press, 1992.
% C***********************************************************************
MJITER=50;
NMAX=100;
R0 = 0.0D0;
RP2 = 0.2D0;
RP5 = 0.5D0;
R1 = 1.0D0;
R100 = 100.0D0;
TOLER=1.D-12;
D = zeros(N,1);
B = zeros(NMAX,1);
Z = zeros(NMAX,1);
      if(N>NMAX)
        return
      end %IF
      V = eye(N);
%       for IP=1:N
%         for IQ=1:N
%           V(IP,IQ)=R0;
%         end
%         V(IP,IP)=R1;
%       end
      for IP=1:N
        B(IP)=A(IP,IP);
        D(IP)=B(IP);
        Z(IP)=R0;
      end
      for I=1:MJITER
        SM=R0;
        for IP=1:N-1
          for IQ=IP+1:N
            SM=SM+abs(A(IP,IQ));
          end
        end
        if(SM<TOLER), break, end
        if(I<4)
          TRESH=RP2*SM/N^2;
        else
          TRESH=R0;
        end %IF
        for IP=1:N-1
          for IQ=IP+1:N
            G=R100*abs(A(IP,IQ));
            if((I>4)&&(abs(D(IP))+G==abs(D(IP)))...
               &&(abs(D(IQ))+G==abs(D(IQ))))
              A(IP,IQ)=R0;
            elseif(abs(A(IP,IQ))>TRESH)
              H=D(IQ)-D(IP);
              if(abs(H)+G==abs(H))
                T=A(IP,IQ)/H;
              else
                THETA=RP5*H/A(IP,IQ);
                T=R1/(abs(THETA)+sqrt(R1+THETA^2));
                if(THETA<R0), T=-T; end
              end %IF
              C=R1/SQRT(R1+T^2);
              S=T*C;
              TAU=S/(R1+C);
              H=T*A(IP,IQ);
              Z(IP)=Z(IP)-H;
              Z(IQ)=Z(IQ)+H;
              D(IP)=D(IP)-H;
              D(IQ)=D(IQ)+H;
              A(IP,IQ)=R0;
              for J=1:IP-1
                G=A(J,IP);
                H=A(J,IQ);
                A(J,IP)=G-S*(H+G*TAU);
                A(J,IQ)=H+S*(G-H*TAU);
              end
              for J=IP+1:IQ-1
                G=A(IP,J);
                H=A(J,IQ);
                A(IP,J)=G-S*(H+G*TAU);
                A(J,IQ)=H+S*(G-H*TAU);
              end
              for J=IQ+1:N
                G=A(IP,J);
                H=A(IQ,J);
                A(IP,J)=G-S*(H+G*TAU);
                A(IQ,J)=H+S*(G-H*TAU);
              end
              for J=1:N
                G=V(J,IP);
                H=V(J,IQ);
                V(J,IP)=G-S*(H+G*TAU);
                V(J,IQ)=H+S*(G-H*TAU);
              end
            end %IF
          end
        end
        for IP=1:N
          B(IP)=B(IP)+Z(IP);
          D(IP)=B(IP);
          Z(IP)=R0;
        end
      end
%       CALL ERRPRT('EE0005')
%   999 CONTINUE
%       RETURN
%       END