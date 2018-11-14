function lamda = evcmp1(k)
% c     ****************************************************************
% c     *                                                              *
% c     *                      subroutine evcmp1                       *
% c     *                                                              *
% c     *                       written by : bh                        *
% c     *                                                              *
% c     *                   last modified : 07/17/90                   *
% c     *                                 : 02/09/94                   *
% c     *                                                              *
% c     *     this subroutine computes the eigenvalues of a 3x3        *
% c     *     symmetric positive definite matrix in stress vector      *
% c     *     form for a block of solid elements.                      *
% c     *                                                              *
% c     ****************************************************************
% c
% c
% c
zero = 0.d0; one = 1.d0; two = 2.d0; jactol = 1.d-8; four = two*two; 
ten = 10.d0; ten_thouth = 0.0001d0;
m = zeros(3,1); lamda = m;
maxswp = 15;
% c
% c
% c              initialize lamda, m, sweep parameters.
% c
% c
      swpnum = 0;
% c
%       do 10 bel= 1,span
% c
         m(1)= one;
         m(2)= one;
         m(3)= one;
         lamda(1) = k(1);
         lamda(2) = k(2);
         lamda(3) = k(3);
% c
% c              scale [k] and ^0m^2 to avoid problems with exponential
% c              overflow and underflow.
% c
% c              find the max and min terms on the diagonal of [k] & ^0m^2
% c
         kj = k(1);
         kj = min( k(2),kj );
         kj = min( k(3),kj );
         ki = k(1);
         ki = max( k(2),ki );
         ki = max( k(3),ki );
         mj = one;
         mi = one;
% c
% c              compute the scale factor and do the scaling
% c
% #sgl         iexp = int( ( log10(kj)+log10(ki)+
         iexp = round( ( log10(kj)+log10(ki)+ ...
                            log10(mj)+log10(mi) ) / four );
         scale = one / ( ten ^ iexp );
         m(1) = m(1) * scale;
         m(2) = m(2) * scale;
         m(3) = m(3) * scale;
         k(1) = k(1) * scale;
         k(4) = k(4) * scale;
         k(2) = k(2) * scale;
         k(6) = k(6) * scale;
         k(5) = k(5) * scale;
         k(3) = k(3) * scale;
% c
%  10   continue
% c
% c              begin a new sweep
% c
      cvgtst = 0;
      
      while swpnum < maxswp && ~cvgtst
      swpnum = swpnum + 1;
      thold = ten_thouth ^ swpnum;
      sqtol = jactol * jactol;
      if( thold < sqtol ) 
          thold = sqtol;
      end
% c
% c              enter sweep loop -- work on lower triangle only
% c                                          ( i > j )
% c
% c              rows are done from top to bottom
% c              columns are done from left to right.
% c
% c
% c           ***************************************
% c           *                                     *
% c           *           row 2 and column 1.       *
% c           *                                     *
% c           ***************************************
% c
% c
%       do 20 bel= 1,span
% c
% c                       check if term is within threshold
% c
         ratiok= (k(4)*k(4))/(k(2)*k(1));
         if ( ratiok >= thold ) %then
% c
% c                      compute the rotatiom matrix:  an identity
% c                      matrix with alpha at position (2,1) and
% c                      gamma at position (1,2).
% c
            kbari= -m(2)*k(4);
            kbarj= -m(1)*k(4);
            kbar = k(2)*m(1)-k(1)*m(2);
            rad  = (kbar*kbar/four)+kbari*kbarj;
% c
%             xsign= one;
%             x= kbar/two+sign(xsign,kbar)* ...
%                     sqrt(rad);
            x= kbar/two+sign(kbar)* ...
                    sqrt(rad);

            if( (abs(x)<jactol*abs(kbarj)) || ...
                (abs(x)<jactol*abs(kbari))    ) %then
               alpha= zero;
               gamma= -k(4)/k(2);
            else
               alpha= kbarj/x;
               gamma= -kbari/x;
            end %if
% c
% c                       perform the rotation.
% c
% c
% c                       row 3, column 2
% c                       row 3, column 1
% c
            ki = k(5);
            kj = k(6);
            k(5)= ki+gamma*kj;
            k(6)= kj+alpha*ki;
% c
% c                       term (2,1) and diagonal terms (2,2) and (1,1).
% c
% c
            kj = k(1);
            mj = m(1);
            ki = k(2);
            mi = m(2);
            k(1)= kj+alpha*alpha*ki+ ...
                      two*alpha*k(4);
            m(1)= mj+alpha*alpha*mi;
            k(2)= ki+gamma*gamma*kj+ ...
                      two*gamma*k(4);
            m(2)= mi+gamma*gamma*mj;
            k(4)= zero;
% c
         end %if
% c
%  20   continue
% c
% c
% c           ***************************************
% c           *                                     *
% c           *           row 3 and column 1.       *
% c           *                                     *
% c           ***************************************
% c
% c
%       do 25 bel= 1,span
% c
% c                       check if term is within threshold
% c
         ratiok = (k(6)*k(6))/(k(3)*k(1));
         if ( ratiok>=thold ) %then
% c
% c                       compute the rotatiom matrix:  an identity
% c                       matrix with alpha at position (3,1) and
% c                       gamma at position (1,3).
% c
            kbari= -m(3)*k(6);
            kbarj= -m(1)*k(6);
            kbar = k(3)*m(1)-k(1)*m(3);
            rad  = (kbar*kbar/four)+kbari*kbarj;
% c
%             xsign= one;
            x= kbar/two+sign(kbar)* ...
                    sqrt(rad);
            if( (abs(x)<jactol*abs(kbarj)) || ...
                (abs(x)<jactol*abs(kbari))    ) %then
               alpha = zero;
               gamma = -k(6)/k(3);
            else
               alpha =  kbarj/x;
               gamma = -kbari/x;
            end %if
% c
% c                       perform the rotation.
% c
% c
% c                       row 3, column 2
% c                       row 2, column 1
% c
            ki  = k(5);
            kj  = k(4);
            k(5) = ki+gamma*kj;
            k(4) = kj+alpha*ki;
% c
% c                       term (3,1) and diagonal terms (3,3) and (1,1).
% c
            kj = k(1);
            mj = m(1);
            ki = k(3);
            mi = m(3);
            k(1)= kj+alpha*alpha*ki+ ...
                      two*alpha*k(6);
            m(1)= mj+alpha*alpha*mi;
            k(3)= ki+gamma*gamma*kj+ ...
                      two*gamma*k(6);
            m(3)= mi+gamma*gamma*mj;
            k(6)= zero;
% c
         end %if
% c
%  25   continue
% c
% c
% c           ***************************************
% c           *                                     *
% c           *           row 3 and column 2.       *
% c           *                                     *
% c           ***************************************
% c
% c
%       do 30 bel= 1,span
% c
% c                       check if term is within threshold
% c
         ratiok = (k(5)*k(5))/(k(3)*k(2));
         if ( ratiok>=thold ) %then
% c
% c                       compute the rotatiom matrix:  an identity
% c                       matrix with alpha at position (3,2) and
% c                       gamma at position (2,3).
% c
            kbari= -m(3)*k(5);
            kbarj= -m(2)*k(5);
            kbar =  k(3)*m(2)-k(2)*m(3);
            rad  = (kbar*kbar/four)+kbari*kbarj;
% c
            xsign= one;
            x= kbar/two+sign(kbar)* ...
                    sqrt(rad);
            if( (abs(x)<jactol*abs(kbarj)) || ...
                (abs(x)<jactol*abs(kbari))    ) %then
               alpha= zero;
               gamma= -k(5)/k(3);
            else
               alpha=  kbarj/x;
               gamma= -kbari/x;
            end %if
% c
% c                       perform the rotation.
% c
% c
% c                       row 3, column 1
% c                       row 2, column 1
% c
            ki = k(6);
            kj = k(4);
            k(6)= ki+gamma*kj;
            k(4)= kj+alpha*ki;
% c
% c                       term (3,2) and diagonal terms (3,3) and (2,2).
% c
            kj = k(2);
            mj = m(2);
            ki = k(3);
            mi = m(3);
            k(2)= kj+alpha*alpha*ki+ ...
                      two*alpha*k(5);
            m(2)= mj+alpha*alpha*mi;
            k(3)= ki+gamma*gamma*kj+ ...
                      two*gamma*k(5);
            m(3)= mi+gamma*gamma*mj;
            k(5)= zero;
% c
         end %if
% c
%  30   continue
% c
% c              end sweep
% c
% c              update eigenvalue vector -- lamda
% c

%       do 35 bel= 1,span
% c
         lamda(1)= k(1)/m(1);
         lamda(2)= k(2)/m(2);
         lamda(3)= k(3)/m(3);
% c
%  35   continue
% c
% c              check off-diagonal elements for convergence
% c
      cvgtst= 1;%.true.
% c
%       do 40 bel= 1,span
% c
         errork= k(4)*k(4)/(k(2)*k(1));
         if (errork>sqtol) 
             cvgtst= 0;%.false.
         end
% c
         errork= k(6)*k(6)/(k(3)*k(1));
         if (errork>sqtol) 
             cvgtst= 0;%.false.
         end
% c
         errork= k(5)*k(5)/(k(3)*k(2));
         if (errork>sqtol) 
             cvgtst= 0;%.false.
         end
% c
%  40   continue
% c
%       if( cvgtst ) go to 45
%       if( swpnum .lt. maxswp ) go to 15
      end
          
% c
% c                       reorder the eigenvalues
% c
%  45   do 50 bel= 1,span
% c
         if(lamda(2) < lamda(1)) %then
            swap= lamda(1);
            lamda(1)= lamda(2);
            lamda(2)= swap;
         end %if
% c
         if(lamda(3) < lamda(1)) %then
            swap= lamda(1);
            lamda(1)= lamda(3);
            lamda(3)= swap;
         end %if
% c
         if(lamda(3) < lamda(2)) %then
            swap= lamda(2);
            lamda(2)= lamda(3);
            lamda(3)= swap;
         end %if
% c
%  50   continue
% c
%       return
      end