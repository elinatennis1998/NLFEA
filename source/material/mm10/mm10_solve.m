% c
% c
% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_solve                        *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 11/26/13                    *
% c     *                                                              *
% c     *     Solve a stress update to prescribed strain state         *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, fail, J] = mm10_solve(props, np1,...
    n, stress, tt, fail, cuts, step)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       logical :: fail
% c
%       double precision, dimension(7) :: R, x, dx, xnew
%       double precision, dimension(7,7) :: J
%       double precision :: nR, inR, atol, rtol, uB, alpha, ls1, ls2,
%      &      nlsx, nRs, c, red
%       integer :: iter, miter, info, ls, mls
%       integer, dimension(7) :: ipiv
%       logical :: debug
[~,max_uhard] = maxparamsCP;
atol1 = 1.0e-6;1;
rtol1 = 1.0e-12;
atol2 = 1e-7;
rtol2 = 1.0e-12;
atol = 1e-9;
if props.miter > 0
    miter = props.miter;
else
    miter = 30;15;250;
end
mmin = 1;
debug = 0;1;
c = 1.0e-4;
red = 0.5;
mls = 4;3;6;10;

% Numerical solving technique; 0=Warp3d, 1=fsolve trust region, 2=fsolve
% finite difference method.
usefsolve = props.usefsolve;%0;2;1;3;
Jtol1 = 1e-4;1;1e-2;
Jtol2 = 1e-4;1e-2;1;
% Note: for finite difference method, only need these subroutines defined:
% -mm10_slipinc
% -mm10_h
% -mm10_ed
% -mm10_dgdd

% Intermediate vectors and arrays to help with pre-computing
vec1 = zeros(max_uhard,1);
vec2 = zeros(max_uhard,1);
arr1 = zeros(max_uhard,max_uhard);
arr2 = zeros(max_uhard,max_uhard);

      x(1:6) = stress;
      x(7:props.num_hard+6) = tt;
% c
       if (debug) 
           disp('Entering solution routine')
       end %if
      
% Prediction of yield stress to initialize the integration algorithm; helps
% for Orowan flow rule type models
  if ( props.h_type > 3 && (props.h_type ~=6) )% new models
      
%       % Module to extrapolate the hardening variables      
%       dt = np1.tinc;
%       % Predict the new hardening variable by extrapolation
%       [props, np1, n, x(1:6), x(7:props.num_hard+6), tt_h] = mm10_formR2(props, np1,...
%           n, x(1:6), x(7:props.num_hard+6));
%       tt_new = tt - tt_h;
%       rate = (tt_new - n.tau_tilde)/dt;
%       x(7:props.num_hard+6) = tt + rate*dt*step;
      
     if usefsolve == 2
          
          
      if norm(stress) == 0
          stress_typical = ones(1,6);
      else
          stress_typical = stress;
      end
      str = x(1:6);
      x(7:props.num_hard+6) = x(7:props.num_hard+6) + n.tt_rate(1:props.num_hard)'*np1.tinc;
      options = optimoptions('fsolve','Jacobian','off','Algorithm','trust-region-dogleg','TolFun',atol,'MaxIter',50,'Display','off','TypicalX',[0.1*norm(stress_typical(1:6))*ones(1,6)]');%,'Display','iter'
      [str,fval,exitflag,output,J] = fsolve(@(str) mm10_formRfsolve1(str,props,np1,n,x(7:props.num_hard+6),vec1,vec2,arr1,arr2), str, options);
      
        iter = output.iterations;
        rrr = rcond(J);
        if rrr < 1e-15
            J;
        end
          nR1 = norm(fval(1:6));
      
        if exitflag == 0 || (exitflag < 0 && nR1 > atol1)
          fail = 1; % since .true. = 1
          return;%error('stress fail')
        end %if
      x(1:6) = str;
        
         
     else
      
         
      % Initial search for current yield stress assuming hardening is fixed
      iter = 0;
      x(7:props.num_hard+6) = x(7:props.num_hard+6) + n.tt_rate(1:props.num_hard)'*np1.tinc;
     [props, np1, n, x(1:6), x(7:props.num_hard+6), vec1,vec2] = mm10_formvecs(props,np1,n,x(1:6), x(7:props.num_hard+6),...
     vec1,vec2);
      [props, np1, n, vec1, vec2, x(1:6), x(7:props.num_hard+6), R] = mm10_formR1(props, np1,...
          n, vec1, vec2, x(1:6), x(7:props.num_hard+6));
      nR1 = sqrt(dot(R(1:6),R(1:6)));
      inR1 = nR1;
      if (debug) 
            fprintf('Iter %i  norm %10.3e \n)',iter, nR1)
      end %if
      while ((nR1 > atol1) && (nR1/inR1 > rtol1))
% c           Jacobian
        if usefsolve == 3 % complex derivative solver
        [props, np1, n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6), J] = mm10_formJ11i(props, np1,...
            n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6));
        elseif usefsolve == 4 % numerically verify real finite differences
        [props, np1, n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6), J] = mm10_formJ11r(props, np1,...
            n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6));
[props, np1, n, x(1:6), x(7:props.num_hard+6), vec1, vec2, arr1, arr2] = mm10_formarrs(props, np1,...
    n, x(1:6), x(7:props.num_hard+6), vec1, vec2, arr1, arr2,1);
        [props, np1, n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6), Jcheck] = mm10_formJ11(props, np1,...
            n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6));
            if norm(J-Jcheck)/norm(J) > Jtol1
                display('finite difference check failed: norm(J-Jcheck)/norm(J) too large')
                keyboard
            else
                J = Jcheck;
            end
        else
[props, np1, n, x(1:6), x(7:props.num_hard+6), vec1, vec2, arr1, arr2] = mm10_formarrs(props, np1,...
    n, x(1:6), x(7:props.num_hard+6), vec1, vec2, arr1, arr2,1);
        [props, np1, n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6), J] = mm10_formJ11(props, np1,...
            n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6));
        end
% c           Increment
        dx = R;
        rrr = rcond(J);
        if rrr < 1e-15
            J;
            keyboard
        end
        dx = -J\dx;%[7, 1, -J, 7, ipiv, dx, 7, info] = DGESV(7, 1, -J,...
            %7, ipiv, dx, 7, info); % http://tinyurl.com/oqk27tm
% c                        % computes the solution to a real system of linear equations A * X = B, where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
% c           Line search
        alpha = 1.0;
        ls11 = 0.5*dot(R(1:6),R(1:6));
        dy = (transpose(J)*R);
        ls21 = c*dot(dx(1:6), dy(1:6));
        ls = 0;
        tr = 1;
        while tr == 1 && ls < mls
          ls = ls + 1;
          nlsx1 = ls11 + ls21*alpha;
          xnew(1:6) = x(1:6) + alpha*dx(1:6)';
          xnew(7:props.num_hard+6) = x(7:props.num_hard+6);
% c             Residual
     [props, np1, n, xnew(1:6), xnew(7:props.num_hard+6), vec1,vec2] = mm10_formvecs(props,np1,n,xnew(1:6), xnew(7:props.num_hard+6),...
     vec1,vec2);
          [props, np1, n,vec1,vec2, xnew(1:6), xnew(7:props.num_hard+6), R] = ...
              mm10_formR1(props, np1, n,vec1,vec2, xnew(1:6), xnew(7:props.num_hard+6));
          nR1 = sqrt(dot(R(1:6),R(1:6)));
          nRs1 = 0.5*dot(R(1:6),R(1:6));
% c
          if ((nRs1 <= nlsx1)|| (ls >= mls))
            x = xnew;
            tr = 0;
            break
          else
            alpha = red*alpha;
          end %if
        end %do
%         if ((jiter > mls))
%           fail = 1; % since .true. = 1
%           error('tr fail')
%         end %if
% c
        iter = iter + 1;
        if (debug) 
            fprintf('Iter %i  norm %10.3e ls %10.3f \n',iter, nR1, alpha)
     
        end %if
% c           Increment and check for failure
        if ((iter > miter) || any(isnan(x)) || any(imag(x)))
          fail = 1; % since .true. = 1
          return;%error('stress fail')
        end %if

      end %do
      
     end
      
  end

      
%% Material update algorithm
      if usefsolve==1 % Use trust region solver
          
%       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-reflective','TolFun',1e-9,'Display','iter','MaxIter',50);
%       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-reflective');
      if norm(stress) == 0
          stress_typical = ones(1,6);
      else
          stress_typical = stress;
      end
      options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-dogleg','TolFun',atol,'MaxIter',50,'Display','off','TypicalX',[0.1*norm(stress_typical(1:6))*ones(1,6) reshape(tt,1,props.num_hard)]');%,'Display','iter','DerivativeCheck','on'
      [x,fval,exitflag,output,J] = fsolve(@(x) mm10_formRfsolve(x, props,np1,n,vec1,vec2,arr1,arr2), x, options);
      
        iter = output.iterations;
        rrr = rcond(J);
        if rrr < 1e-15
            J;
        end
          nR1 = norm(fval(1:6));
          nR2 = norm(fval(7:props.num_hard+6));
      
        if exitflag == 0 || (exitflag < 0 && nR1 > atol1 && nR2 > atol2)
          fail = 1; % since .true. = 1
          return;%error('stress fail')
        end %if
%         % get the hardening variable rates
%       [props, np1, n, vec1, vec2, x(1:6), x(7:props.num_hard+6), R] = mm10_formR(props, np1,...
%           n, vec1, vec2, x(1:6), x(7:props.num_hard+6));
        
        
      elseif usefsolve == 2 % finite difference solver
          
          
      if norm(stress) == 0
          stress_typical = ones(1,6);
      else
          stress_typical = stress;
      end
      options = optimoptions('fsolve','Jacobian','off','Algorithm','trust-region-dogleg','TolFun',atol,'MaxIter',50,'Display','off','TypicalX',[0.1*norm(stress_typical(1:6))*ones(1,6) 0.1*norm(tt)*ones(1,props.num_hard)]');%,'Display','iter'
      [x,fval,exitflag,output,J] = fsolve(@(x) mm10_formRfsolve(x, props,np1,n,vec1,vec2,arr1,arr2), x, options);
      
        iter = output.iterations;
        rrr = rcond(J);
        if rrr < 1e-15
            J;
        end
          nR1 = norm(fval(1:6));
          nR2 = norm(fval(7:props.num_hard+6));
      
        if exitflag == 0 || (exitflag < 0 && nR1 > atol1 && nR2 > atol2)
          fail = 1; % since .true. = 1
          return;%error('stress fail')
        end %if
      
        
      else % Original method ( WARP3D)
          
          
      iter = 0;
      [props, np1, n, vec1, vec2, x(1:6), x(7:props.num_hard+6), R] = mm10_formR(props, np1,...
          n, vec1, vec2, x(1:6), x(7:props.num_hard+6));
      nR1 = norm(R(1:6));
      nR2 = norm(R(7:props.num_hard+6));
      if props.h_type <= 3 || (props.h_type ==6)% DO NOT reset the initial stress residual norm
      inR1 = nR1;
      end
      inR2 = nR2;
      if (debug) 
            fprintf('Iter %i  norm %10.3e \n)',iter, nR1, nR2)
      end %if
      while ((nR1 > atol1) && (nR1/inR1 > rtol1)) || ((nR2 > atol2) && (nR2/inR2 > rtol2)) || iter < mmin
% c           Jacobian
        if usefsolve == 3 % complex derivative solver
        [props, np1, n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6), J] = mm10_formJi(props, np1,...
            n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6));
        elseif usefsolve == 4 % numerically verify real finite differences
        [props, np1, n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6), J] = mm10_formJr(props, np1,...
            n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6));
        [props, np1, n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6), Jcheck] = mm10_formJ(props, np1,...
            n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6));
            if norm(J-Jcheck)/norm(J) > Jtol2
                display('finite difference check failed: norm(J-Jcheck)/norm(J) is too large')
                keyboard
            else
                J = Jcheck;
            end
        else
        [props, np1, n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6), J] = mm10_formJ(props, np1,...
            n,vec1,vec2,arr1,arr2, x(1:6), x(7:props.num_hard+6));
        end
% c           Increment
        dx = R;
        rrr = rcond(J);
        if rrr < 1e-15
            J;
        end
        dx = -J\dx;%[7, 1, -J, 7, ipiv, dx, 7, info] = DGESV(7, 1, -J,...
            %7, ipiv, dx, 7, info); % http://tinyurl.com/oqk27tm
% c                        % computes the solution to a real system of linear equations A * X = B, where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
% c           Line search
        alpha = 1.0;
        ls1 = 0.5*norm(R)^2;
        dy = (transpose(J)*R);
        ls2 = c*dx'*dy;
        ls = 0;
        tr = 1;
        while tr == 1 && ls < mls
          ls = ls + 1;
          nlsx = ls1 + ls2*alpha;
          xnew = x + alpha*dx';
% c             Residual
          [props, np1, n,vec1,vec2, xnew(1:6), xnew(7:props.num_hard+6), R] = ...
              mm10_formR(props, np1, n,vec1,vec2, xnew(1:6), xnew(7:props.num_hard+6));
          nR1 = norm(R(1:6));
          nR2 = norm(R(7:props.num_hard+6));
          nRs = 0.5*norm(R)^2;
% c
          if ((nRs <= nlsx) || (ls >= mls)) % || (cuts > 6 && iter <= 3)
            x = xnew;
            tr = 0;
            break
          else
            alpha = red*alpha;
          end %if
        end %do
%         if ((jiter > mls))
%           fail = 1; % since .true. = 1
%           error('tr fail')
%         end %if
% c
        iter = iter + 1;
        if (debug) 
            fprintf('Iter %i  norm %10.3e ls %10.3f \n',iter, nR1, alpha)
     
        end %if
% c           Increment and check for failure
        if ((iter > miter) || any(isnan(x)) || any(imag(x)))
            if  (cuts > 5 && nR1 < 5 && nR2 < 5e4)
                break
            else
          fail = 1; % since .true. = 1
          return;%error('stress fail')
            end
        end %if

      end %do
      
      end % solver
% 
% c           Set for return
      stress = x(1:6);
      tt = x(7:props.num_hard+6);
% 
%       return
end 