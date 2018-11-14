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
function [props, np1, n, stress, tt, fail] = mm10_solve2(props, np1,...
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
atol1 = 1.0e-6;1;
rtol1 = 1.0e-12;
atol2 = 1e-1;4e4;1.0e-7*max(tt);
rtol2 = 1.0e-12;
atol = 1e-7;
miter = 30;15;250;
mmin = 3;
debug = 0;
c = 1.0e-4;
red = 0.5;
mls = 4;3;6;10;

usefsolve = 0;1;

if props.h_type > 3 && (props.h_type ~=6)% new models except omars
    % c
    if (debug)
        disp('Entering solution routine')
    end %if
    x(1:6) = stress;
    x(7:props.num_hard+6) = tt;
    %       dt = np1.tinc;
    %       % Predict the new hardening variable by extrapolation
    %       [props, np1, n, x(1:6), x(7:props.num_hard+6), tt_h] = mm10_formR2(props, np1,...
    %           n, x(1:6), x(7:props.num_hard+6));
    %       tt_new = tt - tt_h;
    %       rate = (tt_new - n.tau_tilde)/dt;
    %       x(7:props.num_hard+6) = tt + rate*dt*step;
    
    
    % Initial search for stress
    iter = 0;
    [props, np1, n, x(1:6), x(7:props.num_hard+6), R] = mm10_formR1(props, np1,...
        n, x(1:6), x(7:props.num_hard+6));
    nR1 = sqrt(dot(R(1:6),R(1:6)));
    inR1 = nR1;

    while ((nR1 > atol1) && (nR1/inR1 > rtol1))
        % c           Jacobian
        [props, np1, n, x(1:6), x(7:props.num_hard+6), J] = mm10_formJ11(props, np1,...
            n, x(1:6), x(7:props.num_hard+6));
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
            [props, np1, n, xnew(1:6), xnew(7:props.num_hard+6), R] = ...
                mm10_formR1(props, np1, n, xnew(1:6), xnew(7:props.num_hard+6));
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
            fprintf('Iter %i  norm %10.3e ls %10.3f \n',iter, nR, alpha)
            
        end %if
        % c           Increment and check for failure
        if ((iter > miter) || any(isnan(x)) || any(imag(x)))
            fail = 1; % since .true. = 1
            return;%error('stress fail')
        end %if
        
    end %do
    
    
    if usefsolve % Use trust region solver
        
        %       str = x(1:6);
        %       ttn = x(7:props.num_hard+6);
        %       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-reflective','TolFun',1e-9,'Display','off','MaxIter',60);
        % %       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-dogleg','TolFun',tolll,'Display','iter','MaxIter',60,'TypicalX',[tt]');%,'Display','iter','TypicalX',[0.1*tt*ones(1,6) tt]'
        %       [ttn,~,exitflag] = fsolve(@(ttn) mm10_formRfsolve2(ttn, props,np1,n,str), ttn, options);
        %
        %         if exitflag <= 0
        %           fail = 1; % since .true. = 1
        %           return;%error('stress fail')
        %         end %if
        %
        %       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-reflective','TolFun',1e-9,'Display','off');
        % %       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-dogleg','TolFun',1e-9,'Display','iter','TypicalX',[0.1*norm(stress(1:6))*ones(1,6)]');%,'Display','iter','TypicalX',[0.1*tt*ones(1,6) tt]'
        %       [str,~,exitflag] = fsolve(@(str) mm10_formRfsolve1(str, props,np1,n,ttn), str, options);
        %
        %         if exitflag <= 0
        %           fail = 1; % since .true. = 1
        %           return;%error('stress fail')
        %         end %if
        %
        %       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-reflective','TolFun',1e-9,'Display','off','MaxIter',60);
        % %       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-dogleg','TolFun',tolll,'Display','iter','MaxIter',40,'TypicalX',[tt]');%,'Display','iter','TypicalX',[0.1*tt*ones(1,6) tt]'
        %       [ttn,~,exitflag] = fsolve(@(ttn) mm10_formRfsolve2(ttn, props,np1,n,str), ttn, options);
        %
        %         if exitflag <= 0 && exitflag ~= -2
        %           fail = 1; % since .true. = 1
        %           return;%error('stress fail')
        %         end %if
        %
        %       x(1:6) = str;
        %       x(7:props.num_hard+6) = ttn;
        
        %       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-reflective','TolFun',1e-9,'Display','iter','MaxIter',50);
        %       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-reflective');
        options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-dogleg','TolFun',atol,'MaxIter',50,'Display','iter','TypicalX',[0.1*norm(stress(1:6))*ones(1,6) tt]');%,'Display','iter'
        [x,fval,exitflag,output,J] = fsolve(@(x) mm10_formRfsolve(x, props,np1,n), x, options);
        
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
        
    else % Original method
        % c
        %       x(1:6) = stress;
        %       x(7:props.num_hard+6) = tt;
        iter = 0;
        [props, np1, n, x(1:6), x(7:props.num_hard+6), R] = mm10_formR(props, np1,...
            n, x(1:6), x(7:props.num_hard+6));
        %       nR1 = sqrt(dot(R(1:6),R(1:6)));
        %       inR1 = nR1;
        nR2 = sqrt(dot(R(7:props.num_hard+6),R(7:props.num_hard+6)));
        inR2 = nR2;
        if (debug)
            fprintf('Iter %i  norm %10.3e \n)',iter, nR)
        end %if
        while ((nR1 > atol1) && (nR1/inR1 > rtol1)) || ((nR2 > atol2) && (nR2/inR2 > rtol2)) && iter < mmin
            % c           Jacobian
            [props, np1, n, x(1:6), x(7:props.num_hard+6), J] = mm10_formJ(props, np1,...
                n, x(1:6), x(7:props.num_hard+6));
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
            %         ls11 = 0.5*dot(R(1:6),R(1:6));
            %         ls12 = 0.5*dot(R(7:props.num_hard+6),R(7:props.num_hard+6));
            ls1 = 0.5*dot(R,R);
            dy = (transpose(J)*R);
            %         ls21 = c*dot(dx(1:6), dy(1:6));
            %         ls22 = c*dot(dx(7:props.num_hard+6), dy(7:props.num_hard+6));
            ls2 = c*dot(dx, dy);
            ls = 0;
            tr = 1;
            while tr == 1 && ls < mls
                ls = ls + 1;
                %           nlsx1 = ls11 + ls21*alpha;
                %           nlsx2 = ls12 + ls22*alpha;
                nlsx = ls1 + ls2*alpha;
                xnew = x + alpha*dx';
                % c             Residual
                [props, np1, n, xnew(1:6), xnew(7:props.num_hard+6), R] = ...
                    mm10_formR(props, np1, n, xnew(1:6), xnew(7:props.num_hard+6));
                %           nR1 = sqrt(dot(R(1:6),R(1:6)));
                %           nR2 = sqrt(dot(R(7:props.num_hard+6),R(7:props.num_hard+6)));
                %           nRs1 = 0.5*dot(R(1:6),R(1:6));
                %           nRs2 = 0.5*dot(R(7:props.num_hard+6),R(7:props.num_hard+6));
                nR1 = sqrt(dot(R(1:6),R(1:6)));
                nR2 = sqrt(dot(R(7:props.num_hard+6),R(7:props.num_hard+6)));
                nRs = 0.5*dot(R,R);
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
                fprintf('Iter %i  norm %10.3e ls %10.3f \n',iter, nR, alpha)
                
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
    
    
else % Mark's model
    
    
    % c
    if (debug)
        disp('Entering solution routine')
    end %if
    % c
    x(1:6) = stress;
    x(7:props.num_hard+6) = tt;
    
    if usefsolve % Use trust region solver
        
        %       options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-reflective');
        options = optimoptions('fsolve','Jacobian','on','Algorithm','trust-region-dogleg','TolFun',1e-9,'Display','off');%,'Display','iter','TypicalX',[0.1*tt*ones(1,6) tt]'
        [x,~,exitflag] = fsolve(@(x) mm10_formRfsolve(x, props,np1,n), x, options);
        
        if exitflag <= 0
            fail = 1; % since .true. = 1
            return;%error('stress fail')
        end %if
        
    else % Original method
        
        iter = 0;
        [props, np1, n, x(1:6), x(7:props.num_hard+6), R] = mm10_formR(props, np1,...
            n, x(1:6), x(7:props.num_hard+6));
        nR1 = sqrt(dot(R(1:6),R(1:6)));
        inR1 = nR1;
        nR2 = sqrt(dot(R(7:props.num_hard+6),R(7:props.num_hard+6)));
        inR2 = nR2;
        if (debug)
            fprintf('Iter %i  norm %10.3e \n)',iter, nR)
        end %if
        while ((nR1 > atol1) && (nR1/inR1 > rtol1)) || ((nR2 > atol2) && (nR2/inR2 > rtol2)) && iter < mmin
            % c           Jacobian
            [props, np1, n, x(1:6), x(7:props.num_hard+6), J] = mm10_formJ(props, np1,...
                n, x(1:6), x(7:props.num_hard+6));
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
            %         ls11 = 0.5*dot(R(1:6),R(1:6));
            %         ls12 = 0.5*dot(R(7:props.num_hard+6),R(7:props.num_hard+6));
            ls1 = 0.5*dot(R,R);
            dy = (transpose(J)*R);
            %         ls21 = c*dot(dx(1:6), dy(1:6));
            %         ls22 = c*dot(dx(7:props.num_hard+6), dy(7:props.num_hard+6));
            ls2 = c*dot(dx, dy);
            ls = 0;
            tr = 1;
            while tr == 1 && ls < mls
                ls = ls + 1;
                %           nlsx1 = ls11 + ls21*alpha;
                %           nlsx2 = ls12 + ls22*alpha;
                nlsx = ls1 + ls2*alpha;
                xnew = x + alpha*dx';
                % c             Residual
                [props, np1, n, xnew(1:6), xnew(7:props.num_hard+6), R] = ...
                    mm10_formR(props, np1, n, xnew(1:6), xnew(7:props.num_hard+6));
                %           nR1 = sqrt(dot(R(1:6),R(1:6)));
                %           nR2 = sqrt(dot(R(7:props.num_hard+6),R(7:props.num_hard+6)));
                %           nRs1 = 0.5*dot(R(1:6),R(1:6));
                %           nRs2 = 0.5*dot(R(7:props.num_hard+6),R(7:props.num_hard+6));
                nR1 = sqrt(dot(R(1:6),R(1:6)));
                nR2 = sqrt(dot(R(7:props.num_hard+6),R(7:props.num_hard+6)));
                nRs = 0.5*dot(R,R);
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
                fprintf('Iter %i  norm %10.3e ls %10.3f \n',iter, nR, alpha)
                
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
    % c           Set for return
    stress = x(1:6);
    tt = x(7:props.num_hard+6);
    
end
%
%       return
end