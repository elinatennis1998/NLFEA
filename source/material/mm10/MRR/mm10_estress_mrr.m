% c
% c
% c     MRR:
% c
% c           Derivative of hardening fn wrt stress
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, et] = mm10_estress_mrr(props,...
    np1, n,vec1,vec2,arr1,arr2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt
%       double precision, dimension(6) :: et
% c
%       double precision :: mm10_rs
% c
%       double precision :: rs, cta, ct
%       integer :: i
      
      et = zeros(props.num_hard,6);
      % compute derivatives of slip increments with respect to resolved
      % shear stress
%           [~, ~, ~, ~, ~, dslip] = mm10_dgdt_mrr(props, np1,...
%             n, stress, tt);
        dslip = arr1(1:props.num_hard,1)';

        % Load some material parameters
        Qbulk = props.Qbulk;
        k = props.boltzman;
        theta = np1.temp;
        c1 = props.c1;
        c2 = props.c2;
        c3 = props.c3;
        c4 = props.c4;
        c5 = props.c5;
        c6 = props.c6;
        c7 = props.c7;
        c8 = props.c8;
        G = props.mu_0;
        b = props.burgers;
        v = 0.3;%props.nu;
        dt = np1.tinc;
        
        % Compute the shear modulus using Roter's function
        if G < 0
            G = -G;
        else
        G = mm10_mrr_Gtemp(theta);
        end
        % Load the interaction matrices for parallel and forest dislocs
        [Gmat,Hmat] = mm10_mrr_GH(props.s_type);
% c
      for alpha = 1:props.num_hard

          % Get dislocation density
          rho = tt(alpha); % rho^a_SSD

          ms = np1.ms(1:6,alpha);
          rs = stress*ms; % tau^a
          
%           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = Gmat(alpha,:)*reshape(tt(1:props.num_hard),props.num_hard,1);
          rhoP = Hmat(alpha,:)*reshape(tt(1:props.num_hard),props.num_hard,1);
          
          tpass = c1*G*b*sqrt(rhoP); % (16)
          ddipole = sqrt(3)*G*b/(16*pi*(1-v))/(abs(rs)); % (42)
          dddipole = -sqrt(3)*G*b/(16*pi*(1-v))/(abs(rs))^2*sign(rs);
          rhoM = (2*k/(c1*c2*c3*G*b^3))*theta*sqrt(rhoF*rhoP); % (13)
          
          slipinc = vec1(alpha);
          gammadot = abs(slipinc/dt);
          dtdstress = ms;
          dslipinc = dslip(alpha)/dt; % always positive

          % Evaluate the equation
%           h(i) = rho_n + dt*(c4/b*sqrt(rhoP)*gammadot ...
%               + c6*ddipole/b*rhoM*gammadot - c5*rho*gammadot ...
%               - c7*exp(-Qbulk/k/theta)*abs(rs)/(k*theta)*rho^2*gammadot^c8); % (18)
          if gammadot == 0
          badterm = 0;
          else
          badterm = c8*abs(rs)*gammadot^(c8-1)*dslipinc*sign(rs);
          end
          et(alpha,1:6) = dt*(c4/b*sqrt(rhoP)*sign(rs)*dslipinc ...
          + c6/b*rhoM*(ddipole*sign(rs)*dslipinc + dddipole*gammadot) - c5*rho*sign(rs)*dslipinc ...
          - c7*exp(-Qbulk/k/theta)/(k*theta)*rho^2*(badterm + sign(rs)*gammadot^c8)) ...
          *dtdstress; % d(18)/dstress
      
      end
      
% 
% c
%       return
end