% c
% c
% c     halite:
% c
% c           Derivative of hardening fn wrt stress
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, et] = mm10_estress_halite(props,...
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
%           [~, ~, ~, ~, ~, dslip] = mm10_dgdt_halite(props, np1,...
%             n, stress, tt);
        dslip = arr1(1:props.num_hard,1)';

        % Load some material parameters
        if props.s_type == 6
            Qbulk = props.cp_028*ones(12,1);
            c1 = props.cp_001*ones(12,1);
            c2 = props.cp_004*ones(12,1);
            c3 = props.cp_007*ones(12,1);
            c4 = props.cp_010*ones(12,1);
            c5 = props.cp_013*ones(12,1);
            c6 = props.cp_016*ones(12,1);
            c7 = props.cp_019*ones(12,1);
            c8 = props.cp_022*ones(12,1);
        elseif props.s_type == 11
            Qbulk = [props.cp_028*ones(6,1);props.cp_029*ones(6,1);props.cp_030*ones(12,1)];
            c1 = [props.cp_001*ones(6,1);props.cp_002*ones(6,1);props.cp_003*ones(12,1)];
            c2 = [props.cp_004*ones(6,1);props.cp_005*ones(6,1);props.cp_006*ones(12,1)];
            c3 = [props.cp_007*ones(6,1);props.cp_008*ones(6,1);props.cp_009*ones(12,1)];
            c4 = [props.cp_010*ones(6,1);props.cp_011*ones(6,1);props.cp_012*ones(12,1)];
            c5 = [props.cp_013*ones(6,1);props.cp_014*ones(6,1);props.cp_015*ones(12,1)];
            c6 = [props.cp_016*ones(6,1);props.cp_017*ones(6,1);props.cp_018*ones(12,1)];
            c7 = [props.cp_019*ones(6,1);props.cp_020*ones(6,1);props.cp_021*ones(12,1)];
            c8 = [props.cp_022*ones(6,1);props.cp_023*ones(6,1);props.cp_024*ones(12,1)];
        else
            error(' ');
        end
        k = props.boltzman;
        theta = np1.temp;
        G = props.mu_0;
        b = props.burgers;
        v = 0.3;%props.nu;
        dt = np1.tinc;
        
        % Compute the shear modulus using Roter's function
        [G_110, G_100, G_111] = mm10_halite_Gtemp(theta);
        G_all = [G_110*ones(6,1); G_100*ones(6,1); G_111*ones(12,1)];
%         if G < 0
%             G = -G;
%         else
%         G = mm10_halite_Gtemp(theta);
%         end
        % Load the interaction matrices for parallel and forest dislocs
        [Gmat,Hmat] = mm10_halite_GH(props.s_type);
% c
      for alpha = 1:props.num_hard
          
          G = G_all(alpha);
          % Get dislocation density
          rho = tt(alpha); % rho^a_SSD

          ms = np1.ms(1:6,alpha);
          rs = stress*ms; % tau^a
          
%           [rhoF,rhoP] = mm10_rhoFP_halite(props, np1, n, tt, alpha);
          rhoF = Gmat(alpha,:)*reshape(tt(1:props.num_hard),props.num_hard,1);
          rhoP = Hmat(alpha,:)*reshape(tt(1:props.num_hard),props.num_hard,1);
          
          tpass = c1(alpha)*G*b*sqrt(rhoP); % (16)
%================================Ran make ddipole = 0======================
%          ddipole = sqrt(3)*G*b/(16*pi*(1-v))/(abs(rs)); % (42)
%          dddipole = -sqrt(3)*G*b/(16*pi*(1-v))/(abs(rs))^2*sign(rs);
          ddipole = 0;
          dddipole = 0;
%================================Ran make ddipole = 0======================
rhoM = (2*k/(c1(alpha)*c2(alpha)*c3(alpha)*G*b^3))*theta*sqrt(rhoF*rhoP); % (13)
          
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
          badterm = c8(alpha)*abs(rs)*gammadot^(c8(alpha)-1)*dslipinc*sign(rs);
          end
          et(alpha,1:6) = dt*(c4(alpha)/b*sqrt(rhoP)*sign(rs)*dslipinc ...
          + c6(alpha)/b*rhoM*(ddipole*sign(rs)*dslipinc + dddipole*gammadot) - c5(alpha)*rho*sign(rs)*dslipinc ...
          - c7(alpha)*exp(-Qbulk(alpha)/k/theta)/(k*theta)*rho^2*(badterm + sign(rs)*gammadot^c8(alpha))) ...
          *dtdstress; % d(18)/dstress
      
      end
      
% 
% c
%       return
end