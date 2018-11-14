% c
% c
% c     halite:
% c
% c           Derivative of hardening fn wrt hardening
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, etau] = mm10_ehard_halite(props,...
    np1, n,vec1,vec2,arr1,arr2, stress, tt)
%       use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision :: tt, etau
% c
%       double precision :: mm10_slipinc
% c
%       double precision :: ur, s, ct, cta
%       integer :: i
% c
      etau = zeros(props.num_hard,props.num_hard);
      
      % compute derivatives of slip increments with respect to densities
%           [~, ~, ~, ~, ~, dslip] = mm10_dgdh_halite(props, np1,...
%             n, stress, tt);
       dslip = arr2(1:props.nslip,1:props.num_hard);

        % Load some material parameters
        k = props.boltzman;
        theta = np1.temp;
        G = props.mu_0;
        b = props.burgers;
        v = 0.3;%props.nu;
        dt = np1.tinc;
        
        % Compute the shear modulus using Roter's function
%         if G < 0
%             G = -G;
%         else
%         G = mm10_halite_Gtemp(theta);
%         end
        % Load the interaction matrices for parallel and forest dislocs
        [Gmat,Hmat] = mm10_halite_GH(props.s_type);
% c
    % Compute drho_alpha/drho_beta
    % loop over numerator hardening variable
    for alpha = 1:props.num_hard
        % select the material parameter for slip system alpha
        if alpha <=6 || props.s_type == 6
            Qbulk = props.cp_028;
            c1 = props.cp_001;
            c2 = props.cp_004;
            c3 = props.cp_007;
            c4 = props.cp_010;
            c5 = props.cp_013;
            c6 = props.cp_016;
            c7 = props.cp_019;
            c8 = props.cp_022;
            [G, ~, ~] = mm10_halite_Gtemp(theta);
        elseif alpha <= 12
            Qbulk = props.cp_029;
            c1 = props.cp_002;
            c2 = props.cp_005;
            c3 = props.cp_008;
            c4 = props.cp_011;
            c5 = props.cp_014;
            c6 = props.cp_017;
            c7 = props.cp_020;
            c8 = props.cp_023;
            [~, G, ~] = mm10_halite_Gtemp(theta);
        else
            Qbulk = props.cp_030;
            c1 = props.cp_003;
            c2 = props.cp_006;
            c3 = props.cp_009;
            c4 = props.cp_012;
            c5 = props.cp_015;
            c6 = props.cp_018;
            c7 = props.cp_021;
            c8 = props.cp_024;
            [~, ~, G] = mm10_halite_Gtemp(theta);
        end

        % Get dislocation density
        rho = tt(alpha); % rho^a_SSD

        ms = np1.ms(1:6,alpha);
        rs = stress*ms; % tau^a
          
%         [rhoF,rhoP] = mm10_rhoFP_halite(props, np1, n, tt, alpha);
        rhoF = Gmat(alpha,:)*reshape(tt(1:props.num_hard),props.num_hard,1);
        rhoP = Hmat(alpha,:)*reshape(tt(1:props.num_hard),props.num_hard,1);
          
        slipinc = vec1(alpha);
        gammadot = abs(slipinc/dt);
          
        tpass = c1*G*b*sqrt(rhoP); % (16)
%        ddipole = sqrt(3)*G*b/(16*pi*(1-v))/(abs(rs)); % (42)
        rhoM = (2*k/(c1*c2*c3*G*b^3))*theta*sqrt(rhoF*rhoP); % (13)
          
        % loop over denominator hardening variable
        for beta = 1:props.num_hard
          
%           [drhoF,drhoP] = mm10_drhoFP_halite(props, np1, n, tt, alpha, beta);
          drhoF = Gmat(alpha,beta);  % Ran: d(rhoF(alpha))/d(rho(beta))
          drhoP = Hmat(alpha,beta);  % Ran: d(rhoP(alpha))/d(rho(beta))
          
%          dddipole = -sqrt(3)*G*b/(16*pi*(1-v))/(abs(rs))^2 * (-0.5*c1*G*b/sqrt(rhoP)*drhoP);
%================================Ran make ddipole = 0======================
          ddipole = 0;
          dddipole = 0;
%================================Ran make ddipole = 0======================
          drhoM = 0.5*(2*k/(c1*c2*c3*G*b^3))*theta/sqrt(rhoF*rhoP) * (drhoF*rhoP + rhoF*drhoP);
          
          dslipinc = dslip(alpha,beta)/dt;  % Ran: d(gamma_dot(alpha))/d(rho(beta))
          
          if alpha == beta
              deltaij = 1;
          else
              deltaij = 0;
          end

          % Evaluate the equation
          if gammadot == 0
          badterm = 0;
          else
          badterm = c8*gammadot^(c8-1)*dslipinc;
          end
%           h(i) = rho_n + dt*(c4/b*sqrt(rhoP)*gammadot ...
%               + c6*ddipole/b*rhoM*gammadot - c5*rho*gammadot ...
%               - c7*exp(-Qbulk/k/theta)*abs(rs)/(k*theta)*rho^2*gammadot^c8); % (18)
          etau(alpha,beta) = deltaij - dt*(c4/b*(0.5*drhoP/sqrt(rhoP)*gammadot*sign(rs) + sqrt(rhoP)*dslipinc) ...
          + c6/b*(rhoM*ddipole*dslipinc + rhoM*dddipole*gammadot*sign(rs) + drhoM*ddipole*gammadot*sign(rs)) - c5*(rho*dslipinc + deltaij*gammadot*sign(rs)) ...
          - c7*exp(-Qbulk/k/theta)/(k*theta)*abs(rs)*(rho^2*badterm + 2*rho*sign(rs)*deltaij*gammadot^c8))*sign(rs);
      
        end %beta
      
    end %alpha
%       return
end