% c
% c
% c     MRR:
% c
% c           Actual Roters hardening function
function [props, np1, n,vec1,vec2, stress, tt, h] = mm10_h_ornl(props, np1,...
    n,vec1,vec2, stress, tt)

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
        G0 = props.mu_0;
        b = props.burgers;
        v = 0.3;%props.nu;
        dt = np1.tinc;
        
% c     New shear modulus
        G = props.mu_0 - props.D_0 / (exp(props.T_0/...
          np1.temp)- 1.0);
        % Load the interaction matrices for parallel and forest dislocs
        [Gmat,Hmat] = mm10_mrr_GH(props.s_type);

% c
%       h = zeros(1,12);

          % Get dislocation density
          rho = tt; % rho^a_SSD
          rho_n = transpose(n.tau_tilde(1:props.num_hard)); % rho^a_SSD

          rs = transpose(stress*np1.ms); % tau^a
          rs = sign(real(rs)).*rs;
          
%           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = Gmat*reshape(tt(1:props.num_hard),props.num_hard,1);
          rhoP = Hmat*reshape(tt(1:props.num_hard),props.num_hard,1);
          
          tpass = c1*G*b*sqrt(rhoP); % (16)
          ddipole = sqrt(3)*G*b/(16*pi*(1-v))./transpose(rs); % (42)
          rhoM = (2*k/(c1*c2*c3*G*b^3))*theta*transpose(sqrt(rhoF.*rhoP)); % (13)
          
          slipinc = vec1(1:props.nslip);%mm10_slipinc_mrr(props, np1, n, stress, tt, alpha);
%           if real(slipinc) < 0
              slipinc = sign(real(slipinc)).*slipinc;
%           end
          gammadot = transpose(slipinc)/dt;
              
        % Evaluate the hardening equation
          h = rho_n + dt*(c4/b*transpose(sqrt(rhoP)).*gammadot ...
              + c6*ddipole./b.*rhoM.*gammadot - c5*rho.*gammadot ...
              - c7*exp(-Qbulk/k/theta).*transpose(rs)./(k*theta).*rho.^2.*gammadot.^c8); % (18)
%       for alpha = 1:12
% 
%           % Get dislocation density
%           rho = tt(alpha); % rho^a_SSD
%           rho_n = n.tau_tilde(alpha); % rho^a_SSD
% 
%           ms = np1.ms(1:6,alpha);
%           rs = stress*ms; % tau^a
%           if real(rs) < 0
%               rs = -rs;
%           end
%           
% %           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
%           rhoF = Gmat(alpha,:)*reshape(tt(1:12),12,1);
%           rhoP = Hmat(alpha,:)*reshape(tt(1:12),12,1);
%           
%           tpass = c1*G*b*sqrt(rhoP); % (16)
%           ddipole = sqrt(3)*G*b/(16*pi*(1-v))/(rs-tpass); % (42)
%           rhoM = (2*k/(c1*c2*c3*G*b^3))*theta*sqrt(rhoF*rhoP); % (13)
%           
%           slipinc = vec1(alpha);%mm10_slipinc_mrr(props, np1, n, stress, tt, alpha);
%           if real(slipinc) < 0
%               slipinc = -slipinc;
%           end
%           gammadot = slipinc/dt;
%               
%         % Evaluate the hardening equation
%           h(alpha) = rho_n + dt*(c4/b*sqrt(rhoP)*gammadot ...
%               + c6*ddipole/b*rhoM*gammadot - c5*rho*gammadot ...
%               - c7*exp(-Qbulk/k/theta)*rs/(k*theta)*rho^2*gammadot^c8); % (18)
%       end
end