% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdh_mts                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 01/26/15                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. tt for each slip       *
% c     *     system, for use in J12 in material integration. mrr     *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, dgammadtt] = mm10_dgdh_ornl(props, np1,...
            n, stress, tt)
        
      dgammadtt = zeros(props.num_hard,props.num_hard);

        % Load some material parameters
      dt = np1.tinc;
        k = props.boltzman;
        theta = np1.temp;
        G0 = props.mu_0;
        b = props.burgers;
        c1 = props.c1;
        c2 = props.c2;
        tau0 = props.c3;
        p_e = props.p_v;
        q_e = props.q_v;
        Qslip = props.Qslip;
        v_s = props.G_0_y;
        fM = 0.1;
        lamda = c2*b;
        
        
% c     New shear modulus
        G = props.mu_0 - props.D_0 / (exp(props.T_0/...
          np1.temp)- 1.0);

        % Load the interaction matrices for parallel and forest dislocs
        [Gmat,Hmat] = mm10_mrr_GH(props.s_type);
        
      % Compute derivative of slip rate alpha w.r.t. density beta
      % loop over slip rate
        rs = (stress*np1.ms)'; % tau^a
        
        % Parallel and forest dislocations
          rhoF = Gmat*reshape(tt(1:props.num_hard),props.num_hard,1);
          rhoP = Hmat*reshape(tt(1:props.num_hard),props.num_hard,1);
          
        rhoM = fM*tt;
        gamma_0 = rhoM*b*v_s/b*lamda;
        tpass = c1*G*b*sqrt(rhoP); % (16)
        tcut = tau0*G/G0*ones(props.num_hard,1); % (17)
        fract = (sign(real(rs)).*rs-tpass)./tcut;
        
          drhoF = Gmat;
          drhoP = Hmat;
          
             drhoM = fM*eye(props.num_hard);
          dgamma_0 = drhoM*b*v_s/b*lamda; % (15)
          dtpass = 0.5*c1*G*b./(sqrt(rhoP)*ones(1,props.num_hard)).*drhoP; % (16)
          dtcut = 0.0*drhoF; % (17)
          dfract = ((-dtpass).*(tcut*ones(1,props.num_hard)) - ((abs(rs)-tpass)*ones(1,props.num_hard)).*dtcut)./(tcut.^2.d0*ones(1,props.num_hard));
      for alpha = 1:props.num_hard
        for beta = 1:props.num_hard
          % Evaluate the equation
          if fract(alpha) > 1
              % linear extrapolation past the too-high stress (rs) value
              dslipinc = dt * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)^(q_e-1.d0)) ...
                 * (- p_e*(1.d0)^(p_e-1.d0)) * sign(rs(alpha)) * (dgamma_0(alpha,beta)/tcut(alpha) -gamma_0(alpha)*dtcut(alpha,beta)/tcut(alpha)^2.d0);
          elseif fract(alpha) > 0
          slipexp = exp (-(Qslip/k/theta)*(1.d0 - fract(alpha)^p_e)^q_e); %(14)
          dslipinc = dt * (dgamma_0(alpha,beta) * slipexp * sign(rs(alpha)) + gamma_0(alpha) * sign(rs(alpha)) * slipexp * -(Qslip/k/theta)*q_e*(1.d0 - fract(alpha)^p_e)^(q_e-1.d0) ...
              * -p_e*fract(alpha)^(p_e-1.d0) * dfract(alpha,beta)); %gamma_0 * exp (-(Qslip/k/theta)*(1.d0 - fract^p_e)^q_e) * sign(rs); %(14)
          else
              dslipinc = 0;
          end
        
          dgammadtt(alpha,beta) = dslipinc;
        end
      end
        
%       % Compute derivative of slip rate alpha w.r.t. density beta
%       % loop over slip rate
%       for alpha = 1:12
%         
%         ms = np1.ms(1:6,alpha);
%         rs = stress*ms; % tau^a
%         
% %           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
%           rhoF = Gmat(alpha,:)*reshape(tt(1:12),12,1);
%           rhoP = Hmat(alpha,:)*reshape(tt(1:12),12,1);
%           
%         gamma_0 = v_attack*k*theta/(c1*c3*G*b^2)*sqrt(rhoP); % (15)
%         tpass = c1*G*b*sqrt(rhoP); % (16)
%         tcut = Qslip/(c2*c3*b^2)*sqrt(rhoF); % (17)
%         fract = ((abs(rs)-tpass)/tcut);
%         
%         % loop over density
%         for beta = 1:12
%         
% %           [drhoF,drhoP] = mm10_drhoFP_mrr(props, np1, n, tt, alpha, beta);
%           drhoF = Gmat(alpha,beta);
%           drhoP = Hmat(alpha,beta);
%         
%           dgamma_0 = 0.5*v_attack*k*theta/(c1*c3*G*b^2)/sqrt(rhoP)*drhoP; % (15)
%           dtpass = 0.5*c1*G*b/sqrt(rhoP)*drhoP; % (16)
%           dtcut = 0.5*Qslip/(c2*c3*b^2)/sqrt(rhoF)*drhoF; % (17)
%           dfract = ((-dtpass)*tcut - (abs(rs)-tpass)*dtcut)/tcut^2.d0;
% 
%           % Evaluate the equation
%           if fract > 1
%               % linear extrapolation past the too-high stress (rs) value
%               dslipinc = dt * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)^(q_e-1.d0)) ...
%                  * (- p_e*(1.d0)^(p_e-1.d0)) * sign(rs) * (dgamma_0/tcut -gamma_0*dtcut/tcut^2.d0);
%           elseif fract > 0
%           slipexp = exp (-(Qslip/k/theta)*(1.d0 - fract^p_e)^q_e); %(14)
%           dslipinc = dt * (dgamma_0 * slipexp * sign(rs) + gamma_0 * sign(rs) * slipexp * -(Qslip/k/theta)*q_e*(1.d0 - fract^p_e)^(q_e-1.d0) ...
%               * -p_e*fract^(p_e-1.d0) * dfract); %gamma_0 * exp (-(Qslip/k/theta)*(1.d0 - fract^p_e)^q_e) * sign(rs); %(14)
%           else
%               dslipinc = 0;
%           end
%         
%           dgammadtt(alpha,beta) = dslipinc;
% 
%         end %beta
%       
%       end %alpha
        
end