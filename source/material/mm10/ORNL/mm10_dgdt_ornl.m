% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10_dgdt_mts                     *
% c     *                                                              *
% c     *                       written by : tjt                       *
% c     *                                                              *
% c     *                   last modified: 01/26/15                    *
% c     *                                                              *
% c     *     Calculate partial gamma(s) w.r.t. tau(s) for each slip   *
% c     *     system, for use in J11 in material integration. mrr     *
% c     *                                                              *
% c     ****************************************************************
% c
function [props, np1, n, stress, tt, dgdt] = mm10_dgdt_ornl(props, np1,...
            n, stress, tt)
        
      dgdt = zeros(1,12);

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
        
        rs = (stress*np1.ms)'; % tau^a
        
        % Parallel and forest dislocations
          rhoF = Gmat*reshape(tt(1:props.num_hard),props.num_hard,1);
          rhoP = Hmat*reshape(tt(1:props.num_hard),props.num_hard,1);
          
        rhoM = fM*tt;
        gamma_0 = rhoM*b*v_s/b*lamda;
        tpass = c1*G*b*sqrt(rhoP); % (16)
        tcut = tau0*G/G0*ones(props.num_hard,1); % (17)
        fract = (sign(real(rs)).*rs-tpass)./tcut;
        
        dfract = sign(rs)/tcut;
        for s = 1:props.nslip
        if fract(s) > 1
            % linear extrapolation past the too-high stress (rs) value
            b = gamma_0(s);
            m = b * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)^(q_e-1.d0)) ...
               * (- p_e*(1.d0)^(p_e-1.d0)) * sign(rs(s))/tcut(s);
            dslipinc = dt * m;
        elseif fract(s) > 0
        slipexp = exp (-(Qslip/k/theta)*(1.d0 - fract(s)^p_e)^q_e) * sign(rs(s));
        dslipinc = dt * (gamma_0(s) * slipexp * -(Qslip/k/theta)*q_e*(1.d0 - fract(s)^p_e)^(q_e-1.d0) ...
            * -p_e*fract(s)^(p_e-1.d0) * dfract(s)); %(14)
        else
            dslipinc = 0;
        end
        dgdt(s) = dslipinc;
        end
%       for alpha = 1:12
%         
%         ms = np1.ms(1:6,alpha);
%         rs = stress*ms; % tau^a
%         
% %           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
%           rhoF = Gmat(alpha,:)*reshape(tt(1:12),12,1);
%           rhoP = Hmat(alpha,:)*reshape(tt(1:12),12,1);
% 
%         % Compute one dependency
%         gamma_0 = v_attack*k*theta/(c1*c3*G*b^2)*sqrt(rhoP); % (15)
%         tpass = c1*G*b*sqrt(rhoP); % (16)
%         tcut = Qslip/(c2*c3*b^2)*sqrt(rhoF); % (17)
% 
%         % Evaluate the equation
%         fract = ((abs(rs)-tpass)/tcut);
%         dfract = sign(rs)/tcut;
%         if fract > 1
%             % linear extrapolation past the too-high stress (rs) value
%             b = gamma_0;
%             m = b * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)^(q_e-1.d0)) ...
%                * (- p_e*(1.d0)^(p_e-1.d0)) * sign(rs)/tcut;
%             dslipinc = dt * m;
%         elseif fract > 0
%         slipexp = exp (-(Qslip/k/theta)*(1.d0 - fract^p_e)^q_e) * sign(rs);
%         dslipinc = dt * (gamma_0 * slipexp * -(Qslip/k/theta)*q_e*(1.d0 - fract^p_e)^(q_e-1.d0) ...
%             * -p_e*fract^(p_e-1.d0) * dfract); %(14)
%         else
%             dslipinc = 0;
%         end
%         
%         dgdt(alpha) = dslipinc;
% 
%       end
        
end