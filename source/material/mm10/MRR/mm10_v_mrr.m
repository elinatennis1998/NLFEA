% c
% c
% c     MRR:
% c
% c           Form intermediate vectors for faster calculations
function [props, np1, n, stress, tt, vec1, vec2] = mm10_v_mrr(props, np1,...
    n, stress, tt, vec1, vec2)

        % Load some material parameters
        dt = np1.tinc;
        k = props.boltzman;
        theta = np1.temp;
        G = props.mu_0;
        b = props.burgers;
        c1 = props.c1;
        c2 = props.c2;
        c3 = props.c3;
        p_e = props.p_v;
        q_e = props.q_v;
        Qslip = props.Qslip;
        v_attack = props.G_0_y;
        
        % Compute all slip increments
        
        % Compute the shear modulus using Roter's function
        if G < 0
            G = -G;
        else
        G = mm10_mrr_Gtemp(theta);
        end
        % Load the interaction matrices for parallel and forest dislocs
        [Gmat,Hmat] = mm10_mrr_GH(props.s_type);
        
        % Parallel and forest dislocations
          rhoF = Gmat*reshape(tt(1:props.num_hard),props.num_hard,1);
          rhoP = Hmat*reshape(tt(1:props.num_hard),props.num_hard,1);
          rs = transpose(stress*np1.ms); % tau^a
        
%         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, i);

        % Compute some stresses and rates
        gamma_0 = v_attack*k*theta/(c1*c3*G*b^2)*sqrt(rhoP); % (15)
        tpass = c1*G*b*sqrt(rhoP); % (16)
        tcut = Qslip/(c2*c3*b^2)*sqrt(rhoF); % (17)
%           if real(rs) < 0
%         fract = (-rs-tpass)/tcut;
%           else
        fract = (sign(real(rs)).*rs-tpass)./tcut;
%           end
%         fract = ((abs(rs)-tpass)/tcut);

        % Evaluate the slip rate equation
        for s = 1:props.nslip
        if fract(s) > 1
            % linear extrapolation past the too-high stress (rs) value
            b = gamma_0(s);
            x = fract(s);
            m = b * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)^(q_e-1.d0)) ...
               * (- p_e*(1.d0)^(p_e-1.d0)) * sign(real(rs(s)))/tcut(s);
            y = m.*x + b;
            vec1(s) = dt * y;
        elseif fract(s) > 0
        vec1(s) = dt * gamma_0(s) * exp (-(Qslip/k/theta)*(1.d0 - fract(s)^p_e)^q_e) * sign(real(rs(s))); %(14)
        else
            vec1(s) = 0;
        end
        end

end