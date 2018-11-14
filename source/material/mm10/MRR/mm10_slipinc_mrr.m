% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c       from mm10_form.f
% c
% c           Calculate the slip increment along system i  
% c           mrr model
function slipinc = ...
    mm10_slipinc_mrr(props, np1, n, stress, tt, alpha)

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
        
        % Compute the shear modulus using Roter's function
        if G < 0
            G = -G;
        else
        G = mm10_mrr_Gtemp(theta);
        end
        % Load the interaction matrices for parallel and forest dislocs
        [Gmat,Hmat] = mm10_mrr_GH(props.s_type);
        
      ms = np1.ms(1:6,alpha);
      rs = stress*ms; % tau^a
        
%         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, i);
          rhoF = Gmat(alpha,:)*reshape(tt(1:props.num_hard),props.num_hard,1);
          rhoP = Hmat(alpha,:)*reshape(tt(1:props.num_hard),props.num_hard,1);

        % Compute some stresses and rates
        gamma_0 = v_attack*k*theta/(c1*c3*G*b^2)*sqrt(rhoP); % (15)
        tpass = c1*G*b*sqrt(rhoP); % (16)
        tcut = Qslip/(c2*c3*b^2)*sqrt(rhoF); % (17)
          if real(rs) < 0
        fract = (-rs-tpass)/tcut;
          else
        fract = (rs-tpass)/tcut;
          end
%         fract = ((abs(rs)-tpass)/tcut);

        % Evaluate the slip rate equation
        if fract > 1
            % linear extrapolation past the too-high stress (rs) value
            b = gamma_0;
            x = fract;
            m = b * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)^(q_e-1.d0)) ...
               * (- p_e*(1.d0)^(p_e-1.d0)) * sign(real(rs))/tcut;
            y = m*x + b;
            slipinc = dt * y;
        elseif fract > 0
        slipinc = dt * gamma_0 * exp (-(Qslip/k/theta)*(1.d0 - fract^p_e)^q_e) * sign(real(rs)); %(14)
        else
            slipinc = 0;
        end
      
end
