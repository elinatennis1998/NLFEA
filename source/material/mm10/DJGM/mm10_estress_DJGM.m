% c
% c
% c     DJGM:
% c
% c           Derivative of hardening fn wrt stress
function [props, np1, n,vec1,vec2,arr1,arr2, stress, tt, et] = mm10_estress_DJGM(props,...
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
        dslipinc = arr1(1:props.num_hard,1)';

        % Load material parameters
        dt = np1.tinc;
        dslip = dslipinc/dt;
        slipinc = vec1(1:props.nslip);
        slipinc = sign(real(slipinc)).*slipinc;
        gamma_dot = transpose(slipinc)/dt;
        g = tt;
        gtol = props.tau_hat_y;
        
        [h_0, gamma_dot_tilde, g_tilde, rr, nn, ~, g_0, q] = mm10_DJGM_GH(props);
                
        g_s = zeros(props.nslip,1);
        h = zeros(props.nslip,1);
        dh = h;
        
        % Equation [6], g_s equation has to be modified
        for slip_b = 1:props.nslip
            temp = g_tilde(slip_b)*(gamma_dot(slip_b)/gamma_dot_tilde(slip_b))^nn(slip_b);
            if temp > gtol*g_0(slip_b)
                g_s(slip_b) = temp; % threshold to ensure no slip during elastic response
                h(slip_b) = h_0(slip_b)*abs(1-g(slip_b)/g_s(slip_b))^rr(slip_b)*sign(1-g(slip_b)/g_s(slip_b));
                dg_s = g_tilde(slip_b)*nn(slip_b)/gamma_dot(slip_b)*...
                    (gamma_dot(slip_b)/gamma_dot_tilde(slip_b))^nn(slip_b)*...
                    dslip(slip_b);
                dh(slip_b) = h_0(slip_b)*rr(slip_b)*abs(1-g(slip_b)/g_s(slip_b))^(rr(slip_b)-1)...
                    *(dg_s*g(slip_b)/g_s(slip_b)^2);
            else
                g_s(slip_b) = 2/3*g_0(slip_b); % threshold to ensure no slip during elastic response
                h(slip_b) = h_0(slip_b)*abs(1-g(slip_b)/g_s(slip_b))^rr(slip_b)*sign(1-g(slip_b)/g_s(slip_b));
                dh(slip_b) = 0;
            end
        end
        
        % Equation [5]
        dg_dot = zeros(props.nslip,6);
        for slip_a = 1:props.nslip
            for slip_b = 1:props.nslip
                ms = np1.ms(1:6,slip_b);
                dtdstress = ms;
                dg_dot(slip_a,1:6) = dg_dot(slip_a,1:6) + q(slip_a, slip_b)*...
                    (dh(slip_b)*abs(gamma_dot(slip_b)) + ...
                     h(slip_b)*sign(gamma_dot(slip_b))*dslip(slip_b))*dtdstress';
            end

          et(slip_a,1:6) = dt*dg_dot(slip_a,1:6); % d(18)/dstress
        end

% 
% c
%       return
end