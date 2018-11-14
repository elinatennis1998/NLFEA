% c
% c
% c     AADP:
% c
% c           Setup AADP hardening
function [props, np1, n] = mm10_setup_aadp(props, np1, n)
%       use mm10_defs
% init_hard = ?;
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
% c     
%       double precision :: dgc, init_hard
% init_hard=0.1;
% % c
% % c     Continuum effective rate
%       dgc = np1.dg / np1.tinc;
% % c
% % c     New shear modulus
%       np1.mu_harden = props.mu_0 - props.D_0 / (exp(props.T_0/...
%           np1.temp)- 1.0);
% % c
% % c     Get the threshold contributions
%       np1.tau_v = props.tau_hat_v*(1-(props.boltzman*np1.temp...
%           /(np1.mu_harden*(props.burgers^3)*props.G_0_v)*...
%           log(props.eps_dot_0_v/dgc))^(1/props.q_v))...
%           ^(1/props.p_v);
%       np1.tau_y = props.tau_hat_y*(1-(props.boltzman*np1.temp...
%           /(np1.mu_harden*(props.burgers^3)*props.G_0_y)*...
%           log(props.eps_dot_0_y/dgc))^(1/props.q_y))...
%           ^(1/props.p_y);
% % c
% % c     Used existing labels as a convenience, actually get/set the history
%       np1.u(1) = np1.tau_y;
%       if (n.u(1) < 0.0)% then
%         n.tau_y = np1.tau_y;
%       else
%         n.tau_y = n.u(1);
%       end
% 
%       np1.u(2) = np1.mu_harden;
%       if (n.u(2) < 0.0)% then
%         n.mu_harden = np1.mu_harden;
%       else
%         n.mu_harden = n.u(2);
%       end
% % c
% % c     Same here -- check for previous step flag
%       if (n.tau_tilde < 0.0)% then
%         n.tau_tilde = props.tau_a + (np1.mu_harden/props.mu_0)*np1.tau_y...
%             + init_hard;
%       end
% % c     
% %       return
end