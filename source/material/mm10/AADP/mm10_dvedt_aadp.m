function dvdt = mm10_dvedt_aadp(props,np1,n,tau,rho,vbar_e,alpha,dt)
% Function to compute derivative of vbar_e with respect to resolved shear
% stress on system alpha

% Set some properties
dFe = props.G_0_v;
k = props.boltzman;
theta = np1.temp;
s_pe = props.tau_hat_v;
p_e = props.p_v;
q_e = props.q_v;

% Compute one dependency
s_de = mm10_sed_aadp(props, np1, n, rho, alpha); % (A.13)

% Long chain rule
if (abs(tau)/(s_pe + s_de)) > 1
v_oe = props.c1;
dvdt = v_oe * (-q_e*(dFe/k/theta)*(1.d0 - 1.d0)^(q_e-1.d0)) ...
       * (- p_e*(1.d0)^(p_e-1.d0)) * sign(tau)/(s_pe + s_de);
else
dvdt = vbar_e * (-q_e*(dFe/k/theta)*(1.d0 - (abs(tau)/(s_pe + s_de))^p_e)^(q_e-1.d0)) ...
       * (- p_e*(abs(tau)/(s_pe + s_de))^(p_e-1.d0)) * sign(tau)/(s_pe + s_de);
end