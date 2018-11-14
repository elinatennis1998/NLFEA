function dvdt = mm10_dvsdt_aadp(props,np1,n,tau,rho,vbar_s,alpha,dt)
% Function to compute derivative of vbar_s with respect to resolved shear
% stress on system alpha

% Set some properties
dFs = props.G_0_y;
k = props.boltzman;
theta = np1.temp;
s_ps = props.tau_hat_y;
p_s = props.p_y;
q_s = props.q_y;

% Compute one dependency
s_ds = mm10_ssd_aadp(props, np1, n, rho, alpha); % (A.14)

% Long chain rule
if (abs(tau)/(s_ps + s_ds)) > 1
v_os = props.c3;
dvdt = v_os * (-q_s*(dFs/k/theta)*(1.d0 - 1.d0)^(q_s-1.d0)) ...
       * (- p_s*(1.d0)^(p_s-1.d0)) * sign(tau)/(s_ps + s_ds);
else
dvdt = vbar_s * (-q_s*(dFs/k/theta)*(1.d0 - (abs(tau)/(s_ps + s_ds))^p_s)^(q_s-1.d0)) ...
       * (- p_s*(abs(tau)/(s_ps + s_ds))^(p_s-1.d0)) * sign(tau)/(s_ps + s_ds);
end