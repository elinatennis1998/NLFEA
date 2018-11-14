% Compute average screw mobile dislocation velocity for AADP model
function vbar_s = mm10_vbars_aadp(props, np1, n, tau, rho, alpha, dt)

% Set some properties
v_os = props.c3;
dFs = props.G_0_y;
k = props.boltzman;
theta = np1.temp;
s_ps = props.tau_hat_y;
p_s = props.p_y;
q_s = props.q_y;

% Compute one dependency
s_ds = mm10_ssd_aadp(props, np1, n, rho, alpha); % (A.14)

% Evaluate the equation
if (abs(tau)/(s_ps + s_ds)) > 1
    b = v_os;
    x = abs(tau)/(s_ps + s_ds);
%     m = mm10_dvsdt_aadp(props,np1,n,s_ps+s_ds,tt,b,alpha,dt);
    m = b * (-q_s*(dFs/k/theta)*(1.d0 - 1.d0)^(q_s-1.d0)) ...
       * (- p_s*(1.d0)^(p_s-1.d0)) * sign(tau)/(s_ps + s_ds);
    vbar_s = m*x + b;
else
vbar_s = v_os * exp (-(dFs/k/theta)*(1.d0 - (abs(tau)/(s_ps + s_ds))^p_s)^q_s); %(A.12)
end