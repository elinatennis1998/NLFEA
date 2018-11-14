% Compute average edge mobile dislocation velocity for AADP model
function vbar_e = mm10_vbare_aadp(props, np1, n, tau, rho, alpha,dt)

% Set some properties
v_oe = props.c1;
dFe = props.G_0_v;
k = props.boltzman;
theta = np1.temp;
s_pe = props.tau_hat_v;
p_e = props.p_v;
q_e = props.q_v;

% Compute one dependency
s_de = mm10_sed_aadp(props, np1, n, rho, alpha); % (A.13)

% Evaluate the equation
if (abs(tau)/(s_pe + s_de)) > 1
    b = v_oe;
    x = abs(tau)/(s_pe + s_de);
%     m = mm10_dvedt_aadp(props,np1,n,s_pe+s_de,tt,b,alpha,dt);
    m = b * (-q_e*(dFe/k/theta)*(1.d0 - 1.d0)^(q_e-1.d0)) ...
       * (- p_e*(1.d0)^(p_e-1.d0)) * sign(tau)/(s_pe + s_de);
    vbar_e = m*x + b;
else
vbar_e = v_oe * exp (-(dFe/k/theta)*(1.d0 - (abs(tau)/(s_pe + s_de))^p_e)^q_e); %(A.11)
end