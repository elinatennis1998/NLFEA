% exponential term in Roters model
function vbar = mm10_vbar_mrr(props, np1, n, tau, rho, alpha,dt)

% Set some properties
Qslip = props.Qslip;
k = props.boltzman;
theta = np1.temp;
p_e = props.p_v;
q_e = props.q_v;

% Compute one dependency
gamma_0 = mm10_gamma_0_aadp(props, np1, n, rho, alpha); % (15)
tpass = mm10_tpass_aadp(props, np1, n, rho, alpha); % (16)
tcut = mm10_tcut_aadp(props, np1, n, rho, alpha); % (17)

% Evaluate the equation
fract = ((abs(tau)-tpass/tcut));
if fract > 1
%     b = v_oe;
%     x = abs(tau)/(s_pe + s_de);
% %     m = mm10_dvedt_aadp(props,np1,n,s_pe+s_de,tt,b,alpha,dt);
%     m = b * (-q_e*(dFe/k/theta)*(1.d0 - 1.d0)^(q_e-1.d0)) ...
%        * (- p_e*(1.d0)^(p_e-1.d0)) * sign(tau)/(s_pe + s_de);
%     vbar_e = m*x + b;
else
vbar = gamma_0 * exp (-(Qslip/k/theta)*(1.d0 - fract^p_e)^q_e); %(14)
end