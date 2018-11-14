% Compute edge SSD density for slip system i fro AADP model
function gamma_0 = mm10_gamma_0_ornl(props, np1, n, rho, alpha)

% Set some properties
k = props.boltzman;
theta = np1.temp;
G = props.mu_0;
b = props.burgers;
c1 = props.c1;
c3 = props.c3;

% Get parallel dislocations
[~,rhoP] = mm10_rhoFP_ornl(props, np1, n, rho, alpha);

% Evaluate the equation
gamma_0 = k*theta/(c1*c3*G*b^2)*sqrt(rhoP);