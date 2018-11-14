% Compute edge SSD density for slip system i fro AADP model
function gamma_0 = mm10_gamma_0_halite(props, np1, n, rho, alpha)

% Set some properties
k = props.boltzman;
theta = np1.temp;
G = props.mu_0;
b = props.burgers;
if props.s_type == 6
c1 = props.cp_001*ones(12,1);
c3 = props.cp_007*ones(12,1);
elseif props.s_type == 11
c1 = [props.cp_001*ones(6,1);props.cp_002*ones(6,1);props.cp_003*ones(12,1)];
c3 = [props.cp_007*ones(6,1);props.cp_008*ones(6,1);props.cp_009*ones(12,1)];
else
    error(' ');
end

% Get parallel dislocations
[~,rhoP] = mm10_rhoFP_halite(props, np1, n, rho, alpha);

% Evaluate the equation
gamma_0 = k*theta./(c1(alpha).*c3(alpha)*G*b^2).*sqrt(rhoP);