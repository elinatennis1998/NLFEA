function spds0 = sigmaprime(eps,props)

eps0 = props.eps0;
% s0 = props.s0;
% m = props.m;

% spds0 = max(log(eps/eps0),0*eps);
spds0 = log(eps/eps0);