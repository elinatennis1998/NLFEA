function Bsp = Bsprime(eps,props)

B0s = props.B0s;
% eps0 = props.eps0;
m = props.m;
spds0 = sigmaprime(eps,props);
Bsp = B0s*(spds0)^(m);