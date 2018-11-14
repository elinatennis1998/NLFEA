function Kp = Kprime(eps,props)

Blp = Blprime(eps,props);
Bsp = Bsprime(eps,props);
Rp = Rprime(eps,props);
A = props.R0;
N = props.N;
spds0 = sigmaprime(eps,props);
s0 = props.s0;
sp = spds0*s0;
Kp = (sp - (Bsp + Blp + Rp))/(eps/A)^(1/N);