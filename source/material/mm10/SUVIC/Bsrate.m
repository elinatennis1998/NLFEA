function Bsdot = Bsrate(eps,Bs,props)

A1s = props.A1s;
Bsp = Bsprime(eps,props);
A2s = props.A2s;
C = props.C;
q = props.q;
Bsdot = A1s*eps - (A1s/Bsp)*abs(eps)*Bs - A2s*maccauley((abs(Bs)-Bsp)/C)^q*sign(Bs);