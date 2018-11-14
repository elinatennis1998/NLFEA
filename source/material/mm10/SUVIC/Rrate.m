function Rdot = Rrate(eps,R,props)

A3 = props.A3;
Rp = Rprime(eps,props);
A4 = props.A4;
C = props.C;
p = props.p;
Rdot = A3*(Rp-R)/Rp*abs(eps) - A4*maccauley((abs(R)-Rp)/C)^p;