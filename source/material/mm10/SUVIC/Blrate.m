function Bldot = Blrate(eps,Bl,props)

A1l = props.A1l;
Blp = Blprime(eps,props);
A2l = props.A2l;
C = props.C;
q = props.q;
Bldot = A1l*eps - (A1l/Blp)*abs(eps)*Bl - A2l*maccauley((abs(Bl)-Blp)/C)^q*sign(Bl);