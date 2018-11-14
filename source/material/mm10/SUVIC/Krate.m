function Kdot = Krate(eps,K,props)

A5 = props.A5;
Kp = Kprime(eps,props);
A6 = props.A6;
C = props.C;
u = props.u;
Kdot = A5*(Kp-K)/Kp*abs(eps) - A6*maccauley((abs(K)-Kp)/C)^u;