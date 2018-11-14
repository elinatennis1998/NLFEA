function eps = sliprate(s,Bs,Bl,R,K,props)

N = props.N;
A = props.A;
B = Bs + Bl;
eps = A*maccauley((abs(s-B)-R)/K).^N.*sign(s-B);