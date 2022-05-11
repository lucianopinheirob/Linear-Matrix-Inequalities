clc
clear all


P = sdpvar(4,4,'symmetric');
H = sdpvar(4,4,'symmetric');
Z = sdpvar(1);
X = sdpvar(4,4,'full');
Y = sdpvar(4,4,'full');
L = sdpvar(2,4);
F = sdpvar(4,2);
Q = sdpvar(4,4,'full');
R = sdpvar(2,2,'full');
S = sdpvar(4,4,'full');
J = sdpvar(4,4,'full');
u = sdpvar(1);

A = [0 1 0 0; -1 1 0 1; -1 -1 -1 1; 1 0 -1 0];
B1 = [0 1 0; 1 -1 0; 0 -1 -1; 1 1 0];
B2 = [2 -2; -2 0; -2 0; 0 2];
C1 = [0 -1 2 1];
C2 = [-1 0 -1 0; 1 1 2 -2];
D11 = [0 0 0];
D12 = [0 1];
D21 = [0 1 1; 1 1 -1];


des1 = u-trace(Z);

M1 = [P J A*X+B2*L A+B2*R*C2 B1+B2*R*D21];
M2 = [(P)' H Q Y*A+F*C2 Y*B1+F*D21];
M3 = [(J)' (H)' X+X'-P eye(4)+S'-J zeros(4,3)];
M4 = [(A*X+B2*L)' (Q)' (X+X'-P)' Y+Y'-H zeros(4,3)];
M5 = [(B1+B2*R*D21)' (Y*B1+F*D21)' (zeros(4,3))' (zeros(4,3))' eye(3)];
des2 = [M1; M2; M3; M4; M5];

N1 = [Z C1*X+D12*L C1+D12*R*C2];
N2 = [(C1*X+D12*L)' X+X'-P eye(4)+S'-J];
N3 = [(C1*X+D12*L)' (X+X'-P)' Y+Y'-H];
des3 = [N1; N2; N3];

des4 = D11+D12*R*D21;

obj = [u];
prob = [des1>=0; des2>=0; des3>=0 des4==0];
Sett = sdpsettings('solver', 'sedumi');
solvesdp(prob,obj,Sett);

P = double(P);
H = double(H);
Z = double(Z);
X = double(X);
Y = double(Y);
L = double(L);
F = double(F);
R = double(R);
Q = double(Q);
R = double(R);
S = double(S);
J = double(J);
u = double(u);


U = eye(4);
V = S-Y*X;

contr = [inv(V) -inv(V)*Y*B2; zeros(2,4) eye(2)]*[Q-Y*A*X F; L R]*[inv(U) zeros(4,2); -C2*X*inv(U) eye(2)]; 

Ac = contr(1:4, 1:4)
Bc = contr(1:4, 5:6)
Cc = contr(5:6, 1:4)
Dc = contr(5:6, 5:6)
ro = u^0.5



