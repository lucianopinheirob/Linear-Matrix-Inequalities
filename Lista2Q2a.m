clc
clear all

A = [-0.9 1; 0 0.9];
B1 = [-2; 4];
C1 = [1 0];
C2 = [-1 -3];
D11 = [0.1]; 
D21 = [0];


Z = sdpvar(2,2,'symmetric');
X = sdpvar(2,2,'symmetric');
F = sdpvar(1,2);
L = sdpvar(2,1);
G = sdpvar(2,2,'full');
Df = sdpvar(1);
u = sdpvar(1);

M1 = [Z Z Z*A Z*A Z*B1 zeros(2,1)];
M2 = [Z' X X*A+L*C2+G X*A+L*C2 X*B1+L*D21 zeros(2,1)];
M3 = [(Z*A)' (X*A+L*C2+G)' Z Z zeros(2,1) C1'-C2'*Df'-F'];
M4 = [(Z*A)' (X*A+L*C2)' Z' X zeros(2,1) C1'-C2'*Df'];
M5 = [(Z*B1)' (X*B1+L*D21)' (zeros(2,1))' (zeros(2,1))' eye(1) D11'-D21'*Df']
M6 = [(zeros(2,1))' (zeros(2,1))' (C1'-C2'*Df'-F')' (C1'-C2'*Df')' (D11'-D21'*Df')' u*eye(1)]

des1 = [M1; M2; M3; M4; M5; M6];

obj = [u];
prob = [des1>=0];
S = sdpsettings('solver', 'sedumi');
diagnostics = optimize(prob,obj,S);
Z = double(Z);
X = double(X);
F = double(F);
L = double(L);
G = double(G);
Df = double(Df);
u = double(u);

gama = u^0.5;

U = eye(2);
V = eye(2) - X*inv(Z);

Af = inv(U')*G*inv(V*Z);
Bf = inv(U')*L;
Cf = F*inv(V*Z);
Df;

x{1} = [10 -10]';
xf{1} = [0 0]';
w(1) = 1;
z(1) = C1*x{1}+0.1*w(1);
zf(1) = Cf*x{1}+Df*w(1);
e = (z(1)-zf(1))^2;
for i = 2:50
    w(i) = exp(-0.25*i)*cos(10*i);
    x{i} = A*x{i-1} + B1*w(i-1);
    z(i) = C1*x{i} + D11*w(i);
    xf{i} = Af*x{i-1} + Bf*w(i-1);
    zf(i) = Cf*x{i}+Df*w(i); 
    e = e + (z(i)-zf(i))^2;
end

e = e/50;

figure(1)
hold on
stem(z, 'b')
stem(zf, 'g')
stem(w, 'r')
legend('z','zf','w')

Af
Bf
Cf
Df
gama
e



