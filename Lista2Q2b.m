clc
clear all

%Considerando Df = 0

A = [-0.9 1; 0 0.9];
B1 = [-2; 4];
C1 = [1 0];
C2 = [-1 -3];
D11 = [0.1]; 
D21 = [0];


Z = sdpvar(2);
X = sdpvar(2);
F = sdpvar(1,2);
L = sdpvar(2,1);
G = sdpvar(2);
M = sdpvar(1);

N11 = [Z Z Z*B1];
N12 = [Z X X*B1+L*D21];
N13 = [B1'*Z B1'*X+D21'*L' M];
des1 = [N11; N12; N13];

N21 = [Z Z A'*Z A'*X+C2'*L'+G' C1'-F'];
N22 = [Z' X A'*Z A'*X+C2'*L' C1'];
N23 = [(A'*Z)' (A'*Z)' Z Z zeros(2,1)];
N24 = [(A'*X+C2'*L'+G')' (A'*X+C2'*L')' Z' X zeros(2,1)];
N25 = [(C1'-F')' (C1')' (zeros(2,1))' (zeros(2,1))' eye(1)];
des2 = [N21; N22; N23; N24; N25];


obj = trace(M);
prob = [des1>=0; des2>=0];
S = sdpsettings('solver', 'lmilab');
diagnostics = optimize(prob,obj,S);
Z = double(Z);
X = double(X);
F = double(F);
L = double(L);
G = double(G);
M = double(M);

ro = (trace(M))^0.5;

U = eye(2);
V = eye(2) - X*inv(Z);

Af = inv(U')*G*inv(V*Z);
Bf = inv(U')*L;
Cf = F*inv(V*Z);
Df = 0;

x{1} = [10 -10]';
xf{1} = [0 0]';
w(1) = 1;
z(1) = C1*x{1}+0.1*w(1);
zf(1) = Cf*xf{1}+Df*w(1);
e = (z(1)-zf(1))^2;
for i = 2:50
    w(i) = exp(-0.25*i)*cos(10*i);
    x{i} = A*x{i-1} + B1*w(i-1);
    xf{i} = Af*x{i-1} + Bf*w(i-1);
    z(i) = C1*x{i} + D11*w(i);
    zf(i) = Cf*xf{i}+Df*w(i); 
    e = e + (z(i)-zf(i))^2;
end

e = e/50;

figure(1)
hold on
stem(z, 'b')
stem(zf, 'g')
stem(w, 'r')
legend('z','zf','w')


