clc
clear all

A = [1 1; 5 -2];
B = [3; 2];
C = [1 0];

z11 = sdpvar(1);

w11 = sdpvar(1);
w22 = sdpvar(1);


Z = [z11 0];
W = [w11 0; 0 w22];

des1 = A*W+W*A'+B*Z+Z'*B';
des2 = W;

obj = [];
prob = [des1<=0; des2>=0];
S = sdpsettings('solver', 'lmilab');
diagnostics = optimize(prob,obj,S);

Z = double(Z);
W = double(W);

K = Z*inv(W);

L = K(1)

