clc
clear all

A = [-4 -7 2; 0 -4 -2; -3 2 0];
B = [1; 0; 1];
C = [1 0 0];

z11 = sdpvar(1);

w11 = sdpvar(1);
w22 = sdpvar(1);
w33 = sdpvar(1);
w23 = sdpvar(1);

Z = [z11 0 0];
W = [w11 0 0; 0 w22 w23; 0 w23 w33];

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

