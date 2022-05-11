clc
clear all

A = [1 -2 6; 7 2 -5; -9 7 4];
B = [2 0; -3 -1; 2 1];

z11 = sdpvar(1);
z22 = sdpvar(1);
z13 = sdpvar(1);

w11 = sdpvar(1);
w13 = sdpvar(1);
w22 = sdpvar(1);
w31 = sdpvar(1);
w33 = sdpvar(1);

Z = [z11 0 z13; 0 z22 0];
W = [w11 0 w13; 0 w22 0; w31 0 w33];

des1 = A*W+W*A'+B*Z+Z'*B';
des2 = W;

obj = [];
prob = [des1<=0; des2>=0];
S = sdpsettings('solver', 'lmilab');
diagnostics = optimize(prob,obj,S);

Z = double(Z);
W = double(W);

K = Z*inv(W)


