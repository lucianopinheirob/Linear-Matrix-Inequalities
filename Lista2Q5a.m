clc
clear all

A = [1 -2 6; 7 2 -5; -9 7 4];
B = [2 0; -3 -1; 2 1];

Z = sdpvar(2,3);
W = sdpvar(3);


des1 = A*W+W*A'+B*Z+Z'*B';
des2 = W;

prob = [des1<=0; des2>=0];
S = sdpsettings('solver', 'sedumi');
diagnostics = optimize(prob,obj,S);

Z = double(Z);
W = double(W);

K = Z*inv(W)


